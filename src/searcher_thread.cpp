// Copyright 2012 Evrytania LLC (http://www.evrytania.com)
//
// Written by James Peroulas <james@evrytania.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>
#include <itpp/itbase.h>
#include <itpp/signal/transforms.h>
#include <itpp/stat/misc_stat.h>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/thread.hpp>
#include <boost/thread/condition.hpp>
#include <list>
#include <sstream>
#include <signal.h>
#include <queue>
//#include <valgrind/callgrind.h>
#include <sys/syscall.h>
#include <sys/types.h>
#include <curses.h>
#include "rtl-sdr.h"
#include "common.h"
#include "macros.h"
#include "lte_lib.h"
#include "constants.h"
#include "capbuf.h"
#include "itpp_ext.h"
#include "searcher.h"
#include "dsp.h"
#include "LTE-Tracker.h"
#include "filter_coef.h"

#ifdef HAVE_HACKRF
#include "hackrf.h"
#endif

#ifdef HAVE_BLADERF
#include <libbladeRF.h>
#endif // HAVE_BLADERF

using namespace itpp;
using namespace std;

#define DS_COMB_ARM 2

// This is the searcher process. It requests captured data from the main
// thread and launches a new thread for every cell it finds. Each new
// cell thread then requests sample data from the main thread.
void searcher_thread(
  capbuf_sync_t & capbuf_sync,
  global_thread_data_t & global_thread_data,
  tracked_cell_list_t & tracked_cell_list
) {
  if (verbosity>=1) {
    cout << "Searcher process has been launched." << endl;
  }

  global_thread_data.searcher_thread_id=syscall(SYS_gettid);

  if (nice(20)==-1) {
    cerr << "Error: could not reduce searcher process priority" << endl;
    ABORT(-1);
  }

  // Keep track of serial numbers to be used when launching a new
  // tracker thread.
  ivec serial_num(504);
  serial_num=1;

  // Shortcut
  const double & fc_requested=global_thread_data.fc_requested;
  const double & fc_programmed=global_thread_data.fc_programmed;
  const double & fs_programmed=global_thread_data.fs_programmed;
//  double freq_correction;

  const bool sampling_carrier_twist = global_thread_data.sampling_carrier_twist();
  double k_factor = global_thread_data.k_factor();
//  double correction = global_thread_data.correction();

  vec coef(( sizeof( chn_6RB_filter_coef )/sizeof(float) ));
  for (uint16 i=0; i<length(coef); i++) {
    coef(i) = chn_6RB_filter_coef[i];
  }

  // Calculate the threshold vector
  const uint8 thresh1_n_nines=12;
  double rx_cutoff=(6*12*15e3/2+4*15e3)/(FS_LTE/16/2);

  // for PSS correlate
  //cout << "DS_COMB_ARM override!!!" << endl;
#define DS_COMB_ARM 2
  mat xc_incoherent_collapsed_pow;
  imat xc_incoherent_collapsed_frq;
  vector <mat> xc_incoherent_single(3);
  vector <mat> xc_incoherent(3);
  vec sp_incoherent;
  vector <mat> xc(3);
  vec sp;

  // for SSS detection
#define THRESH2_N_SIGMA 3
  vec sss_h1_np_est_meas;
  vec sss_h2_np_est_meas;
  cvec sss_h1_nrm_est_meas;
  cvec sss_h2_nrm_est_meas;
  cvec sss_h1_ext_est_meas;
  cvec sss_h2_ext_est_meas;
  mat log_lik_nrm;
  mat log_lik_ext;

  // for time frequency grid
  // Extract time and frequency grid
  cmat tfg;
  vec tfg_timestamp;
  // Compensate for time and frequency offsets
  cmat tfg_comp;
  vec tfg_comp_timestamp;

  vec period_ppm;
  double xcorr_pss_time;

  Real_Timer tt; // for profiling

  // Get the current frequency offset (because it won't change anymore after main thread launches this thread, so move it outside loop)
  vec f_search_set(1);
  cmat pss_fo_set;// pre-generate frequencies offseted pss time domain sequence
    // because it is already included in global_thread_data.frequency_offset();
  f_search_set(0)=global_thread_data.frequency_offset();
  pss_fo_set_gen(f_search_set, pss_fo_set);

//  freq_correction = fc_programmed*(correction-1)/correction;

//  cout << opencl_platform << " " << opencl_device <<  " " << sampling_carrier_twist << "\n";
  uint16 opencl_platform = global_thread_data.opencl_platform();
  uint16 opencl_device = global_thread_data.opencl_device();
  lte_opencl_t lte_ocl(opencl_platform, opencl_device);

  #ifdef USE_OPENCL
  uint16 filter_workitem = global_thread_data.filter_workitem();
  uint16 xcorr_workitem = global_thread_data.xcorr_workitem();
  #ifdef FILTER_MCHN_SIMPLE_KERNEL
  lte_ocl.setup_filter_mchn((string)"filter_mchn_simple_kernel.cl", CAPLENGTH, length(f_search_set)*3, pss_fo_set.cols(), xcorr_workitem);
  #else
  lte_ocl.setup_filter_mchn((string)"filter_mchn_kernels.cl", CAPLENGTH, length(f_search_set)*3, pss_fo_set.cols(), xcorr_workitem);
  #endif
  lte_ocl.setup_filter_my((string)"filter_my_kernels.cl", CAPLENGTH, filter_workitem);
  #endif

  // Loop forever.
  while (true) {
    // Used to measure searcher cycle time.
    tt.tic();

    // Request data.
    {
      boost::mutex::scoped_lock lock(capbuf_sync.mutex);
      capbuf_sync.request=true;

      // Wait for data to become ready.
      capbuf_sync.condition.wait(lock);
    }

    // Results are stored in this vector.
    list<Cell> detected_cells;

    // Local reference to the capture buffer.
    cvec &capbuf=capbuf_sync.capbuf;

    capbuf = capbuf - mean(capbuf); // remove DC

    #ifdef USE_OPENCL
      lte_ocl.filter_my(capbuf); // be careful! capbuf.zeros() will slow down the xcorr part pretty much!
    #else
      filter_my(coef, capbuf);
    #endif

    // Correlate
    uint16 n_comb_xc;
    uint16 n_comb_sp;
    if (verbosity>=2) {
      cout << "  Calculating PSS correlations" << endl;
    }

    sampling_ppm_f_search_set_by_pss(lte_ocl, 0, capbuf, pss_fo_set, 1, 0, f_search_set, period_ppm, xc, xcorr_pss_time);

    xcorr_pss(capbuf,f_search_set,DS_COMB_ARM,fc_requested,fc_programmed,fs_programmed,xc,xc_incoherent_collapsed_pow,xc_incoherent_collapsed_frq,xc_incoherent_single,xc_incoherent,sp_incoherent,sp,n_comb_xc,n_comb_sp,sampling_carrier_twist,k_factor);

    // Calculate the threshold vector
    double R_th1=chi2cdf_inv(1-pow(10.0,-thresh1_n_nines),2*n_comb_xc*(2*DS_COMB_ARM+1));
    vec Z_th1=R_th1*sp_incoherent/rx_cutoff/137/n_comb_xc/(2*DS_COMB_ARM+1);

    // Search for the peaks
    if (verbosity>=2) {
      cout << "  Searching for and examining correlation peaks..." << endl;
    }
    peak_search(xc_incoherent_collapsed_pow,xc_incoherent_collapsed_frq,Z_th1,f_search_set,fc_requested,fc_programmed,xc_incoherent_single,DS_COMB_ARM, sampling_carrier_twist, (const double)k_factor, detected_cells);

    // Loop and check each peak
    list<Cell>::iterator iterator=detected_cells.begin();
    int tdd_flag = 1;
    while (iterator!=detected_cells.end()) {
      tdd_flag = !tdd_flag;

      // Detect SSS if possible
#define THRESH2_N_SIGMA 3
      (*iterator)=sss_detect((*iterator),capbuf,THRESH2_N_SIGMA,fc_requested,fc_programmed,fs_programmed,sss_h1_np_est_meas,sss_h2_np_est_meas,sss_h1_nrm_est_meas,sss_h2_nrm_est_meas,sss_h1_ext_est_meas,sss_h2_ext_est_meas,log_lik_nrm,log_lik_ext,sampling_carrier_twist,tdd_flag);
      if ((*iterator).n_id_1!=-1) {
        if (verbosity>=2) {
          cout << "Detected PSS/SSS correspoding to cell ID: " << (*iterator).n_id_cell() << endl;
        }

        // Check to see if this cell has already been detected previously.
        bool match=false;
        {
          boost::mutex::scoped_lock lock(tracked_cell_list.mutex);
          list<tracked_cell_t *>::iterator tci=tracked_cell_list.tracked_cells.begin();
          match=false;
          while (tci!=tracked_cell_list.tracked_cells.end()) {
            if ((*(*tci)).n_id_cell==(*iterator).n_id_cell()) {
              match=true;
              break;
            }
            ++tci;
          }
        }
        if (match) {
          if (verbosity>=2) {
            cout << "Cell already being tracked..." << endl;
          }
          ++iterator;
          continue;
        }

        // Fine FOE
        (*iterator)=pss_sss_foe((*iterator),capbuf,fc_requested,fc_programmed,fs_programmed,sampling_carrier_twist,tdd_flag);

        // Extract time and frequency grid
        extract_tfg((*iterator),capbuf,fc_requested,fc_programmed,fs_programmed,tfg,tfg_timestamp,sampling_carrier_twist);

        // Create object containing all RS
        RS_DL rs_dl((*iterator).n_id_cell(),6,(*iterator).cp_type);

        // Compensate for time and frequency offsets
        (*iterator)=tfoec((*iterator),tfg,tfg_timestamp,fc_requested,fc_programmed,rs_dl,tfg_comp,tfg_comp_timestamp,sampling_carrier_twist);

        // Finally, attempt to decode the MIB
        (*iterator)=decode_mib((*iterator),tfg_comp,rs_dl);
        //cout << (*iterator) << endl << endl;
        //sleep(100000);
        if ((*iterator).n_rb_dl==-1) {
          // No MIB could be successfully decoded.
          iterator=detected_cells.erase(iterator);
          continue;
        }

        /*
        if (verbosity>=1) {
          cout << "Detected a new cell!" << endl;
          cout << "  cell ID: " << (*iterator).n_id_cell() << endl;
          cout << "  RX power level: " << db10((*iterator).pss_pow) << " dB" << endl;
          cout << "  residual frequency offset: " << (*iterator).freq_superfine << " Hz" << endl;
          cout << "  frame start: " << (*iterator).frame_start << endl;
        }
        */

  //      cout << ((*iterator).frame_start) <<  " " << k_factor << " " << capbuf_sync.late << (*iterator).k_factor <<  "\n";
  //      (*iterator).frame_start = 3619.95;
  //      capbuf_sync.late = 0;
        // Launch a cell tracker process!
        //tracked_cell_t * new_cell = new tracked_cell_t((*iterator).n_id_cell(),(*iterator).n_ports,(*iterator).cp_type,(*iterator).frame_start/k_factor+capbuf_sync.late,serial_num((*iterator).n_id_cell()));
        tracked_cell_t * new_cell = new tracked_cell_t(
          (*iterator).n_id_cell(),
          (*iterator).n_ports,
          (*iterator).duplex_mode,
          (*iterator).cp_type,
          (*iterator).n_rb_dl,
          (*iterator).phich_duration,
          (*iterator).phich_resource,
          (*iterator).frame_start*(FS_LTE/16)/(fs_programmed*k_factor)+capbuf_sync.late,
          serial_num((*iterator).n_id_cell())//,
  //        (*iterator).freq_superfine
        );

        serial_num((*iterator).n_id_cell())++;
        // Cannot launch thread here. If thread was launched here, it would
        // have the same (low) priority as the searcher thread.
        {
          boost::mutex::scoped_lock lock(tracked_cell_list.mutex);
          tracked_cell_list.tracked_cells.push_back(new_cell);
        }
        //CALLGRIND_START_INSTRUMENTATION;
  #define MAX_DETECTED 1e6
        static uint32 n_found=0;
        n_found++;
        if (n_found==MAX_DETECTED) {
          cout << "Searcher thread has stopped!" << endl;
          sleep(1000000);
        }

        ++iterator;

      } else {
        // No SSS detected.
        iterator=detected_cells.erase(iterator);
        continue;
      }
    }

    global_thread_data.searcher_cycle_time(tt.toc());
//    global_thread_data.searcher_cycle_time(xcorr_pss_time);
  }
  // Will never reach here...
}

