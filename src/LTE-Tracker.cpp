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
#include <itpp/itbase.h>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/thread.hpp>
#include <boost/thread/condition.hpp>
#include <list>
#include <sstream>
#include <signal.h>
#include <queue>
#include "rtl-sdr.h"
#include "common.h"
#include "macros.h"
#include "lte_lib.h"
#include "constants.h"
#include "capbuf.h"
#include "itpp_ext.h"
#include "searcher.h"
#include "dsp.h"
#include "rtl-sdr.h"

using namespace itpp;
using namespace std;

uint8 verbosity=1;
// Declared as global so the sig handler can have access to it.
rtlsdr_dev_t * dev=NULL;

static void sighandler(
  int signum
) {
  cerr << "Error: caught signal, exiting!" << endl;
  if (dev!=NULL) {
    rtlsdr_close(dev);
  }
  exit(-1);
}

// Simple usage screen.
void print_usage() {
  cout << "LTE-Tracker v" << MAJOR_VERSION << "." << MINOR_VERSION << "." << PATCH_LEVEL << " (" << BUILD_TYPE << ") help screen" << endl << endl;
  cout << "LTE-Tracker -f frequency [optional_parameters]" << endl;
  cout << "  Basic options" << endl;
  cout << "    -h --help" << endl;
  cout << "      print this help screen" << endl;
  cout << "    -v --verbose" << endl;
  cout << "      increase status messages from program" << endl;
  cout << "    -b --brief" << endl;
  cout << "      reduce status messages from program" << endl;
  cout << "    -i --device-index N" << endl;
  cout << "      specify which attached RTLSDR dongle to use" << endl;
  cout << "  Frequency search options:" << endl;
  cout << "    -f --freq fc" << endl;
  cout << "      frequency where cells are located" << endl;
  cout << "  Dongle LO correction options:" << endl;
  cout << "    -p --ppm ppm" << endl;
  cout << "      crystal remaining PPM error" << endl;
  cout << "    -c --correction c" << endl;
  cout << "      crystal correction factor" << endl;
  cout << "See CellSearch for documentation on c and ppm." << endl;
}

// Parse the command line arguments and return optional parameters as
// variables.
// Also performs some basic sanity checks on the parameters.
void parse_commandline(
  // Inputs
  const int & argc,
  char * const argv[],
  // Outputs
  double & fc,
  double & ppm,
  double & correction,
  int & device_index
) {
  // Default values
  fc=-1;
  ppm=100;
  correction=1;
  device_index=-1;

  while (1) {
    static struct option long_options[] = {
      {"help",         no_argument,       0, 'h'},
      {"verbose",      no_argument,       0, 'v'},
      {"brief",        no_argument,       0, 'b'},
      {"freq",         required_argument, 0, 'f'},
      {"ppm",          required_argument, 0, 'p'},
      {"correction",   required_argument, 0, 'c'},
      {"device-index", required_argument, 0, 'i'},
      {0, 0, 0, 0}
    };
    /* getopt_long stores the option index here. */
    int option_index = 0;
    int c = getopt_long (argc, argv, "hvbf:p:c:i:",
                     long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c) {
      char * endp;
      case 0:
        // Code should only get here if a long option was given a non-null
        // flag value.
        cout << "Check code!" << endl;
        exit(-1);
        break;
      case 'h':
        print_usage();
        exit(-1);
        break;
      case 'v':
        verbosity=2;
        break;
      case 'b':
        verbosity=0;
        break;
      case 'f':
        fc=strtod(optarg,&endp);
        if ((optarg==endp)||(*endp!='\0')) {
          cerr << "Error: could not parse frequency" << endl;
          exit(-1);
        }
        break;
      case 'p':
        ppm=strtod(optarg,&endp);
        if ((optarg==endp)||(*endp!='\0')) {
          cerr << "Error: could not parse ppm value" << endl;
          exit(-1);
        }
        break;
      case 'c':
        correction=strtod(optarg,&endp);
        if ((optarg==endp)||(*endp!='\0')) {
          cerr << "Error: could not parse correction factor" << endl;
          exit(-1);
        }
        break;
      case 'i':
        device_index=strtol(optarg,&endp,10);
        if ((optarg==endp)||(*endp!='\0')) {
          cerr << "Error: could not parse device index" << endl;
          exit(-1);
        }
        if (device_index<0) {
          cerr << "Error: device index cannot be negative" << endl;
          exit(-1);
        }
        break;
      case '?':
        /* getopt_long already printed an error message. */
        exit(-1);
      default:
        exit(-1);
    }
  }

  // Error if extra arguments are found on the command line
  if (optind < argc) {
    cerr << "Error: unknown/extra arguments specified on command line" << endl;
    exit(-1);
  }

  // Second order command line checking. Ensure that command line options
  // are consistent.
  if (fc==-1) {
    cerr << "Error: must specify a frequency. (Try --help)" << endl;
    exit(-1);
  }
  // Start and end frequencies should be on a 100kHz raster.
  if (fc<1e6) {
    cerr << "Error: frequency must be greater than 1MHz" << endl;
    exit(-1);
  }
  if (fc/100e3!=itpp::round(fc/100e3)) {
    fc=itpp::round(fc/100e3)*100e3;
    cout << "Warning: frequency has been rounded to the nearest multiple of 100kHz" << endl;
  }
  // PPM values should be positive an most likely less than 200 ppm.
  if (ppm<0) {
    cerr << "Error: ppm value must be positive" << endl;
    exit(-1);
  }
  if (ppm>200) {
    cout << "Warning: ppm value appears to be set unreasonably high" << endl;
  }
  // Warn if correction factor is greater than 1000 ppm.
  if (abs(correction-1)>1000e-6) {
    cout << "Warning: crystal correction factor appears to be unreasonable" << endl;
  }

  if (verbosity>=1) {
    cout << "LTE CellSearch v" << MAJOR_VERSION << "." << MINOR_VERSION << "." << PATCH_LEVEL << " (" << BUILD_TYPE << ") beginning" << endl;
    cout << "  Search frequency: " << fc/1e6 << " MHz" << endl;
    cout << "  PPM: " << ppm << endl;
    stringstream temp;
    temp << setprecision(20) << correction;
    cout << "  correction: " << temp.str() << endl;
  }
}

// In high SNR environments, a cell may be detected on different carrier
// frequencies and with different frequency offsets. Keep only the cell
// with the highest received power.
void dedup(
  const list <Cell> & detected_cells,
  list <Cell> & cells_final
) {
  cells_final.clear();
  list <Cell>::const_iterator it_n=detected_cells.begin();
  while (it_n!=detected_cells.end()) {
    bool match=false;
    list <Cell>::iterator it_f=cells_final.begin();
    while (it_f!=cells_final.end()) {
      // Do these detected cells match?
      if ((*it_n).n_id_cell()==(*it_f).n_id_cell()) {
        match=true;
        // Keep either the new cell or the old cell, but not both.
        if ((*it_n).pss_pow>(*it_f).pss_pow) {
          (*it_f)=(*it_n);
        }
        break;
      }
      ++it_f;
    }
    if (!match) {
      // This cell does not match any previous cells. Add this to the
      // final list of cells.
      cells_final.push_back((*it_n));
    }
    ++it_n;
  }
}

// Helper function to assist in formatting frequency offsets.
string freq_formatter(
  const double & freq
) {
  stringstream temp;
  if (abs(freq)<998.0) {
    temp << setw(5) << setprecision(3) << freq << "h";
  } else if (abs(freq)<998000.0) {
    temp << setw(5) << setprecision(3) << freq/1e3 << "k";
  } else if (abs(freq)<998000000.0) {
    temp << setw(5) << setprecision(3) << freq/1e6 << "m";
  } else if (abs(freq)<998000000000.0) {
    temp << setw(5) << setprecision(3) << freq/1e9 << "g";
  } else if (abs(freq)<998000000000000.0) {
    temp << setw(5) << setprecision(3) << freq/1e12 << "t";
  } else {
    temp << freq;
  }
  return temp.str();
}

// Open the USB device
void config_usb(
  const double & correction,
  const int32 & device_index_cmdline,
  const double & fc
) {
  int32 device_index=device_index_cmdline;

  int8 n_rtlsdr=rtlsdr_get_device_count();
  if (n_rtlsdr==0) {
    cerr << "Error: no RTL-SDR USB devices found..." << endl;
    exit(-1);
  }

  // Choose which device to use
  if ((n_rtlsdr==1)&&(device_index==-1)) {
    device_index=0;
  }
  if ((device_index<0)||(device_index>=n_rtlsdr)) {
    cerr << "Error: must specify which USB device to use with --device-index" << endl;
    cerr << "Found the following USB devices:" << endl;
    char vendor[256],product[256],serial[256];
    for (uint8 t=0;t<n_rtlsdr;t++) {
      rtlsdr_get_device_usb_strings(t,vendor,product,serial);
      cerr << "Device index " << t << ": [Vendor: " << vendor << "] [Product: " << product << "] [Serial#: " << serial << "]" << endl;
    }
    exit(-1);
  }

  // Open device
  if (rtlsdr_open(&dev,device_index)<0) {
    cerr << "Error: unable to open RTLSDR device" << endl;
    exit(-1);
  }

  // Sampling frequency
  if (rtlsdr_set_sample_rate(dev,itpp::round(1920000*correction))<0) {
    cerr << "Error: unable to set sampling rate" << endl;
    exit(-1);
  }

  // Center frequency
  if (rtlsdr_set_center_freq(dev,itpp::round(fc*correction))<0) {
    cerr << "Error: unable to set center frequency" << endl;
    exit(-1);
  }

  // Turn on AGC
  if (rtlsdr_set_tuner_gain_mode(dev,0)<0) {
    cerr << "Error: unable to enter AGC mode" << endl;
    exit(-1);
  }

  // Reset the buffer
  if (rtlsdr_reset_buffer(dev)<0) {
    cerr << "Error: unable to reset RTLSDR buffer" << endl;
    exit(-1);
  }

  // Discard about 1.5s worth of data to give the AGC time to converge
  if (verbosity>=2) {
    cout << "Waiting for AGC to converge..." << endl;
  }
  uint32 n_read=0;
  int n_read_current;
#define BLOCK_SIZE 16*16384
  uint8 * buffer=(uint8 *)malloc(BLOCK_SIZE*sizeof(uint8));
  while (true) {
    if (rtlsdr_read_sync(dev,buffer,BLOCK_SIZE,&n_read_current)<0) {
      cerr << "Error: synchronous read failed" << endl;
      break;
    }
    if (n_read_current<BLOCK_SIZE) {
      cerr << "Error: short read; samples lost" << endl;
      break;
    }
    n_read+=n_read_current;
    if (n_read>2880000)
      break;
  }
  free(buffer);
}

// Structure to describe a cell which is currently being tracked.
typedef struct {
  boost::mutex mutex;
  boost::condition condition;
  uint16 n_id_cell;
  int8 n_ports;
  double frame_timing;
  queue <cvec> fifo;
  bool kill_me;
} tracked_cell_t;
// Structure that is used to record all the tracked cells.
typedef struct {
  // List of cells which are currently being tracked.
  // Only the searcher can add elements to this list.
  // Only the main thread can remove elements from this list.
  boost::mutex mutex;
  list <tracked_cell_t *> tracked_cells;
} tracked_cell_list_t;
// Global data shared by all threads
typedef struct {
  // The frequency offset of the dongle. This value will be updated
  // continuously.
  boost::mutex frequency_offset_mutex;
  double frequency_offset;
  // This value will never change.
  double fc;
} global_thread_data_t;
// IPC between main thread and searcher thread covering data capture issues.
typedef struct {
  boost::mutex mutex;
  boost::condition condition;
  bool request;
  cvec capbuf;
  double late;
} capbuf_sync_t;

// Wrapper around the USB device that returns samples one at a time.
class rtl_wrap {
  public:
    // Constructor
    rtl_wrap(rtlsdr_dev_t * dev);
    // Destructor
    ~rtl_wrap();
    complex <double> get_samp();
  private:
    uint8 * buffer;
    uint32 offset;
};
rtl_wrap::rtl_wrap(rtlsdr_dev_t * dev) {
  if (rtlsdr_reset_buffer(dev)<0) {
    cerr << "Error: unable to reset RTLSDR buffer" << endl;
    exit(-1);
  }
  buffer=(uint8 *)malloc(BLOCK_SIZE*sizeof(uint8));
  offset=BLOCK_SIZE;
}
rtl_wrap::~rtl_wrap() {
  free(buffer);
}
complex <double> rtl_wrap::get_samp() {
  if (offset==BLOCK_SIZE) {
    offset=0;
    int n_read;
    if (rtlsdr_read_sync(dev,buffer,BLOCK_SIZE,&n_read)<0) {
      cerr << "Error: synchronous read failed" << endl;
      exit(-1);
    }
    if (n_read<BLOCK_SIZE) {
      cerr << "Error: short read; samples lost" << endl;
      exit(-1);
    }
  }
  complex <double> samp=complex<double>((buffer[offset]-127.0)/128.0,(buffer[offset+1]-127.0)/128.0);
  offset+=2;
  return samp;
}

// Perform an initial cell search solely for the purpose of calibrating
// the oscillator.
// This code performs a full cell search on one particular center frequency.
// Future cell searches, after calibration, perform someone limited cell
// searches and this it is possible for this search to find cells that
// subsequent searches cannot.
double kalibrate(
  const double & fc,
  const double & ppm,
  const double & correction
) {
  if (verbosity>=1) {
    cout << "Calibrating local oscillator." << endl;
  }

  // Generate a list of frequency offsets that should be searched for each
  // center frequency.
  const uint16 n_extra=floor_i((fc*ppm/1e6+2.5e3)/5e3);
  const vec f_search_set=to_vec(itpp_ext::matlab_range(-n_extra*5000,5000,n_extra*5000));
  // Results are stored in this vector.
  list <Cell> detected_cells;
  // Loop until a cell is found
  while (detected_cells.size()<1) {
    // Fill capture buffer
    cvec capbuf;
    capture_data(fc,correction,false,false,".",capbuf);

    // Correlate
#define DS_COMB_ARM 2
    mat xc_incoherent_collapsed_pow;
    imat xc_incoherent_collapsed_frq;
    vf3d xc_incoherent_single;
    vf3d xc_incoherent;
    vec sp_incoherent;
    vcf3d xc;
    vec sp;
    uint16 n_comb_xc;
    uint16 n_comb_sp;
    if (verbosity>=2) {
      cout << "  Calculating PSS correlations" << endl;
    }
    xcorr_pss(capbuf,f_search_set,DS_COMB_ARM,fc,xc_incoherent_collapsed_pow,xc_incoherent_collapsed_frq,xc_incoherent_single,xc_incoherent,sp_incoherent,xc,sp,n_comb_xc,n_comb_sp);

    // Calculate the threshold vector
    const uint8 thresh1_n_nines=12;
    double R_th1=chi2cdf_inv(1-pow(10.0,-thresh1_n_nines),2*n_comb_xc*(2*DS_COMB_ARM+1));
    double rx_cutoff=(6*12*15e3/2+4*15e3)/(FS_LTE/16/2);
    vec Z_th1=R_th1*sp_incoherent/rx_cutoff/137/2/n_comb_xc/(2*DS_COMB_ARM+1);

    // Search for the peaks
    if (verbosity>=2) {
      cout << "  Searching for and examining correlation peaks..." << endl;
    }
    peak_search(xc_incoherent_collapsed_pow,xc_incoherent_collapsed_frq,Z_th1,f_search_set,fc,detected_cells);

    // Loop and check each peak
    list<Cell>::iterator iterator=detected_cells.begin();
    while (iterator!=detected_cells.end()) {
      // Detect SSS if possible
      vec sss_h1_np_est_meas;
      vec sss_h2_np_est_meas;
      cvec sss_h1_nrm_est_meas;
      cvec sss_h2_nrm_est_meas;
      cvec sss_h1_ext_est_meas;
      cvec sss_h2_ext_est_meas;
      mat log_lik_nrm;
      mat log_lik_ext;
#define THRESH2_N_SIGMA 3
      (*iterator)=sss_detect((*iterator),capbuf,THRESH2_N_SIGMA,fc,sss_h1_np_est_meas,sss_h2_np_est_meas,sss_h1_nrm_est_meas,sss_h2_nrm_est_meas,sss_h1_ext_est_meas,sss_h2_ext_est_meas,log_lik_nrm,log_lik_ext);
      if ((*iterator).n_id_1==-1) {
        // No SSS detected.
        iterator=detected_cells.erase(iterator);
        continue;
      }

      // Fine FOE
      (*iterator)=pss_sss_foe((*iterator),capbuf,fc);

      // Extract time and frequency grid
      cmat tfg;
      vec tfg_timestamp;
      extract_tfg((*iterator),capbuf,fc,tfg,tfg_timestamp);

      // Create object containing all RS
      RS_DL rs_dl((*iterator).n_id_cell(),6,(*iterator).cp_type);

      // Compensate for time and frequency offsets
      cmat tfg_comp;
      vec tfg_comp_timestamp;
      (*iterator)=tfoec((*iterator),tfg,tfg_timestamp,fc,rs_dl,tfg_comp,tfg_comp_timestamp);

      // Finally, attempt to decode the MIB
      (*iterator)=decode_mib((*iterator),tfg_comp,rs_dl);
      if ((*iterator).n_rb_dl==-1) {
        // No MIB could be successfully decoded.
        iterator=detected_cells.erase(iterator);
        continue;
      }

      if (verbosity>=2) {
        cout << "  Detected a cell!" << endl;
        cout << "    cell ID: " << (*iterator).n_id_cell() << endl;
        cout << "    RX power level: " << db10((*iterator).pss_pow) << " dB" << endl;
        cout << "    residual frequency offset: " << (*iterator).freq_superfine << " Hz" << endl;
      }

      ++iterator;
    }

    if ((detected_cells.size()==0)&&(verbosity>=1)) {
      cout << "Calibration failed (no cells detected). Trying again..." << endl;
    }
  }

  // Generate final list of detected cells. (Remove duplicates).
  list <Cell> cells_final;
  dedup(detected_cells,cells_final);

  // Search for the strongest cell. It should always be cells_final[0], but
  // perform a search just in case...
  Cell best;
  best.pss_pow=-INFINITY;
  list <Cell>::iterator it=cells_final.begin();
  while (it!=cells_final.end()) {
    if ((*it).pss_pow>best.pss_pow) {
      best=(*it);
    }
    ++it;
  }

  // Calculate the correction factor.
  // This is where we know the carrier is located
  const double true_location=best.fc;
  // We can calculate the RTLSDR's actual frequency
  const double crystal_freq_actual=best.fc-best.freq_superfine;
  // Calculate correction factors
  const double correction_residual=true_location/crystal_freq_actual;
  const double correction_new=correction*correction_residual;

  if (verbosity>=1) {
    cout << "Calibration succeeded!" << endl;
    cout << "   Residual frequency offset: " << best.freq_superfine << endl;
    cout << "   New correction factor: ";
    stringstream ss;
    ss << setprecision(20) << correction_new;
    cout << ss.str() << endl;
  }

  return best.freq_superfine;
}

// This is the searcher process. It requests captured data from the main
// thread and launches a new thread for every cell it finds. Each new
// cell thread then requests sample data from the main thread.
void searcher_proc(
  capbuf_sync_t & capbuf_sync,
  global_thread_data_t & global_thread_data,
  tracked_cell_list_t & tracked_cell_list
) {
  if (verbosity>=1) {
    cout << "Searcher process has been launched." << endl;
  }

  // Make a reference out of the pointer to clean up the code.
  double & fc=global_thread_data.fc;

  // Loop forever.
  while (true) {
    // Request data.
    {
      boost::mutex::scoped_lock lock(capbuf_sync.mutex);
      capbuf_sync.request=true;

      // Wait for data to become ready.
      capbuf_sync.condition.wait(lock);
    }

    // Get the current frequency offset
    double k_factor;
    vec f_search_set(1);
    {
      boost::mutex::scoped_lock lock(global_thread_data.frequency_offset_mutex);
      k_factor=(fc-global_thread_data.frequency_offset)/fc;
      f_search_set(0)=global_thread_data.frequency_offset;
    }

    // Results are stored in this vector.
    list<Cell> detected_cells;

    // Local reference to the capture buffer.
    cvec &capbuf=capbuf_sync.capbuf;

    // Correlate
#define DS_COMB_ARM 2
    mat xc_incoherent_collapsed_pow;
    imat xc_incoherent_collapsed_frq;
    vf3d xc_incoherent_single;
    vf3d xc_incoherent;
    vec sp_incoherent;
    vcf3d xc;
    vec sp;
    uint16 n_comb_xc;
    uint16 n_comb_sp;
    if (verbosity>=2) {
      cout << "  Calculating PSS correlations" << endl;
    }
    xcorr_pss(capbuf,f_search_set,DS_COMB_ARM,fc,xc_incoherent_collapsed_pow,xc_incoherent_collapsed_frq,xc_incoherent_single,xc_incoherent,sp_incoherent,xc,sp,n_comb_xc,n_comb_sp);

    // Calculate the threshold vector
    const uint8 thresh1_n_nines=12;
    double R_th1=chi2cdf_inv(1-pow(10.0,-thresh1_n_nines),2*n_comb_xc*(2*DS_COMB_ARM+1));
    double rx_cutoff=(6*12*15e3/2+4*15e3)/(FS_LTE/16/2);
    vec Z_th1=R_th1*sp_incoherent/rx_cutoff/137/2/n_comb_xc/(2*DS_COMB_ARM+1);

    // Search for the peaks
    if (verbosity>=2) {
      cout << "  Searching for and examining correlation peaks..." << endl;
    }
    peak_search(xc_incoherent_collapsed_pow,xc_incoherent_collapsed_frq,Z_th1,f_search_set,fc,detected_cells);

    // Loop and check each peak
    list<Cell>::iterator iterator=detected_cells.begin();
    while (iterator!=detected_cells.end()) {
      // Detect SSS if possible
      vec sss_h1_np_est_meas;
      vec sss_h2_np_est_meas;
      cvec sss_h1_nrm_est_meas;
      cvec sss_h2_nrm_est_meas;
      cvec sss_h1_ext_est_meas;
      cvec sss_h2_ext_est_meas;
      mat log_lik_nrm;
      mat log_lik_ext;
#define THRESH2_N_SIGMA 3
      (*iterator)=sss_detect((*iterator),capbuf,THRESH2_N_SIGMA,fc,sss_h1_np_est_meas,sss_h2_np_est_meas,sss_h1_nrm_est_meas,sss_h2_nrm_est_meas,sss_h1_ext_est_meas,sss_h2_ext_est_meas,log_lik_nrm,log_lik_ext);
      if ((*iterator).n_id_1==-1) {
        // No SSS detected.
        iterator=detected_cells.erase(iterator);
        continue;
      }
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
      (*iterator)=pss_sss_foe((*iterator),capbuf,fc);

      // Extract time and frequency grid
      cmat tfg;
      vec tfg_timestamp;
      extract_tfg((*iterator),capbuf,fc,tfg,tfg_timestamp);

      // Create object containing all RS
      RS_DL rs_dl((*iterator).n_id_cell(),6,(*iterator).cp_type);

      // Compensate for time and frequency offsets
      cmat tfg_comp;
      vec tfg_comp_timestamp;
      (*iterator)=tfoec((*iterator),tfg,tfg_timestamp,fc,rs_dl,tfg_comp,tfg_comp_timestamp);

      // Finally, attempt to decode the MIB
      (*iterator)=decode_mib((*iterator),tfg_comp,rs_dl);
      if ((*iterator).n_rb_dl==-1) {
        // No MIB could be successfully decoded.
        iterator=detected_cells.erase(iterator);
        continue;
      }

      if (verbosity>=1) {
        cout << "Detected a new cell!" << endl;
        cout << "  cell ID: " << (*iterator).n_id_cell() << endl;
        cout << "  RX power level: " << db10((*iterator).pss_pow) << " dB" << endl;
        cout << "  residual frequency offset: " << (*iterator).freq_superfine << " Hz" << endl;
      }

      // Launch a cell tracker process!
      tracked_cell_t * new_cell = new(tracked_cell_t);
      (*new_cell).n_id_cell=(*iterator).n_id_cell();
      (*new_cell).frame_timing=(*iterator).frame_start*k_factor+capbuf_sync.late;
      (*new_cell).n_ports=(*iterator).n_ports;
      (*new_cell).kill_me=false;
      /*if (pthread_create(&new_cell.thread,NULL,tracker_proc,(void *)(&tracked_cells))) {
        cerr << "Error: unable to launch cell tracker thread" << endl;
        exit(-1);

      }
      */
      {
        boost::mutex::scoped_lock lock(tracked_cell_list.mutex);
        tracked_cell_list.tracked_cells.push_back(new_cell);
      }

      ++iterator;
    }
  }

  // Will never reach here...
  //pthread_exit(NULL);
}

// A packet of information that is sent from the main thread to each
// tracker thread.
typedef struct {
  cvec data;
  uint8 slot_num;
  uint8 sym_num;
  double late;
} fifo_pdu;

//typedef struct {
//} tracker_shared_data_t;

// Process that tracks a cell that has been found by the searcher.
/*
void * tracker_proc(
  void * shared_temp
) {
  ivec cn=concat(itpp_ext::matlab_range(-36,-1),itpp_ext::matlab_range(1,36));

  uint8 slot_num=0;
  uint8 sym_num=0;
  // Each iteration of this loop processes one OFDM symbol.
  while (true) {
    pthread_mutex_lock(&fifo_mutex);
    if (fifo.size()==0) {
      pthread_cond_wait(&fifo_cond,&fifo_mutex);
    }
    if ((fifo.front().slot_num!=slot_num)||(fifo.front().sym_num!=sym_num)) {
      // We should never get here...
      cerr << "Error: cell tracker synchronization error! Check code!" << endl;
      exit(-1);
    }

    // Convert to frequency domain and extract 6 center RB's.
    cvec dft_out=dft(fifo.front().data);
    cvec fd=concat(dft_out.right(36),dft_out.mid(1,36));
    // Compensate for the fact that the DFT was located improperly.
    elem_mult(fd,exp((-J*2*pi*fifo.front()late/128)*cn));
    // POP the fifo
    fifo.pop();
    pthread_mutex_unlock(&fifo_mutex);

    // Increase the local counter
    sym_num=mod(sym_num+1,n_rb_dl);
    if (sym_num==0) {
      slot_num=mod(slot_num+1,20);
    }
  }
}
*/

// Main cell search routine.
int main(
  const int argc,
  char * const argv[]
) {
  // This is so that CTRL-C properly closes the rtl-sdr device before exiting
  // the program.
  struct sigaction sigact;
  sigact.sa_handler=sighandler;
  sigemptyset(&sigact.sa_mask);
  sigact.sa_flags=0;
  sigaction(SIGINT,&sigact,NULL);
  sigaction(SIGTERM,&sigact,NULL);
  sigaction(SIGQUIT,&sigact,NULL);
  sigaction(SIGPIPE,&sigact,NULL);

  // Command line parameters are stored here.
  double fc;
  double ppm;
  double correction;
  int32 device_index;

  // Get search parameters from user
  parse_commandline(argc,argv,fc,ppm,correction,device_index);

  // Open the USB device.
  config_usb(correction,device_index,fc);

  // Data shared between threads
  tracked_cell_list_t tracked_cells;
  capbuf_sync_t capbuf_sync;
  global_thread_data_t global_thread_data;

  // Calibrate the dongle's oscillator. This is similar to running the
  // program CellSearch with only one center frequency. All information
  // is discarded except for the frequency offset.
  global_thread_data.frequency_offset=kalibrate(fc,ppm,correction);

  // Start the cell searcher thread.
  // Now that the oscillator has been calibrated, we can perform
  // a 'real' search.
  capbuf_sync.request=false;
  capbuf_sync.capbuf.set_size(19200*8);
  global_thread_data.fc=fc;
  boost::thread searcher_thread(searcher_proc,boost::ref(capbuf_sync),boost::ref(global_thread_data),boost::ref(tracked_cells));
  //if (pthread_create(&searcher_thread,NULL,searcher_proc,(void *)(&tracked_cells))) {
  //  cerr << "Error: could not create searcher thread" << endl;
  //  exit(-1);
  //}

  // Wrap the USB device to simplify obtaining complex samples one by one.
  double sample_time=0;
  rtl_wrap sample_source(dev);

  // Main loop which distributes data to the appropriate subthread.
  bool searcher_capbuf_filling=false;
  uint32 searcher_capbuf_idx;
  while (true) {
    double k_factor=(fc-global_thread_data.frequency_offset)/fc;
    complex <double> sample=sample_source.get_samp();
    sample_time+=k_factor;
    sample_time=WRAP(sample_time,0.0,19200.0);

    // Should we begin filling the capture buffer?
    {
      boost::mutex::scoped_lock lock(capbuf_sync.mutex);
      if ((capbuf_sync.request)&&(!searcher_capbuf_filling)&&(abs(WRAP(sample_time-0,-19200.0/2,19200.0/2))<0.5)) {
        capbuf_sync.request=false;
        searcher_capbuf_filling=true;
        searcher_capbuf_idx=0;
        capbuf_sync.late=WRAP(sample_time-0,-19200.0/2,19200.0/2);
      }
    }

    // Populate the capture buffer
    if (searcher_capbuf_filling) {
      capbuf_sync.capbuf(searcher_capbuf_idx++)=sample;
      //cout << searcher_capbuf_idx << " " << capbuf_sync.capbuf.size() << endl;
      if (searcher_capbuf_idx==(unsigned)capbuf_sync.capbuf.size()) {
        // Buffer is full. Signal the searcher thread.
        searcher_capbuf_filling=false;
        boost::mutex::scoped_lock lock(capbuf_sync.mutex);
        capbuf_sync.condition.notify_one();
      }
    }
  }

  // Successful exit.
  exit (0);
}

