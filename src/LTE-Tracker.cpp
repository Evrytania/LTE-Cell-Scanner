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
#include <itpp/signal/transforms.h>
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
  // Hidden options. Only useful for debugging.
  //cout << "  Capture buffer options:" << endl;
  //cout << "    -l --load filename" << endl;
  //cout << "      read data from file repeatedly instead of using live data" << endl;
  //cout << "    -n --noise-power" << endl;
  //cout << "      add AWGN noise at the specified power level (in dB)
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
  int & device_index,
  bool & use_recorded_data,
  string & filename,
  double & noise_power
) {
  // Default values
  fc=-1;
  ppm=100;
  correction=1;
  device_index=-1;
  use_recorded_data=false;
  noise_power=NAN;

  while (1) {
    static struct option long_options[] = {
      {"help",         no_argument,       0, 'h'},
      {"verbose",      no_argument,       0, 'v'},
      {"brief",        no_argument,       0, 'b'},
      {"freq",         required_argument, 0, 'f'},
      {"ppm",          required_argument, 0, 'p'},
      {"correction",   required_argument, 0, 'c'},
      {"device-index", required_argument, 0, 'i'},
      {"load",         required_argument, 0, 'l'},
      {"noise-power",  required_argument, 0, 'n'},
      {0, 0, 0, 0}
    };
    /* getopt_long stores the option index here. */
    int option_index = 0;
    int c = getopt_long (argc, argv, "hvbf:p:c:i:l:n:",
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
      case 'l':
        use_recorded_data=true;
        filename=optarg;
        break;
      case 'n':
        noise_power=strtod(optarg,&endp);
        noise_power=udb10(noise_power);
        if ((optarg==endp)||(*endp!='\0')) {
          cerr << "Error: could not parse noise power" << endl;
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

//
// Data structures used to communicate between threads.
//
// A packet of information that is sent from the main thread to each
// tracker thread.
typedef struct {
  cvec data;
  uint8 slot_num;
  uint8 sym_num;
  double late;
} fifo_pdu;
// Structure to describe a cell which is currently being tracked.
class tracked_cell_t {
  public:
    // Initializer
    tracked_cell_t(
      const uint16 & n_id_cell,
      const int8 & n_ports,
      const cp_type_t::cp_type_t & cp_type,
      const double & frame_timing
    ) : n_id_cell(n_id_cell), n_ports(n_ports), cp_type(cp_type), frame_timing(frame_timing) {
      kill_me=false;
      sym_num=0;
      slot_num=0;
      target_cap_start_time=(n_symb_dl()==7)?10:32;
      filling=0;
      buffer.set_size(128);
      buffer_offset=0;
    }
    uint8 const n_symb_dl() const {
      return (cp_type==cp_type_t::NORMAL)?7:((cp_type==cp_type_t::EXTENDED)?6:-1);
    }
    boost::mutex mutex;
    boost::condition condition;
    boost::thread thread;
    // These are not allowed to change
    const uint16 n_id_cell;
    const int8 n_ports;
    const cp_type_t::cp_type_t cp_type;
    // These are constantly changing
    double frame_timing;
    queue <fifo_pdu> fifo;
    bool kill_me;
    // Calculated values returned by the cell tracker process.
    // Changing constantly.
    double power;
    // The producer process (main) will use the private members as
    // local variables...
    friend int main(const int argc,char * const argv[]);
    uint8 sym_num;
    uint8 slot_num;
    uint32 target_cap_start_time;
    double late;
    bool filling;
    cvec buffer;
    uint16 buffer_offset;
  private:
};
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
    rtl_wrap(rtlsdr_dev_t * dev,const bool & use_recorded_data,const string & filename,const double & noise_power);
    // Destructor
    ~rtl_wrap();
    complex <double> get_samp();
  private:
    bool use_recorded_data;
    cvec sig_tx;
    uint8 * buffer;
    uint32 offset;
    double noise_power_sqrt;
};
rtl_wrap::rtl_wrap(
  rtlsdr_dev_t * dev,
  const bool & urd,
  const string & filename,
  const double & noise_power
) {
  buffer=(uint8 *)malloc(BLOCK_SIZE*sizeof(uint8));
  use_recorded_data=urd;
  if (isfinite(noise_power)) {
    noise_power_sqrt=sqrt(noise_power);
  }else {
    noise_power_sqrt=NAN;
  }
  if (use_recorded_data) {
    it_ifile itf(filename);
    itf.seek("sig_tx");
    itf>>sig_tx;
    offset=0;
  } else {
    if (rtlsdr_reset_buffer(dev)<0) {
      cerr << "Error: unable to reset RTLSDR buffer" << endl;
      exit(-1);
    }
    offset=BLOCK_SIZE;
  }
}
rtl_wrap::~rtl_wrap() {
  free(buffer);
}
complex <double> rtl_wrap::get_samp() {
  complex <double> samp;
  if (use_recorded_data) {
    samp=sig_tx(offset);
    offset=mod(offset+1,length(sig_tx));
  } else {
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
    samp=complex<double>((buffer[offset]-127.0)/128.0,(buffer[offset+1]-127.0)/128.0);
    offset+=2;
  }
  if (isfinite(noise_power_sqrt)) {
    return samp+blnoise(1).get(0)*noise_power_sqrt;
  } else {
    return samp;
  }
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
  const double & correction,
  const bool & use_recorded_data,
  const string & filename,
  const double & noise_power
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
    // Fill capture buffer either from a file or from live data.
    cvec capbuf(153600);
    if (use_recorded_data) {
      cvec cbload;
      it_ifile itf(filename);
      itf.seek("sig_tx");
      itf>>cbload;
      if (length(cbload)>length(capbuf)) {
        capbuf=cbload(1,length(capbuf));
      } else {
        for (int32 t=0;t<length(capbuf);t++) {
          capbuf(t)=cbload(mod(t,length(cbload)));
        }
      }
    } else {
      capture_data(fc,correction,false,false,".",capbuf);
    }
    if (isfinite(noise_power))
      capbuf+=blnoise(length(capbuf))*sqrt(noise_power);

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

// Process that tracks a cell that has been found by the searcher.
// Data structure stored in the 'data' fifo.
typedef struct{
  uint8 slot_num;
  uint8 sym_num;
  cvec fd;
} data_fifo_pdu_t;
// Data structure used to store the 'raw' channel estimates
typedef struct {
  double shift;
  uint8 slot_num;
  uint8 sym_num;
  cvec ce;
} ce_raw_fifo_pdu_t;
// Data structure used to store the filtered channel estimates
typedef struct {
  double shift;
  uint8 slot_num;
  uint8 sym_num;
  double sp;
  double np;
  cvec ce_filt;
} ce_filt_fifo_pdu_t;
void tracker_proc(
  tracked_cell_t & tracked_cell,
  global_thread_data_t & global_thread_data
) {
  ivec cn=concat(itpp_ext::matlab_range(-36,-1),itpp_ext::matlab_range(1,36));
  RS_DL rs_dl(tracked_cell.n_id_cell,6,tracked_cell.cp_type);

  uint8 slot_num=0;
  uint8 sym_num=0;
  double bulk_phase_offset=0;
  deque <data_fifo_pdu_t> data_fifo;
  vector <deque <ce_raw_fifo_pdu_t> > ce_raw_fifo(tracked_cell.n_ports);
  vector <deque <ce_filt_fifo_pdu_t> > ce_filt_fifo(tracked_cell.n_ports);
  // Each iteration of this loop processes one OFDM symbol.
  while (true) {
    cvec fd;
    double frequency_offset;
    double k_factor;
    {
      boost::mutex::scoped_lock lock(tracked_cell.mutex);
      if (tracked_cell.fifo.empty()) {
        tracked_cell.condition.wait(lock);
      }
      ASSERT(!tracked_cell.fifo.empty());
      if ((tracked_cell.fifo.front().slot_num!=slot_num)||(tracked_cell.fifo.front().sym_num!=sym_num)) {
        // We should never get here...
        cerr << "Error: cell tracker synchronization error! Check code!" << endl;
        exit(-1);
      }

      // Convert to frequency domain and extract 6 center RB's.
      {
        boost::mutex::scoped_lock lock(global_thread_data.frequency_offset_mutex);
        frequency_offset=global_thread_data.frequency_offset;
        k_factor=(global_thread_data.fc-frequency_offset)/global_thread_data.fc;
      }
      // How many time samples have passed since the previous DFT?
      // Also perform FOC to remove ICI
      cvec dft_in=fshift(tracked_cell.fifo.front().data,-frequency_offset,FS_LTE/16);
      // Remove the 2 sample delay
      dft_in=concat(dft_in(2,-1),dft_in(0,1));
      cvec dft_out=dft(fshift(tracked_cell.fifo.front().data,-frequency_offset,FS_LTE/16));
      fd=concat(dft_out.right(36),dft_out.mid(1,36));
      // Compensate for the fact that the DFT was located improperly and also
      // for the bulk phase offset due to frequency errors.
      uint8 n_samp_elapsed;
      if (tracked_cell.cp_type==cp_type_t::EXTENDED) {
        n_samp_elapsed=128+32;
      } else {
        n_samp_elapsed=(sym_num==0)?128+10:128+9;
      }
      bulk_phase_offset=WRAP(bulk_phase_offset+2*pi*n_samp_elapsed*k_factor*(1/(FS_LTE/16))*-frequency_offset,-pi,pi);
      fd=exp(J*bulk_phase_offset)*elem_mult(fd,exp((-J*2*pi*tracked_cell.fifo.front().late/128)*cn));
      // POP the fifo
      tracked_cell.fifo.pop();
    }
    // At this point, we have the frequency domain data for this slot and
    // this symbol number. FOC and TOC has already been performed.

    // Save this information into the data fifo for further processing
    // once the channel estimates are ready for this ofdm symbol.
    data_fifo_pdu_t dfp;
    dfp.slot_num=slot_num;
    dfp.sym_num=sym_num;
    dfp.fd=fd;
    data_fifo.push_back(dfp);

    // Extract any RS that might be present and perform TOE.
    for (uint8 port_num=0;port_num<tracked_cell.n_ports;port_num++) {
      double shift=rs_dl.get_shift(slot_num,sym_num,port_num);
      if (isnan(shift))
        continue;
      cvec rs_raw=fd(itpp_ext::matlab_range(round_i(shift),6,71));
      cvec ce_raw=elem_mult(rs_raw,conj(rs_dl.get_rs(slot_num,sym_num)));
      ce_raw_fifo_pdu_t cerp;
      cerp.shift=shift;
      cerp.slot_num=slot_num;
      cerp.sym_num=sym_num;
      cerp.ce=ce_raw;
      ce_raw_fifo[port_num].push_back(cerp);
    }

    // Filter and perform FOE and TOE on the raw channel estimates
    for (uint8 port_num=0;port_num<tracked_cell.n_ports;port_num++) {
      if (ce_raw_fifo[port_num].size()!=3)
        continue;

      // Perform primitive filtering by averaging nearby samples.
      cvec ce_filt(12);
      for (uint8 t=0;t<12;t++) {
        complex <double> total=0;
        uint8 n_total=0;
        ivec ind;
        ind=itpp_ext::matlab_range(t-1,t+1);
        del_oob(ind);
        total=sum(ce_raw_fifo[port_num][1].ce.get(ind));
        n_total=length(ind);
        if (ce_raw_fifo[port_num][0].shift<ce_raw_fifo[port_num][1].shift) {
          ind=itpp_ext::matlab_range(t,t+1);
        } else {
          ind=itpp_ext::matlab_range(t-1,t);
        }
        del_oob(ind);
        total=total+sum(ce_raw_fifo[port_num][0].ce.get(ind));
        total=total+sum(ce_raw_fifo[port_num][2].ce.get(ind));
        n_total=2*length(ind);
        ce_filt(t)=total/n_total;
      }
      // Store filtered channel estimates.
      ce_filt_fifo_pdu_t pdu;
      pdu.shift=ce_raw_fifo[port_num][1].shift;
      pdu.slot_num=ce_raw_fifo[port_num][1].slot_num;
      pdu.sym_num=ce_raw_fifo[port_num][1].sym_num;
      pdu.sp=sigpower(ce_filt);
      pdu.np=sigpower(ce_raw_fifo[port_num][1].ce-ce_filt);
      pdu.ce_filt=ce_filt;
      ce_filt_fifo[port_num].push_back(pdu);

      // FOE
      cvec foe=elem_mult(conj(ce_raw_fifo[port_num][0].ce),ce_raw_fifo[port_num][2].ce);
      // Calculate the noise on each FOE estimate.
      vec foe_np=pdu.np*pdu.np+2*pdu.np*sqr(ce_filt);
      // Calculate the weight to use for each estimate
      vec weight=elem_div(sqr(ce_filt),foe_np);
      // MRC
      complex <double> foe_comb=sum(elem_mult(foe,to_cvec(weight)));
      double foe_comb_np=sum(elem_mult(foe_np,weight,weight));
      // Scale. Only necessary for NP.
      double scale=1/sum(elem_mult(sqr(ce_filt),weight));
      foe_comb=foe_comb*scale;
      foe_comb_np=foe_comb_np*scale*scale;

      // Update system frequency offset.
      double residual_f=arg(foe_comb)/(2*pi)/(k_factor*.0005);
      double residual_f_np=foe_comb_np/2;
      {
        boost::mutex::scoped_lock lock(global_thread_data.frequency_offset_mutex);
        global_thread_data.frequency_offset=(
          global_thread_data.frequency_offset*(1/.0001)+
          (global_thread_data.frequency_offset+residual_f)*(1/residual_f_np)
        )/( 1/.0001+1/residual_f_np);
        frequency_offset=global_thread_data.frequency_offset;
        k_factor=(global_thread_data.fc-frequency_offset)/global_thread_data.fc;
        cout << "FO: " << frequency_offset << endl;
      }

      // TOE
      cvec toe=concat(
        elem_mult(conj(ce_raw_fifo[port_num][1].ce(0,4)),ce_raw_fifo[port_num][1].ce(1,5)),
        elem_mult(conj(ce_raw_fifo[port_num][1].ce(6,10)),ce_raw_fifo[port_num][1].ce(7,11))
      );
      // Calculate the noise on each TOE estimate.
      vec toe_np=pdu.np*pdu.np+2*pdu.np*sqr(concat(ce_filt(0,4),ce_filt(6,10)));
      // Calculate the weight to use for each estimate
      weight=elem_div(sqr(concat(ce_filt(0,4),ce_filt(6,10))),toe_np);
      // MRC
      complex <double> toe_comb=sum(elem_mult(toe,to_cvec(weight)));
      double toe_comb_np=sum(elem_mult(toe_np,weight,weight));
      // Scale. Only necessary for NP.
      scale=1/sum(elem_mult(sqr(concat(ce_filt(0,4),ce_filt(6,10))),weight));
      toe_comb=toe_comb*scale;
      toe_comb_np=toe_comb_np*scale*scale;
      double delay=-arg(toe_comb)/6/(2*pi/128);
      double delay_np=toe_comb_np/2;

      // Update frame timing based on TOE
      {
        boost::mutex::scoped_lock lock(tracked_cell.mutex);
        tracked_cell.frame_timing=(
          tracked_cell.frame_timing*(1/.0001)+
          (tracked_cell.frame_timing+delay)*(1/delay_np)
        )/( 1/.0001+1/delay_np);
        tracked_cell.frame_timing=itpp_ext::matlab_mod(tracked_cell.frame_timing,19200.0);
        cout << "TO: " << tracked_cell.frame_timing << endl;
      }

      // Finished with the raw channel estimates.
      ce_raw_fifo[port_num].pop_front();

      // Interpolate filtered channel estimates
      while (!ce_filt_fifo[port_num].empty())
        ce_filt_fifo[port_num].pop_front();

    }

    // Process data now that channel estimates are available for every
    // data sample.
    while (!data_fifo.empty())
      data_fifo.pop_front();

    // Increase the local counter
    sym_num=mod(sym_num+1,tracked_cell.n_symb_dl());
    if (sym_num==0) {
      slot_num=mod(slot_num+1,20);
    }
  }
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
      tracked_cell_t * new_cell = new tracked_cell_t((*iterator).n_id_cell(),(*iterator).n_ports,(*iterator).cp_type,(*iterator).frame_start*k_factor+capbuf_sync.late);
      (*new_cell).thread=boost::thread(tracker_proc,boost::ref(*new_cell),boost::ref(global_thread_data));
      {
        boost::mutex::scoped_lock lock(tracked_cell_list.mutex);
        tracked_cell_list.tracked_cells.push_back(new_cell);
      }
      cout << "Only one cell is allowed to be detected!!!" << endl;
      sleep(1000000);

      ++iterator;
    }
  }
  // Will never reach here...
}

// Main routine.
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
  bool use_recorded_data;
  string filename;
  double noise_power;

  // Get search parameters from the user
  parse_commandline(argc,argv,fc,ppm,correction,device_index,use_recorded_data,filename,noise_power);

  // Open the USB device.
  if (!use_recorded_data)
    config_usb(correction,device_index,fc);

  // Data shared between threads
  tracked_cell_list_t tracked_cell_list;
  capbuf_sync_t capbuf_sync;
  global_thread_data_t global_thread_data;

  // Calibrate the dongle's oscillator. This is similar to running the
  // program CellSearch with only one center frequency. All information
  // is discarded except for the frequency offset.
  global_thread_data.frequency_offset=kalibrate(fc,ppm,correction,use_recorded_data,filename,noise_power);

  // Start the cell searcher thread.
  // Now that the oscillator has been calibrated, we can perform
  // a 'real' search.
  capbuf_sync.request=false;
  capbuf_sync.capbuf.set_size(19200*8);
  global_thread_data.fc=fc;
  boost::thread searcher_thread(searcher_proc,boost::ref(capbuf_sync),boost::ref(global_thread_data),boost::ref(tracked_cell_list));

  // Wrap the USB device to simplify obtaining complex samples one by one.
  double sample_time=0;
  rtl_wrap sample_source(dev,use_recorded_data,filename,noise_power);

  // Main loop which distributes data to the appropriate subthread.
  bool searcher_capbuf_filling=false;
  uint32 searcher_capbuf_idx;
  while (true) {
    // Each iteration of this loop processes one sample.
    double k_factor;
    {
      boost::mutex::scoped_lock lock(global_thread_data.frequency_offset_mutex);
      k_factor=(fc-global_thread_data.frequency_offset)/fc;
    }
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
      if (searcher_capbuf_idx==(unsigned)capbuf_sync.capbuf.size()) {
        // Buffer is full. Signal the searcher thread.
        searcher_capbuf_filling=false;
        boost::mutex::scoped_lock lock(capbuf_sync.mutex);
        capbuf_sync.condition.notify_one();
      }
    }

    // Loop for each tracked cell and save data, if necessary.
    {
      boost::mutex::scoped_lock lock(tracked_cell_list.mutex);
      list <tracked_cell_t *>::iterator it=tracked_cell_list.tracked_cells.begin();
      while (it!=tracked_cell_list.tracked_cells.end()) {
        tracked_cell_t & tracked_cell=(*(*it));
        boost::mutex::scoped_lock lock2(tracked_cell.mutex);

        // See if we should start filling the buffer.
        if (!tracked_cell.filling) {
          double tdiff=WRAP(sample_time-(tracked_cell.frame_timing+tracked_cell.target_cap_start_time),-19200.0/2,19200.0/2);
          if (
            // Ideal start time is 0.5 samples away from current time
            (abs(tdiff)<0.5) ||
            // It's possible for the frame timing to change between iterations
            // of the outer sample loop and because of this, it's possible that
            // we missed the best start. Start capturing anyways.
            ((tdiff>0)&&(tdiff<2))
          ) {
            // Configure parameters for this capture
            tracked_cell.filling=true;
            tracked_cell.late=tdiff;
            tracked_cell.buffer_offset=0;
          }
        }

        // Save this sample if our state indicates we are filling the
        // buffer.
        if (tracked_cell.filling) {
          tracked_cell.buffer(tracked_cell.buffer_offset++)=sample;
          if (tracked_cell.buffer_offset==128) {
            // Buffer is full!
            // Send PDU
            fifo_pdu p;
            p.data=tracked_cell.buffer;
            p.slot_num=tracked_cell.slot_num;
            p.sym_num=tracked_cell.sym_num;
            p.late=tracked_cell.late;
            tracked_cell.fifo.push(p);
            tracked_cell.condition.notify_one();
            // Calculate trigger parameters of next capture
            tracked_cell.filling=false;
            if (tracked_cell.n_symb_dl()==6) {
              tracked_cell.target_cap_start_time+=32+128;
            } else {
              tracked_cell.target_cap_start_time+=(tracked_cell.sym_num==6)?128+10:128+9;
            }
            tracked_cell.target_cap_start_time=mod(tracked_cell.target_cap_start_time,19200);
            tracked_cell.sym_num=mod(tracked_cell.sym_num+1,tracked_cell.n_symb_dl());
            if (tracked_cell.sym_num==0) {
              tracked_cell.slot_num=mod(tracked_cell.slot_num+1,20);
            }
          }
        }
        ++it;
      }
    }
  }

  // Successful exit.
  exit (0);
}

