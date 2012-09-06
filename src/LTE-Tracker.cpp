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

// Global variables that can be set by the command line. Used for debugging.
double global_1=0;
double global_2=0;
double global_3=0;
double global_4=0;
double global_5=0;
double global_6=0;
double global_7=0;
double global_8=0;
double global_9=0;

/*
static void sighandler(
  int signum
) {
  cerr << "Error: caught signal, exiting!" << endl;
  if (dev!=NULL) {
    rtlsdr_close(dev);
  }
  exit(-1);
}
*/

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
  //cout << "    -r --repeat" << endl;
  //cout << "      cyclically repeat the data read from the file forever" << endl;
  //cout << "    -s --rtl_sdr" << endl;
  //cout << "      data file was created by the rtl_sdr program. Defaults to an itpp file." << endl;
  //cout << "    -n --noise-power value" << endl;
  //cout << "      add AWGN noise at the specified power level (in dB)
  //cout << "  Global variables used for testing" << endl;
  //cout << "    -1 --g1 value" << endl;
  //cout << "    -2 --g2 value" << endl;
  //cout << "    -3 --g3 value" << endl;
  //cout << "    -4 --g4 value" << endl;
  //cout << "    -5 --g5 value" << endl;
  //cout << "    -6 --g6 value" << endl;
  //cout << "    -7 --g7 value" << endl;
  //cout << "    -8 --g8 value" << endl;
  //cout << "    -9 --g9 value" << endl;
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
  bool & repeat,
  bool & rtl_sdr_format,
  double & noise_power
) {
  // Default values
  fc=-1;
  ppm=100;
  correction=1;
  device_index=-1;
  use_recorded_data=false;
  repeat=false;
  rtl_sdr_format=false;
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
      {"repeat",       no_argument,       0, 'r'},
      {"rtl_sdr",      no_argument,       0, 's'},
      {"noise-power",  required_argument, 0, 'n'},
      {"g1",           required_argument, 0, '1'},
      {"g2",           required_argument, 0, '2'},
      {"g3",           required_argument, 0, '3'},
      {"g4",           required_argument, 0, '4'},
      {"g5",           required_argument, 0, '5'},
      {"g6",           required_argument, 0, '6'},
      {"g7",           required_argument, 0, '7'},
      {"g8",           required_argument, 0, '8'},
      {"g9",           required_argument, 0, '9'},
      {0, 0, 0, 0}
    };
    /* getopt_long stores the option index here. */
    int option_index = 0;
    int c = getopt_long (argc, argv, "hvbf:p:c:i:l:rsn:123456789",
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
      case 'r':
        repeat=true;
        break;
      case 's':
        rtl_sdr_format=true;
        break;
      case 'n':
        noise_power=strtod(optarg,&endp);
        noise_power=udb10(noise_power);
        if ((optarg==endp)||(*endp!='\0')) {
          cerr << "Error: could not parse noise power" << endl;
          exit(-1);
        }
        break;
      case '1':
        global_1=strtod(optarg,&endp);
        if ((optarg==endp)||(*endp!='\0')) {
          cerr << "Error: could not global variable 1" << endl;
          exit(-1);
        }
        break;
      case '2':
        global_2=strtod(optarg,&endp);
        if ((optarg==endp)||(*endp!='\0')) {
          cerr << "Error: could not global variable 2" << endl;
          exit(-1);
        }
        break;
      case '3':
        global_3=strtod(optarg,&endp);
        if ((optarg==endp)||(*endp!='\0')) {
          cerr << "Error: could not global variable 3" << endl;
          exit(-1);
        }
        break;
      case '4':
        global_4=strtod(optarg,&endp);
        if ((optarg==endp)||(*endp!='\0')) {
          cerr << "Error: could not global variable 4" << endl;
          exit(-1);
        }
        break;
      case '5':
        global_5=strtod(optarg,&endp);
        if ((optarg==endp)||(*endp!='\0')) {
          cerr << "Error: could not global variable 5" << endl;
          exit(-1);
        }
        break;
      case '6':
        global_6=strtod(optarg,&endp);
        if ((optarg==endp)||(*endp!='\0')) {
          cerr << "Error: could not global variable 6" << endl;
          exit(-1);
        }
        break;
      case '7':
        global_7=strtod(optarg,&endp);
        if ((optarg==endp)||(*endp!='\0')) {
          cerr << "Error: could not global variable 7" << endl;
          exit(-1);
        }
        break;
      case '8':
        global_8=strtod(optarg,&endp);
        if ((optarg==endp)||(*endp!='\0')) {
          cerr << "Error: could not global variable 8" << endl;
          exit(-1);
        }
        break;
      case '9':
        global_9=strtod(optarg,&endp);
        if ((optarg==endp)||(*endp!='\0')) {
          cerr << "Error: could not global variable 9" << endl;
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
  if (optind<argc) {
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
    cout << "LTE Tracker v" << MAJOR_VERSION << "." << MINOR_VERSION << "." << PATCH_LEVEL << " (" << BUILD_TYPE << ") beginning" << endl;
    cout << "  Search frequency: " << fc/1e6 << " MHz" << endl;
    cout << "  PPM: " << ppm << endl;
    stringstream temp;
    temp << setprecision(20) << correction;
    cout << "  correction: " << temp.str() << endl;
  }
}

// Small helper function to increment the slot number and the symbol number.
void slot_sym_inc(
  const uint8 n_symb_dl,
  uint8 & slot_num,
  uint8 & sym_num
) {
  sym_num=mod(sym_num+1,n_symb_dl);
  if (sym_num==0)
    slot_num=mod(slot_num+1,20);
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

  // Discard about 3s worth of data to give the AGC time to converge
  if (verbosity>=2) {
    cout << "Waiting for AGC to converge..." << endl;
  }
  uint32 n_read=0;
  int n_read_current;
#define BLOCK_SIZE (16*16384)
  uint8 * buffer=(uint8 *)malloc(BLOCK_SIZE*sizeof(uint8));
  while (true) {
    if (rtlsdr_read_sync(dev,buffer,BLOCK_SIZE,&n_read_current)<0) {
      cerr << "Error: synchronous read failed" << endl;
      exit(-1);
    }
    if (n_read_current<BLOCK_SIZE) {
      cerr << "Error: short read; samples lost" << endl;
      exit(-1);
    }
    n_read+=n_read_current;
    if (n_read>2880000*2)
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
  double frequency_offset;
  double frame_timing;
} td_fifo_pdu_t;
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
      target_cap_start_time=(cp_type==cp_type_t::NORMAL)?10:32;
      filling=0;
      buffer.set_size(128);
      buffer_offset=0;
      ac_fd.set_size(12);
      ac_fd=complex <double> (0,0);
      bulk_phase_offset=0;
      tracker_proc_ready=false;
    }
    uint8 const n_symb_dl() const {
      return (cp_type==cp_type_t::NORMAL)?7:((cp_type==cp_type_t::EXTENDED)?6:-1);
    }
    boost::mutex mutex;
    boost::condition condition;
    boost::thread thread;
    // Indicates that the tracker process is ready to receive data.
    bool tracker_proc_ready;
    // These are not allowed to change
    const uint16 n_id_cell;
    const int8 n_ports;
    const cp_type_t::cp_type_t cp_type;
    // These are constantly changing
    double frame_timing;
    queue <td_fifo_pdu_t> fifo;
    bool kill_me;
    // Calculated values returned by the cell tracker process.
    // Changing constantly.
    cvec ac_fd;
    // Tracker process uses these as local variables.
    double bulk_phase_offset;
    // The producer process (main) will use these members as
    // local variables...
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
// IPC between main thread and producer thread.
typedef struct {
  boost::mutex mutex;
  boost::condition condition;
  deque <uint8> fifo;
} sampbuf_sync_t;

void read_datafile(
  const string & filename,
  const bool & rtl_sdr_format,
  cvec & sig_tx
) {
  if (!rtl_sdr_format) {
    cout << "Trying itpp format..." << endl;
    it_ifile itf(filename);
    itf.seek("sig_tx");
    itf>>sig_tx;
    cout << "read itpp format" << endl;
    return;
  } else {
    cvec sig_tx_pre;
    itpp_ext::rtl_sdr_to_cvec(filename,sig_tx_pre);
    sig_tx.set_size(length(sig_tx_pre)*2-10);
    for (int32 t=0;t<length(sig_tx);t+=2) {
      // FIXME: Do proper interpolation
      sig_tx(t)=sig_tx_pre(t>>1);
      sig_tx(t+1)=(sig_tx_pre(t>>1)+sig_tx_pre((t>>1)+1))/2;
    }
    // Drop several seconds while AGC converges.
    sig_tx=sig_tx(FS_LTE/16*4,-1);
  }
  if (length(sig_tx)==0) {
    cerr << "Error: no data in file!" << endl;
    exit(-1);
  }
}

// Wrapper around the USB device that returns samples one at a time.
class rtl_wrap {
  public:
    // Constructor
    rtl_wrap(rtlsdr_dev_t * dev,const bool & use_recorded_data,const string & filename,const bool & repeat,const bool & rtl_sdr_format,const double & noise_power);
    // Destructor
    ~rtl_wrap();
    complex <double> get_samp();
    //void reset();
  private:
    complex <double> get_samp_pre();
    bool use_recorded_data;
    cvec sig_tx;
    uint8 * buffer;
    uint32 offset;
    double noise_power_sqrt;
    bool repeat;
    bool rtl_sdr_format;
    bool phase_even;
    complex <double> samp_d1;
    complex <double> samp_d2;
    //FILE * file;
};
rtl_wrap::rtl_wrap(
  rtlsdr_dev_t * dev,
  const bool & urd,
  const string & filename,
  const bool & rpt,
  const bool & rsdf,
  const double & noise_power
) {
  buffer=(uint8 *)malloc(BLOCK_SIZE*sizeof(uint8));
  use_recorded_data=urd;
  repeat=rpt;
  rtl_sdr_format=rsdf;
  phase_even=true;
  samp_d1=complex <double> (0,0);
  samp_d2=complex <double> (0,0);
  //cout << "Opening log file" << endl;
  //file=fopen("tracker_log.dat","wb");
  //if (!file) {
  //  cerr << "Error: could not open ddata capture log file" << endl;
  //  exit(-1);
  //}
  //cout << "Log file opened!" << endl;
  if (isfinite(noise_power)) {
    noise_power_sqrt=sqrt(noise_power);
  } else {
    noise_power_sqrt=NAN;
  }
  if (use_recorded_data) {
    // Note that the entire file is read into memory and stored as a complex
    // double!
    read_datafile(filename,rtl_sdr_format,sig_tx);
    //it_ifile itf(filename);
    //itf.seek("sig_tx");
    //itf>>sig_tx;
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
complex <double> rtl_wrap::get_samp_pre() {
  if (offset==BLOCK_SIZE) {
    offset=0;
    int n_read;
    if (rtlsdr_read_sync(dev,buffer,BLOCK_SIZE,&n_read)<0) {
      cerr << "Error: synchronous read failed" << endl;
      exit(-1);
    }
    //if (fwrite(buffer, 1, n_read, file)!=(size_t)n_read) {
    //  cerr<<"Error: Short write, samples lost, exiting!" << endl;
    //  exit(-1);
    //}
    if (n_read<BLOCK_SIZE) {
      cerr << "Error: short read; samples lost" << endl;
      exit(-1);
    }
  }
  complex <double>samp=complex<double>((buffer[offset]-127.0)/128.0,(buffer[offset+1]-127.0)/128.0);
  offset+=2;
  return samp;
}
complex <double> rtl_wrap::get_samp() {
  static long long cnt=0;

  complex <double> samp;
  if (use_recorded_data) {
    //samp=sig_tx(offset);
    samp=1*sig_tx(offset)+(0.0*J)*exp(J*((double)(cnt++))*2*pi/10000000)*sig_tx(mod(offset-5000,length(sig_tx)));;
    offset=mod(offset+1,length(sig_tx));
    if ((offset==0)&&(!repeat)) {
      // Not that if N complex samples are read from the file and repeat is
      // false, N-1 samples will be produced by this class before program
      // execution terminates.
      cerr << "Error: no more sample data in file!" << endl;
      exit(-1);
    }
  } else {
    // This code is wrong! FIXME
    if (phase_even==true) {
      samp_d2=samp_d1;
      samp_d1=get_samp_pre();
      samp=samp_d1;
    } else {
      samp=(samp_d1+samp_d2)/2;
    }
    phase_even=!phase_even;
    return samp;
  }
  if (isfinite(noise_power_sqrt)) {
    return samp+blnoise(1).get(0)*noise_power_sqrt;
  } else {
    return samp;
  }
}
/*
complex <double> rtl_wrap::get_samp() {
}
*/
/*
void rtl_wrap::reset() {
  offset=0;
}
*/

// Perform an initial cell search solely for the purpose of calibrating
// the oscillator.
// This code can probably be rolled into the searcher so as to eliminate
// some duplicated code.
double kalibrate(
  const double & fc,
  const double & ppm,
  const double & correction,
  const bool & use_recorded_data,
  const string & filename,
  const bool & rtl_sdr_format,
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
    cvec capbuf;
    // Fill capture buffer either from a file or from live data.
    if (use_recorded_data) {
      capbuf.set_size(153600);
      cvec cbload;
      read_datafile(filename,rtl_sdr_format,cbload);
      //it_ifile itf(filename);
      //itf.seek("sig_tx");
      //itf>>cbload;
      if (length(cbload)>=length(capbuf)) {
        capbuf=cbload(1,length(capbuf));
      } else {
        for (int32 t=0;t<length(capbuf);t++) {
          capbuf(t)=cbload(mod(t,length(cbload)));
        }
      }
    } else {
      capture_data(fc,correction,false,false,".",capbuf);
    }
    //cout << "Capbuf power: " << db10(sigpower(capbuf)) << " dB" << endl;
    if (isfinite(noise_power))
      capbuf+=blnoise(length(capbuf))*sqrt(noise_power);

    // Correlate
    //cout << "DS_COMB_ARM override!!!" << endl;
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
    peak_search(xc_incoherent_collapsed_pow,xc_incoherent_collapsed_frq,Z_th1,f_search_set,fc,xc_incoherent_single,DS_COMB_ARM,detected_cells);

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

// Data structure stored in the 'data' fifo.
typedef struct{
  uint8 slot_num;
  uint8 sym_num;
  cvec syms;
} data_fifo_pdu_t;
// Data structure used to store the 'raw' channel estimates
typedef struct {
  double shift;
  uint8 slot_num;
  uint8 sym_num;
  cvec ce;
  double frequency_offset;
  double frame_timing;
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
typedef struct {
  uint8 slot_num;
  uint8 sym_num;
  double sp;
  double np;
  cvec ce_interp;
} ce_interp_fifo_pdu_t;
typedef struct {
  cvec syms;
  cmat ce;
  vec sp;
  vec np;
} mib_fifo_pdu_t;

void get_fd(
  tracked_cell_t & tracked_cell,
  const double & fc,
  const uint8 & slot_num,
  const uint8 & sym_num,
  const ivec & cn,
  double & bulk_phase_offset,
  cvec & syms,
  double & frequency_offset,
  double & frame_timing
) {
  // Lock the tracked_cell data until we obtain one element from
  // the fifo and convert to the frequency domain.
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
  frequency_offset=tracked_cell.fifo.front().frequency_offset;
  frame_timing=tracked_cell.fifo.front().frame_timing;
  double k_factor=(fc-frequency_offset)/fc;
  // Also perform FOC to remove ICI
  cvec dft_in=fshift(tracked_cell.fifo.front().data,-frequency_offset,FS_LTE/16*k_factor);
  // Remove the 2 sample delay
  dft_in=concat(dft_in(2,-1),dft_in(0,1));
  cvec dft_out=dft(dft_in);
  syms=concat(dft_out.right(36),dft_out.mid(1,36));
  // Compensate for the fact that the DFT was located improperly and also
  // for the bulk phase offset due to frequency errors.
  uint8 n_samp_elapsed;
  // How many time samples have passed since the previous DFT?
  if (tracked_cell.cp_type==cp_type_t::EXTENDED) {
    n_samp_elapsed=128+32;
  } else {
    n_samp_elapsed=(sym_num==0)?128+10:128+9;
  }
  bulk_phase_offset=WRAP(bulk_phase_offset+2*pi*n_samp_elapsed*(1/(FS_LTE/16))*-frequency_offset,-pi,pi);
  syms=exp(J*bulk_phase_offset)*elem_mult(syms,exp((-J*2*pi*tracked_cell.fifo.front().late/128)*cn));
  // POP the fifo
  tracked_cell.fifo.pop();
  // At this point, we have the frequency domain data for this slot and
  // this symbol number. FOC and TOC has already been performed.
}

cvec filter_ce(
  const ce_raw_fifo_pdu_t & rs_prev,
  const ce_raw_fifo_pdu_t & rs_curr,
  const ce_raw_fifo_pdu_t & rs_next
) {
  cvec ce_filt(12);
  for (uint8 t=0;t<12;t++) {
    complex <double> total=0;
    uint8 n_total=0;
    ivec ind;
    ind=itpp_ext::matlab_range(t-1,t+1);
    del_oob(ind);
    total=sum(rs_curr.ce.get(ind));
    n_total=length(ind);
    if (rs_prev.shift<rs_curr.shift) {
      ind=itpp_ext::matlab_range(t,t+1);
    } else {
      ind=itpp_ext::matlab_range(t-1,t);
    }
    del_oob(ind);
    total=total+sum(rs_prev.ce.get(ind));
    total=total+sum(rs_next.ce.get(ind));
    n_total+=2*length(ind);
    ce_filt(t)=total/n_total;
  }
  return ce_filt;
}

void do_foe(
  global_thread_data_t & global_thread_data,
  const ce_raw_fifo_pdu_t & rs_prev,
  const ce_raw_fifo_pdu_t & rs_next,
  const double & rs_curr_np,
  const cvec & ce_filt
) {
  cvec foe=elem_mult(conj(rs_prev.ce),rs_next.ce);
  // Calculate the noise on each FOE estimate.
  vec foe_np=rs_curr_np*rs_curr_np+2*rs_curr_np*sqr(ce_filt);
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
  double frequency_offset=rs_prev.frequency_offset;
  //double k_factor=(global_thread_data.fc-frequency_offset)/global_thread_data.fc;
  double residual_f=arg(foe_comb)/(2*pi)/0.0005;
  double residual_f_np=MAX(foe_comb_np/2,.001);
  //residual_f+=1000;
  //cout << residual_f << " f " << db10(residual_f_np) << endl;
  //residual_f=0;
  {
    boost::mutex::scoped_lock lock(global_thread_data.frequency_offset_mutex);
    global_thread_data.frequency_offset=(
      global_thread_data.frequency_offset*(1/.0001)+
      (frequency_offset+residual_f)*(1/residual_f_np)
    )/(1/.0001+1/residual_f_np);
    //cout << "FO: " << frequency_offset << endl;
  }
}

void do_toe_v2(
  tracked_cell_t & tracked_cell,
  const ce_raw_fifo_pdu_t & rs_prev,
  const ce_raw_fifo_pdu_t & rs_curr,
  const double & rs_curr_sp,
  const double & rs_curr_np
) {
  complex <double> toe1;
  complex <double> toe2;
  if (rs_prev.shift<rs_curr.shift) {
    toe1=sum(elem_mult(conj(rs_prev.ce),rs_curr.ce))/12;
    toe2=(
      sum(elem_mult(conj(rs_curr.ce(0,4)),rs_prev.ce(1,5)))+
      sum(elem_mult(conj(rs_curr.ce(6,10)),rs_prev.ce(7,11)))
    )/10;
  } else {
    toe1=sum(elem_mult(conj(rs_curr.ce),rs_prev.ce))/12;
    toe2=(
      sum(elem_mult(conj(rs_prev.ce(0,4)),rs_curr.ce(1,5)))+
      sum(elem_mult(conj(rs_prev.ce(6,10)),rs_curr.ce(7,11)))
    )/10;
  }
  toe1=toe1/sqrt(rs_curr_sp);
  toe2=toe2/sqrt(rs_curr_sp);
  double delay=-(arg(toe1)+arg(toe2))/2/3/(2*pi/128);
  double delay_np=MAX(rs_curr_np/rs_curr_sp/2/12,.001);

  // Update frame timing based on TOE
  {
    boost::mutex::scoped_lock lock(tracked_cell.mutex);
    double diff=WRAP((rs_curr.frame_timing+delay)-tracked_cell.frame_timing,-19200.0/2,19200.0/2);
    diff=(0*(1/.0001)+diff*(1/delay_np))/(1/.0001+1/delay_np);
    tracked_cell.frame_timing=itpp_ext::matlab_mod(tracked_cell.frame_timing+diff,19200.0);
    //cout << "TO: " << setprecision(15) << tracked_cell.frame_timing << endl;
  }
}

void do_toe(
  tracked_cell_t & tracked_cell,
  const ce_raw_fifo_pdu_t & rs_curr,
  const cvec & rs_curr_filt,
  const double & rs_curr_np
) {
  cvec toe=concat(
    elem_mult(conj(rs_curr.ce(0,4)),rs_curr.ce(1,5)),
    elem_mult(conj(rs_curr.ce(6,10)),rs_curr.ce(7,11))
  );
  // Calculate the noise on each TOE estimate.
  vec toe_np=rs_curr_np*rs_curr_np+2*rs_curr_np*sqr(concat(rs_curr_filt(0,4),rs_curr_filt(6,10)));
  // Calculate the weight to use for each estimate
  vec weight=elem_div(sqr(concat(rs_curr_filt(0,4),rs_curr_filt(6,10))),toe_np);
  // MRC
  complex <double> toe_comb=sum(elem_mult(toe,to_cvec(weight)));
  double toe_comb_np=sum(elem_mult(toe_np,weight,weight));
  // Scale. Only necessary for NP.
  double scale=1.0/sum(elem_mult(sqr(concat(rs_curr_filt(0,4),rs_curr_filt(6,10))),weight));
  toe_comb=toe_comb*scale;
  toe_comb_np=toe_comb_np*scale*scale;
  double delay=-arg(toe_comb)/6/(2*pi/128);
  double delay_np=MAX(toe_comb_np/2,.001);
  //cout << delay << " " << db10(delay_np) << endl;
  //delay=0;
  //cout << delay_np << endl;

  // Update frame timing based on TOE
  {
    boost::mutex::scoped_lock lock(tracked_cell.mutex);
    double diff=WRAP((rs_curr.frame_timing+delay)-tracked_cell.frame_timing,-19200.0/2,19200.0/2);
    diff=(0*(1/.0001)+diff*(1/delay_np))/(1/.0001+1/delay_np);
    tracked_cell.frame_timing=itpp_ext::matlab_mod(tracked_cell.frame_timing+diff,19200.0);
    //cout << "TO: " << setprecision(15) << tracked_cell.frame_timing << endl;
  }
}

void do_fd_ac(
  tracked_cell_t & tracked_cell,
  const ce_raw_fifo_pdu_t & rs_curr,
  const double & rs_curr_sp,
  const double & rs_curr_np
) {
  cvec ac_fd(12);
  ac_fd=complex <double> (0,0);
  for (uint8 d=0;d<12;d++) {
    for (uint8 t=0;t<12-d;t++) {
      ac_fd(d)+=conj(rs_curr.ce(t))*rs_curr.ce(t+d);
    }
    ac_fd(d)=ac_fd(d)/(12-d);
  }
  // Normalize
  //ac_fd=ac_fd/ac_fd(0);
  ac_fd=ac_fd/rs_curr_sp;
  vec ac_fd_np=(rs_curr_np*rs_curr_np/(rs_curr_sp*rs_curr_sp)+2*rs_curr_np/rs_curr_sp)/itpp_ext::matlab_range(12.0,-1.0,1.0);
  {
    boost::mutex::scoped_lock lock(tracked_cell.mutex);
    tracked_cell.ac_fd=elem_div(tracked_cell.ac_fd*(1/.001)+elem_mult(ac_fd,to_cvec(1.0/ac_fd_np)),to_cvec(1/.001+1.0/ac_fd_np));
  }
}

void interp2d(
  const tracked_cell_t & tracked_cell,
  const ce_filt_fifo_pdu_t & rs_prev,
  const ce_filt_fifo_pdu_t & rs_curr,
  const uint8 & port_num,
  deque <ce_interp_fifo_pdu_t> & ce_interp_fifo,
  uint8 & ce_interp_fifo_initialized
) {
  // Interpolate in the frequency domain.
  vec X=itpp_ext::matlab_range(rs_prev.shift,6.0,71.0);
  cvec Y=rs_prev.ce_filt;
  vec x=itpp_ext::matlab_range(0.0,71.0);
  cvec rs_prev_interp=interp1(X,Y,x);
  X=itpp_ext::matlab_range(rs_curr.shift,6.0,71.0);
  Y=rs_curr.ce_filt;
  cvec rs_curr_interp=interp1(X,Y,x);

  // Interpolate in the time domain and push onto FIFO
  uint8 slot_num=rs_prev.slot_num;
  uint8 sym_num=rs_prev.sym_num;
  // Time difference between the current and previous channel estimates.
  double time_diff;
  if (port_num>2) {
    time_diff=0.0005;
  } else {
    if (tracked_cell.cp_type==cp_type_t::EXTENDED) {
      time_diff=3*(128+32);
    } else {
      if (rs_prev.sym_num==0) {
        time_diff=4*(128+9);
      } else {
        time_diff=2*(128+9)+(128+10);
      }
    }
    time_diff=time_diff*(1/(FS_LTE/16));
  }

  double time_offset=0;
  while ((slot_num!=rs_curr.slot_num)||(sym_num!=rs_curr.sym_num)) {
    // Interpolate in the time domain.
    cvec rs_mid=rs_prev_interp+(rs_curr_interp-rs_prev_interp)*(time_offset/time_diff);
    double rs_mid_sp=rs_prev.sp+(rs_curr.sp-rs_prev.sp)*(time_offset/time_diff);
    double rs_mid_np=rs_prev.np+(rs_curr.np-rs_prev.np)*(time_offset/time_diff);

    // Push onto the interpolated CE fifo.
    ce_interp_fifo_pdu_t pdu;
    pdu.ce_interp=rs_mid;
    pdu.sp=rs_mid_sp;
    pdu.np=rs_mid_np;
    if (!ce_interp_fifo_initialized) {
      ce_interp_fifo_initialized=true;
      uint8 tsy=0;
      uint8 tsl=0;
      while ((tsy!=sym_num)||(tsl!=slot_num)) {
        pdu.sym_num=tsy;
        pdu.slot_num=tsl;
        ce_interp_fifo.push_back(pdu);
        tsy=mod(tsy+1,tracked_cell.n_symb_dl());
        if (tsy==0) {
          tsl=mod(tsl+1,20);
        }
      }
    }
    pdu.slot_num=slot_num;
    pdu.sym_num=sym_num;
    ce_interp_fifo.push_back(pdu);

    // Increment counters.
    if (tracked_cell.cp_type==cp_type_t::EXTENDED) {
      time_offset+=(128+32)*(1/(FS_LTE/16));
    } else {
      if (sym_num==6) {
        time_offset+=(128+10)*(1/(FS_LTE/16));
      } else {
        time_offset+=(128+9)*(1/(FS_LTE/16));
      }
    }
    slot_sym_inc(tracked_cell.n_symb_dl(),slot_num,sym_num);
    //sym_num=mod(sym_num+1,tracked_cell.n_symb_dl());
    //if (sym_num==0) {
    //  slot_num=mod(slot_num+1,20);
    //}
  }
}

// Small helper function returns true if all the fifos contain data.
bool ce_ready(
  const vector <deque <ce_interp_fifo_pdu_t> > & ce_interp_fifo
) {
  bool ready=true;
  for (uint8 t=0;t<ce_interp_fifo.size();t++) {
    ready=ready&&(!ce_interp_fifo[t].empty());
    if (!ready)
      break;
  }
  return ready;
}

// Similar to pbch_extract but this one works only on MIB samples. The
// other function has access to all the OFDM symbols.
void pbch_extract_rt(
  const tracked_cell_t & tracked_cell,
  const deque <mib_fifo_pdu_t> & mib_fifo,
  cvec & pbch_sym,
  cmat & pbch_ce,
  mat & np
) {
  // Shortcuts
  const int16 & n_id_cell=tracked_cell.n_id_cell;
  const uint8 & n_ports=tracked_cell.n_ports;
  const cp_type_t::cp_type_t & cp_type=tracked_cell.cp_type;

  uint16 n_syms=(cp_type==cp_type_t::NORMAL)?(1920/2):(1728/2);
  pbch_sym.set_size(n_syms);
  pbch_ce.set_size(n_ports,n_syms);
  np.set_size(n_ports,n_syms);

  const uint8 v_shift_m3=mod(n_id_cell,3);
  uint16 idx=0;
  for (uint8 fr=0;fr<4;fr++) {
    for (uint8 symn=0;symn<4;symn++) {
      for (uint16 sc=0;sc<72;sc++) {
        // Skip if there might be an RS occupying this position.
        if ((mod(sc,3)==v_shift_m3)&&((symn==0)||(symn==1)||((symn==3)&&(cp_type==cp_type_t::EXTENDED)))) {
          continue;
        }
        pbch_sym(idx)=mib_fifo[fr*4+symn].syms(sc);
        pbch_ce.set_col(idx,mib_fifo[fr*4+symn].ce.get_col(sc));
        np.set_col(idx,mib_fifo[fr*4+symn].np);
        //cout << mib_fifo[fr*4+symn].np << endl;
        idx++;
      }
    }
  }
  ASSERT(idx==n_syms);
}

int8 do_mib_decode(
  tracked_cell_t & tracked_cell,
  const cvec & syms,
  const cmat & ce,
  const vec & sp,
  const vec & np,
  const uint8 & data_slot_num,
  const uint8 & data_sym_num,
  deque <mib_fifo_pdu_t> & mib_fifo,
  bool & mib_fifo_synchronized,
  double & mib_fifo_decode_failures
) {
  // Assemble symbols for MIB decoding.
  if ((data_slot_num==1)&&(data_sym_num<=3)) {
    mib_fifo_pdu_t pdu;
    pdu.syms=syms;
    pdu.ce=ce;
    pdu.sp=sp;
    pdu.np=np;
    mib_fifo.push_back(pdu);
  }

  // Does the MIB fifo have enough data to attempt MIB decoding?
  //cout << mib_fifo.size() << endl;
  if (mib_fifo.size()==16) {
    //cout << "MIB decode starting" << endl;
    cvec pbch_sym;
    cmat pbch_ce;
    mat np_pre;
    pbch_extract_rt(tracked_cell,mib_fifo,pbch_sym,pbch_ce,np_pre);
    //cout << abs(pbch_sym) << endl;
    //cout << arg(pbch_sym) << endl;
    //exit(-1);
    //cout << pbch_sym << endl;
    //cout << "pbch_sym: " << pbch_sym << endl;
    //cout << "pbch_ce: " << pbch_ce << endl;

    vec np_mib;
    cvec syms_mib;
    // Much of the following code is copied from the searcher...
    // TODO: wrap copied code into a function.
    // Perform channel compensation and also estimate noise power in each
    // symbol.
    if (tracked_cell.n_ports==1) {
      cvec gain=conj(elem_div(pbch_ce.get_row(0),to_cvec(sqr(pbch_ce.get_row(0)))));
      syms_mib=elem_mult(pbch_sym,gain);
      np_mib=elem_mult(np_pre.get_row(0),sqr(gain));
      //cout << db10(np_mib) << endl;
      //cout << np_mib << endl;
      //cout << abs(gain) << endl;
      //cout << np_pre.get_row(0) << endl;
      //cout << (np_pre.get_row(0)<0) << endl;
      //exit(-1);
    } else {
      syms_mib.set_size(length(pbch_sym));
      np_mib.set_size(length(pbch_sym));
#ifndef NDEBUG
      syms_mib=NAN;
      np_mib=NAN;
#endif
      for (int32 t=0;t<length(syms_mib);t+=2) {
        // Simple zero-forcing
        // http://en.wikipedia.org/wiki/Space-time_block_coding_based_transmit_diversity
        complex <double> h1,h2;
        double np_temp;
        if (tracked_cell.n_ports==2) {
          h1=(pbch_ce(0,t)+pbch_ce(0,t+1))/2;
          h2=(pbch_ce(1,t)+pbch_ce(1,t+1))/2;
          np_temp=(np_pre(0,t)+np_pre(1,t))/2;
        } else {
          if (mod(t,4)==0) {
            h1=(pbch_ce(0,t)+pbch_ce(0,t+1))/2;
            h2=(pbch_ce(2,t)+pbch_ce(2,t+1))/2;
            np_temp=(np_pre(0,t)+np_pre(2,t))/2;
          } else {
            h1=(pbch_ce(1,t)+pbch_ce(1,t+1))/2;
            h2=(pbch_ce(3,t)+pbch_ce(3,t+1))/2;
            np_temp=(np_pre(1,t)+np_pre(3,t))/2;
          }
        }
        complex <double> x1=pbch_sym(t);
        complex <double> x2=pbch_sym(t+1);
        double scale=pow(h1.real(),2)+pow(h1.imag(),2)+pow(h2.real(),2)+pow(h2.imag(),2);
        syms_mib(t)=(conj(h1)*x1+h2*conj(x2))/scale;
        syms_mib(t+1)=conj((-conj(h2)*x1+h1*conj(x2))/scale);
        np_mib(t)=(pow(abs(h1)/scale,2)+pow(abs(h2)/scale,2))*np_temp;
        np_mib(t+1)=np_mib(t);
      }
      // 3dB factor comes from precoding for transmit diversity
      syms_mib*=pow(2,0.5);
    }
    //cout << abs(syms_mib) << endl;
    //cout << arg(syms_mib) << endl;
    //cout << db10(np_mib) << endl;
    //exit(-1);
    //cout << "syms_mib: " << syms_mib << endl;

    // Extract the bits from the complex modulated symbols.
    vec e_est=lte_demodulate(syms_mib,np_mib,modulation_t::QAM);
    // Unscramble
    bvec scr=lte_pn(tracked_cell.n_id_cell,length(e_est));
    for (int32 t=0;t<length(e_est);t++) {
      if (scr(t)) e_est(t)=-e_est(t);
    }
    // Undo ratematching
    mat d_est=lte_conv_deratematch(e_est,40);
    // Decode
    bvec c_est=lte_conv_decode(d_est);
    // Calculate received CRC
    bvec crc_est=lte_calc_crc(c_est(0,23),CRC16);
    // Apply CRC mask
    if (tracked_cell.n_ports==2) {
      for (uint8 t=0;t<16;t++) {
        crc_est(t)=1-((int)crc_est(t));
      }
    } else if (tracked_cell.n_ports==4) {
      for (uint8 t=1;t<length(crc_est);t+=2) {
        crc_est(t)=1-((int)crc_est(t));
      }
    }
    // Did we find it?
    if (crc_est==c_est(24,-1)) {
      // YES!
      cout << "Cell ID " << tracked_cell.n_id_cell << " MIB SUCCESS!" << endl;
      mib_fifo_synchronized=1;
      mib_fifo_decode_failures=0;
      for (uint8 t=0;t<16;t++) {
        mib_fifo.pop_front();
      }
    } else {
      // No :(
      cout << "Cell ID " << tracked_cell.n_id_cell << " MIB failure!" << endl;
      if (mib_fifo_synchronized) {
        mib_fifo_decode_failures++;
        for (uint8 t=0;t<16;t++) {
          mib_fifo.pop_front();
        }
      } else {
        mib_fifo_decode_failures+=0.25;
        for (uint8 t=0;t<4;t++) {
          mib_fifo.pop_front();
        }
      }
    }

    // 10ms of time increases mib_fifo_decode_failures by 0.25.
    // After several seconds of MIB decoding failures, drop the cell.
    if (mib_fifo_decode_failures>=100) {
      cout << "Dropped a cell!" << endl;
      boost::mutex::scoped_lock lock(tracked_cell.mutex);
      tracked_cell.kill_me=true;
      return -1;
    }
    //cout << "MIB decode finished" << endl;
  }

  return 0;
}

// Process that tracks a cell that has been found by the searcher.
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
  vector <deque <ce_interp_fifo_pdu_t> > ce_interp_fifo(tracked_cell.n_ports);
  // Cannot use bool here because all the bits would be packed into bytes
  // and could not be passed by reference.
  vector <uint8> ce_interp_fifo_initialized(tracked_cell.n_ports,0);
  deque <mib_fifo_pdu_t> mib_fifo;
  bool mib_fifo_synchronized=false;
  double mib_fifo_decode_failures=0;
  // Now that everything has been initialized, indicate to the producer
  // thread that we are ready for data.
  {
    boost::mutex::scoped_lock lock(tracked_cell.mutex);
    tracked_cell.tracker_proc_ready=true;
  }
  // Each iteration of this loop processes one OFDM symbol.
  while (true) {
    // Get the next frequency domain sample from the fifo.
    cvec syms;
    double frequency_offset;
    double frame_timing;
    get_fd(tracked_cell,global_thread_data.fc,slot_num,sym_num,cn,bulk_phase_offset,syms,frequency_offset,frame_timing);

    // Save this information into the data fifo for further processing
    // once the channel estimates are ready for this ofdm symbol.
    data_fifo_pdu_t dfp;
    dfp.slot_num=slot_num;
    dfp.sym_num=sym_num;
    dfp.syms=syms;
    data_fifo.push_back(dfp);

    // Extract any RS that might be present.
    for (uint8 port_num=0;port_num<tracked_cell.n_ports;port_num++) {
      double shift=rs_dl.get_shift(slot_num,sym_num,port_num);
      if (isnan(shift))
        continue;
      //cout << "S" << shift << endl;
      cvec rs_raw=syms(itpp_ext::matlab_range(round_i(shift),6,71));
      //cout << "A" << rs_dl.get_rs(slot_num,sym_num) << endl;
      //cout << slot_num << " x " << sym_num << endl;
      //cout << "B" << rs_raw << endl;
      cvec ce_raw=elem_mult(rs_raw,conj(rs_dl.get_rs(slot_num,sym_num)));
      ce_raw_fifo_pdu_t cerp;
      cerp.shift=shift;
      cerp.slot_num=slot_num;
      cerp.sym_num=sym_num;
      cerp.ce=ce_raw;
      cerp.frequency_offset=frequency_offset;
      cerp.frame_timing=frame_timing;
      ce_raw_fifo[port_num].push_back(cerp);
    }

    // For each port, filter and perform FOE, TOE, and interpolation on the raw
    // channel estimates. Also perform some measurements.
    // All tasks that need access to the raw channel estimates should
    // go in this loop.
    for (uint8 port_num=0;port_num<tracked_cell.n_ports;port_num++) {
      // In order to filter the raw channel estimates for OFDM symbol n,
      // we need the raw channel estimates for OFDM symbols n-1, n, and n+1.
      if (ce_raw_fifo[port_num].size()!=3)
        continue;

      // Shortcuts
      ce_raw_fifo_pdu_t & rs_prev=ce_raw_fifo[port_num][0];
      ce_raw_fifo_pdu_t & rs_curr=ce_raw_fifo[port_num][1];
      ce_raw_fifo_pdu_t & rs_next=ce_raw_fifo[port_num][2];

      // Perform primitive filtering by averaging nearby samples.
      const cvec rs_curr_filt=filter_ce(rs_prev,rs_curr,rs_next);
      const double rs_curr_sp=sigpower(rs_curr_filt);
      const double rs_curr_np=sigpower(rs_curr.ce-rs_curr_filt);
      //cout << "SP RSC " << db10(sigpower(rs_curr.ce)) << endl;
      //cout << "SP/NP " << db10(rs_curr_sp) << " / " << db10(rs_curr_np) << endl;
      // Store filtered channel estimates.
      ce_filt_fifo_pdu_t pdu;
      pdu.shift=rs_curr.shift;
      pdu.slot_num=rs_curr.slot_num;
      pdu.sym_num=rs_curr.sym_num;
      pdu.sp=rs_curr_sp;
      pdu.np=rs_curr_np;
      //cout << "Cell " << tracked_cell.n_id_cell << " port" << port_num << " CRS SNR " << db10(rs_curr_sp/rs_curr_np) << " dB" << endl;
      pdu.ce_filt=rs_curr_filt;
      //pdu.frequency_offset=rs_curr.frequency_offset;
      //pdu.frame_timing=rs_curr.frame_timing;
      ce_filt_fifo[port_num].push_back(pdu);

      // FOE
      do_foe(global_thread_data,rs_prev,rs_next,rs_curr_np,rs_curr_filt);

      // TOE
      //do_toe(tracked_cell,rs_curr,rs_curr_filt,rs_curr_np);
      do_toe_v2(tracked_cell,rs_prev,rs_curr,rs_curr_sp,rs_curr_np);

      // Estimate frequency domain autocorrelations.
      do_fd_ac(tracked_cell,rs_curr,rs_curr_sp,rs_curr_np);

      // Estimate the time domain autocorrelations.
      //do_td_ac();

      // Finished working with the raw channel estimates.
      ce_raw_fifo[port_num].pop_front();
    }

    // Tasks that need access to the filtered channel estimates should
    // go in this loop.
    for (uint8 port_num=0;port_num<tracked_cell.n_ports;port_num++) {
      // For interpolation, we need two OFDM symbols.
      if (ce_filt_fifo[port_num].size()!=2)
        continue;

      // Shortcuts
      ce_filt_fifo_pdu_t & rs_prev=ce_filt_fifo[port_num][0];
      ce_filt_fifo_pdu_t & rs_curr=ce_filt_fifo[port_num][1];

      interp2d(tracked_cell,rs_prev,rs_curr,port_num,ce_interp_fifo[port_num],ce_interp_fifo_initialized[port_num]);

      // Finished working with the filtered channel estimates.
      ce_filt_fifo[port_num].pop_front();
    }

    // Process data if channel estimates are available for each antenna and for
    // every data sample.
    while ((!data_fifo.empty())&&ce_ready(ce_interp_fifo)) {
      // Synchronization check.
      for (uint8 t=0;t<tracked_cell.n_ports;t++) {
        if (
          (data_fifo.front().slot_num!=ce_interp_fifo[t].front().slot_num)||
          (data_fifo.front().sym_num!=ce_interp_fifo[t].front().sym_num)
        ) {
          cerr << "Error: synchronization error! Check code!" << endl;
          exit(-1);
        }
      }

      // Shortcuts
      cvec & syms=data_fifo.front().syms;
      cmat ce(tracked_cell.n_ports,72);
      vec sp(tracked_cell.n_ports);
      vec np(tracked_cell.n_ports);
      uint8 data_slot_num=data_fifo.front().slot_num;
      uint8 data_sym_num=data_fifo.front().sym_num;
      for (uint8 t=0;t<tracked_cell.n_ports;t++) {
        ce.set_row(t,ce_interp_fifo[t].front().ce_interp);
        sp(t)=ce_interp_fifo[t].front().sp;
        np(t)=ce_interp_fifo[t].front().np;
      }

      // Measure signal power and noise power on PSS/SSS
      //do_pss_sss_sigpwer();

      // Perform MIB decoding
      if (do_mib_decode(tracked_cell,syms,ce,sp,np,data_slot_num,data_sym_num,mib_fifo,mib_fifo_synchronized,mib_fifo_decode_failures)==-1) {
        // We have failed to detect an MIB for a long time. Exit this
        // thread.
        cout << "Tracker thread exiting..." << endl;
        return;
      }

      // Done processing data. Pop data vector and CE vectors.
      data_fifo.pop_front();
      for (uint8 t=0;t<tracked_cell.n_ports;t++) {
        ce_interp_fifo[t].pop_front();
      }
    }

    // Increase the local counter
    slot_sym_inc(tracked_cell.n_symb_dl(),slot_num,sym_num);
    //sym_num=mod(sym_num+1,tracked_cell.n_symb_dl());
    //if (sym_num==0) {
    //  slot_num=mod(slot_num+1,20);
    //}
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

  // Shortcut
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
    peak_search(xc_incoherent_collapsed_pow,xc_incoherent_collapsed_frq,Z_th1,f_search_set,fc,xc_incoherent_single,DS_COMB_ARM,detected_cells);

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
        cout << "  frame start: " << (*iterator).frame_start << endl;
      }

      // Launch a cell tracker process!
      k_factor=k_factor;
      //cout << "Timing error is purposely introduced here!!!" << endl;
      //tracked_cell_t * new_cell = new tracked_cell_t((*iterator).n_id_cell(),(*iterator).n_ports,(*iterator).cp_type,(*iterator).frame_start/k_factor+capbuf_sync.late+global_1);
      tracked_cell_t * new_cell = new tracked_cell_t((*iterator).n_id_cell(),(*iterator).n_ports,(*iterator).cp_type,(*iterator).frame_start/k_factor+capbuf_sync.late);
      (*new_cell).thread=boost::thread(tracker_proc,boost::ref(*new_cell),boost::ref(global_thread_data));
      {
        boost::mutex::scoped_lock lock(tracked_cell_list.mutex);
        tracked_cell_list.tracked_cells.push_back(new_cell);
      }
      //cout << "Only one cell is allowed to be detected!!!" << endl;
      //sleep(1000000);

      ++iterator;
    }
  }
  // Will never reach here...
}

// Process that takes samples and distributes them to the appropriate
// process.
void producer_proc(
  sampbuf_sync_t & sampbuf_sync,
  capbuf_sync_t & capbuf_sync,
  global_thread_data_t & global_thread_data,
  tracked_cell_list_t & tracked_cell_list,
  double & fc
) {
  // Main loop which distributes data to the appropriate subthread.
  //Real_Timer tt;
  double sample_time=-1;
  bool searcher_capbuf_filling=false;
  uint32 searcher_capbuf_idx=0;
  uint32 n_samps_read=0;
  unsigned long long int sample_number=0;
  // Elevate privileges of the producer thread.
  //int retval=nice(-10);
  //if (retval==-1) {
  //  cerr << "Error: could not set elevated privileges" << endl;
  //  exit(-1);
  //}
  //int policy=SCHED_RR;
  //struct sched_param param;
  //param.sched_priority=50;
  ////pthread_getschedparam(pthread_self(), &policy, &param);
  //if (pthread_setschedparam(pthread_self(),policy,&param)) {
  //  cerr << "Error: could not elevate main thread priority" << endl;
  //  exit(-1);
  //}
  //tt.tic();
  while (true) {
    // Each iteration of this loop processes one sample.
    double k_factor;
    double frequency_offset;
    {
      boost::mutex::scoped_lock lock(global_thread_data.frequency_offset_mutex);
      frequency_offset=global_thread_data.frequency_offset;
      k_factor=(fc-frequency_offset)/fc;
      if (mod(sample_number++,19200*100)==0) {
        //tt.toc_print();
        //cout << sample_number << endl;
        // Status message displayed every second.
        cout << "System frequency offset is currently: " << frequency_offset << endl;
        boost::mutex::scoped_lock lock(tracked_cell_list.mutex);
        list <tracked_cell_t *>::iterator it=tracked_cell_list.tracked_cells.begin();
        while (it!=tracked_cell_list.tracked_cells.end()) {
          cout << "Cell ID " << (*(*it)).n_id_cell << " TO: " << setprecision(10) << (*(*it)).frame_timing << endl;
          ++it;
        }
      }
    }

    // Get the next sample
    complex <double> sample;
    {
      boost::mutex::scoped_lock lock(sampbuf_sync.mutex);
      while (sampbuf_sync.fifo.size()<2) {
        sampbuf_sync.condition.wait(lock);
      }
      uint8 real=sampbuf_sync.fifo[0];
      uint8 imag=sampbuf_sync.fifo[1];
      sample=complex <double>((real-127.0)/128.0,(imag-127.0)/128.0);
      sampbuf_sync.fifo.pop_front();
      sampbuf_sync.fifo.pop_front();
    }
    n_samps_read++;
    sample_time+=1.0/k_factor;
    sample_time=WRAP(sample_time,0.0,19200.0);

    // Should we begin filling the capture buffer?
    {
      boost::mutex::scoped_lock lock(capbuf_sync.mutex);
      if ((capbuf_sync.request)&&(!searcher_capbuf_filling)&&(abs(WRAP(sample_time-0,-19200.0/2,19200.0/2))<0.5)) {
        //cout << "searcher data cap beginning" << endl;
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
        //cout << "searcher data cap finished" << endl;
      }
    }

    // Loop for each tracked cell and save data, if necessary. Also delete
    // threads that may have lost lock.
    {
      boost::mutex::scoped_lock lock(tracked_cell_list.mutex);
      list <tracked_cell_t *>::iterator it=tracked_cell_list.tracked_cells.begin();
      while (it!=tracked_cell_list.tracked_cells.end()) {
        tracked_cell_t & tracked_cell=(*(*it));
        boost::mutex::scoped_lock lock2(tracked_cell.mutex);

        // Delete the tracker if lock has been lost.
        if (tracked_cell.kill_me) {
          // TODO: Fix memory leak here!
          //delete (*it);
          it=tracked_cell_list.tracked_cells.erase(it);
          continue;
        }

        // See if we should start filling the buffer.
        if (tracked_cell.tracker_proc_ready&&!tracked_cell.filling) {
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
            td_fifo_pdu_t p;
            p.data=tracked_cell.buffer;
            p.slot_num=tracked_cell.slot_num;
            p.sym_num=tracked_cell.sym_num;
            p.late=tracked_cell.late;
            // Record the frequency offset and frame timing as they were
            // during the capture.
            p.frequency_offset=frequency_offset;
            p.frame_timing=tracked_cell.frame_timing;
            tracked_cell.fifo.push(p);
            tracked_cell.condition.notify_one();
            //cout << "Sleeping..." << endl;
            sleep(0.0005/8);
            //cout << "fifo size: " << tracked_cell.fifo.size() << endl;
            // Calculate trigger parameters of next capture
            tracked_cell.filling=false;
            if (tracked_cell.cp_type==cp_type_t::EXTENDED) {
              tracked_cell.target_cap_start_time+=32+128;
            } else {
              tracked_cell.target_cap_start_time+=(tracked_cell.sym_num==6)?128+10:128+9;
            }
            tracked_cell.target_cap_start_time=mod(tracked_cell.target_cap_start_time,19200);
            slot_sym_inc(tracked_cell.n_symb_dl(),tracked_cell.slot_num,tracked_cell.sym_num);
          }
        }
        ++it;
      }
    }
  }
}

static void rtlsdr_callback(
  unsigned char * buf,
  uint32_t len,
  void * ctx
) {
  sampbuf_sync_t & sampbuf_sync=*((sampbuf_sync_t *)ctx);

  //cout << "Callback with " << len << " samples" << endl;

  if (len==0) {
    cerr << "Received 'zero' samples from USB..." << endl;
    exit(-1);
  }

  boost::mutex::scoped_lock lock(sampbuf_sync.mutex);
  for (uint32 t=0;t<len;t++) {
    sampbuf_sync.fifo.push_back(buf[t]);
  }
  sampbuf_sync.condition.notify_one();
}

// Main routine.
int main(
  const int argc,
  char * const argv[]
) {
  // This is so that CTRL-C properly closes the rtl-sdr device before exiting
  // the program.
  //struct sigaction sigact;
  //sigact.sa_handler=sighandler;
  //sigemptyset(&sigact.sa_mask);
  //sigact.sa_flags=0;
  //sigaction(SIGINT,&sigact,NULL);
  //sigaction(SIGTERM,&sigact,NULL);
  //sigaction(SIGQUIT,&sigact,NULL);
  //sigaction(SIGPIPE,&sigact,NULL);

  // Command line parameters are stored here.
  double fc;
  double ppm;
  double correction;
  int32 device_index;
  bool use_recorded_data;
  string filename;
  bool repeat;
  bool rtl_sdr_format;
  double noise_power;

  // Get search parameters from the user
  parse_commandline(argc,argv,fc,ppm,correction,device_index,use_recorded_data,filename,repeat,rtl_sdr_format,noise_power);

  // Open the USB device.
  if (!use_recorded_data)
    config_usb(correction,device_index,fc);

  // Data shared between threads
  sampbuf_sync_t sampbuf_sync;
  tracked_cell_list_t tracked_cell_list;
  capbuf_sync_t capbuf_sync;
  global_thread_data_t global_thread_data;

  // Calibrate the dongle's oscillator. This is similar to running the
  // program CellSearch with only one center frequency. All information
  // is discarded except for the frequency offset.
  global_thread_data.frequency_offset=kalibrate(fc,ppm,correction,use_recorded_data,filename,rtl_sdr_format,noise_power);

  // Start the cell searcher thread.
  // Now that the oscillator has been calibrated, we can perform
  // a 'real' search.
  capbuf_sync.request=false;
  capbuf_sync.capbuf.set_size(19200*8);
  global_thread_data.fc=fc;
  boost::thread searcher_thread(searcher_proc,boost::ref(capbuf_sync),boost::ref(global_thread_data),boost::ref(tracked_cell_list));

  // Wrap the USB device to simplify obtaining complex samples one by one.
  rtl_wrap sample_source(dev,use_recorded_data,filename,repeat,rtl_sdr_format,noise_power);

  // Start the producer thread.
  boost::thread producer_thread(producer_proc,boost::ref(sampbuf_sync),boost::ref(capbuf_sync),boost::ref(global_thread_data),boost::ref(tracked_cell_list),boost::ref(fc));

  // Start the async read process. This should never return.
  rtlsdr_read_async(dev,rtlsdr_callback,(void *)&sampbuf_sync,0,0);

  // Successful exit.
  exit (0);
}

