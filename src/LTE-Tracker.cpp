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
#include <sys/stat.h>
#include <sys/syscall.h>
#include <sys/types.h>
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
#include "LTE-Tracker.h"

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
    // sync mode is unreliable in that it may drop samples. However, in
    // this case we can still use it because it will gurantee a minimum
    // number of samples have been discarded.
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
    //cvec sig_tx_pre;
    itpp_ext::rtl_sdr_to_cvec(filename,sig_tx);
    /*
    sig_tx.set_size(length(sig_tx_pre)*2-10);
    for (int32 t=0;t<length(sig_tx);t+=2) {
      // FIXME: Do proper interpolation
      sig_tx(t)=sig_tx_pre(t>>1);
      sig_tx(t+1)=(sig_tx_pre(t>>1)+sig_tx_pre((t>>1)+1))/2;
    }
    // Drop several seconds while AGC converges.
    */
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

static void rtlsdr_callback(
  unsigned char * buf,
  uint32_t len,
  void * ctx
) {
  sampbuf_sync_t & sampbuf_sync=*((sampbuf_sync_t *)ctx);

  //cout << "Callback with " << len << " samples" << endl;

  if (len==0) {
    cerr << "Error: received no samples from USB device..." << endl;
    exit(-1);
  }

  boost::mutex::scoped_lock lock(sampbuf_sync.mutex);
  for (uint32 t=0;t<len;t++) {
    sampbuf_sync.fifo.push_back(buf[t]);
  }
  sampbuf_sync.fifo_peak_size=MAX(sampbuf_sync.fifo.size(),sampbuf_sync.fifo_peak_size);
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
  global_thread_data_t global_thread_data(fc);
  global_thread_data.main_thread_id=syscall(SYS_gettid);

  // Calibrate the dongle's oscillator. This is similar to running the
  // program CellSearch with only one center frequency. All information
  // is discarded except for the frequency offset.
  global_thread_data.frequency_offset(kalibrate(fc,ppm,correction,use_recorded_data,filename,rtl_sdr_format,noise_power));

  // Start the cell searcher thread.
  // Now that the oscillator has been calibrated, we can perform
  // a 'real' search.
  capbuf_sync.request=false;
  capbuf_sync.capbuf.set_size(19200*8);
  //global_thread_data.fc=fc;
  boost::thread searcher_thr(searcher_thread,boost::ref(capbuf_sync),boost::ref(global_thread_data),boost::ref(tracked_cell_list));

  // Wrap the USB device to simplify obtaining complex samples one by one.
  rtl_wrap sample_source(dev,use_recorded_data,filename,repeat,rtl_sdr_format,noise_power);

  // Start the producer thread.
  boost::thread producer_thr(producer_thread,boost::ref(sampbuf_sync),boost::ref(capbuf_sync),boost::ref(global_thread_data),boost::ref(tracked_cell_list),boost::ref(fc));

  sampbuf_sync.fifo_peak_size=0;

  // Launch the display thread
  boost::thread display_thr(display_thread,boost::ref(sampbuf_sync),boost::ref(global_thread_data),boost::ref(tracked_cell_list));

  if (use_recorded_data) {
    if (!rtl_sdr_format) {
      cerr << "Error: only rtl_sdr format supported currently." << endl;
      exit(-1);
    }
    // Get filesize

    struct stat filestatus;
    stat(filename.c_str(),&filestatus);
    uint32 file_length=filestatus.st_size;
    //cout << "file length: " << file_length << " bytes\n";
    if (floor(file_length/2.0)!=file_length/2.0) {
      cout << "Warning: file contains an odd number of samples" << endl;
    }

    // Open file
    FILE *file;
    file=fopen(filename.c_str(),"rb");
    if (!file) {
      cerr << "Error: could not open input file" << endl;
      exit(-1);
    }

    // Read entire file, all at once!
    uint8 * buffer=(uint8 *)malloc(file_length*sizeof(uint8));
    uint32 n_read=fread(buffer,1,file_length,file);
    if (n_read!=file_length) {
      cerr << "Error: error while reading file" << endl;
      exit(-1);
    }

    uint32 offset=0;
    while (true) {
      {
        boost::mutex::scoped_lock lock(sampbuf_sync.mutex);
        for (uint32 t=offset;t<file_length;t++) {
          sampbuf_sync.fifo.push_back(buffer[offset++]);
          if (mod(offset,1920000*2)==0)
            break;
        }
        sampbuf_sync.condition.notify_one();
      }
      sleep(1);
      if (offset==file_length)
        break;
    }
    free(buffer);

    //sleep(2<<30);
    sleep(500*60);
    exit(-1);
  } else {
    // Start the async read process. This should never return.
    rtlsdr_read_async(dev,rtlsdr_callback,(void *)&sampbuf_sync,0,0);
  }

  // Successful exit. (Should never get here!)
  exit (0);
}

