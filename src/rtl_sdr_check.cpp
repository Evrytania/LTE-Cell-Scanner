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
#include <iostream>
#include <itpp/itbase.h>
#include <boost/math/special_functions/gamma.hpp>
#include <curses.h>
#include "rtl-sdr.h"
#include "common.h"
#include "lte_lib.h"
#include "itpp_ext.h"
#include "constants.h"
#include "macros.h"
#include "dsp.h"

using namespace std;
using namespace itpp;

double global_1=0;
double global_2=0;
double global_3=0;
double global_4=0;
double global_5=0;
double global_6=0;
double global_7=0;
double global_8=0;
double global_9=0;

// Simple usage screen.
void print_usage() {
  cout << "rtl_sdr_check v" << MAJOR_VERSION << "." << MINOR_VERSION << "." << PATCH_LEVEL << " (" << BUILD_TYPE << ") help screen" << endl << endl;
  cout << "rtl_sdr_check -d cap_data.dat -f fc -s fs -i cell_id [optional_parameters]" << endl;
  cout << "  -h --help" << endl;
  cout << "    This help screen" << endl;
  cout << "  -d --data-file filename" << endl;
  cout << "    capture filename created by rtl_sdr" << endl;
  cout << "  -f --freq-center fc" << endl;
  cout << "    center frequency provided to rtl_sdr" << endl;
  cout << "  -o --freq-offset fo" << endl;
  cout << "    frequency offset of LTE carrier" << endl;
  cout << "  -s --sampling-freq fs" << endl;
  cout << "    sampling frequency provided to rtl_sdr" << endl;
  cout << "  -i --cell-id cell_id" << endl;
  cout << "    cell id to search for" << endl;
  // Hidden options. Only useful for debugging.
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
}

// Parse the command line arguments and return optional parameters as
// variables.
// Also performs some basic sanity checks on the parameters.
void parse_commandline(
  // Inputs
  const int & argc,
  char * const argv[],
  // Outputs
  string & filename,
  double & fc,
  double & f_off,
  double & fs,
  uint16 & n_id_cell
) {
  // Default values
  bool filename_set=false;
  bool fc_set=false;
  bool f_off_set=false;
  bool fs_set=false;
  bool n_id_cell_set=false;

  while (1) {
    static struct option long_options[] = {
      {"help",          no_argument,       0, 'h'},
      {"data-file",     required_argument, 0, 'd'},
      {"freq-center",   required_argument, 0, 'f'},
      {"freq-offset",   required_argument, 0, 'o'},
      {"sampling-freq", required_argument, 0, 's'},
      {"cell_id",       required_argument, 0, 'i'},
      {"g1",            required_argument, 0, '1'},
      {"g2",            required_argument, 0, '2'},
      {"g3",            required_argument, 0, '3'},
      {"g4",            required_argument, 0, '4'},
      {"g5",            required_argument, 0, '5'},
      {"g6",            required_argument, 0, '6'},
      {"g7",            required_argument, 0, '7'},
      {"g8",            required_argument, 0, '8'},
      {"g9",            required_argument, 0, '9'},
      {0, 0, 0, 0}
    };
    /* getopt_long stores the option index here. */
    int option_index = 0;
    int c = getopt_long (argc, argv, "hd:f:o:s:i:1:2:3:4:5:6:7:8:9:",
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
        ABORT(-1);
        break;
      case 'h':
        print_usage();
        ABORT(-1);
        break;
      case 'd':
        filename=optarg;
        filename_set=true;
        break;
      case 'f':
        fc=strtod(optarg,&endp);
        fc_set=true;
        if ((optarg==endp)||(*endp!='\0')) {
          cerr << "Error: could not parse center frequency" << endl;
          ABORT(-1);
        }
        break;
      case 'o':
        f_off=strtod(optarg,&endp);
        f_off_set=true;
        if ((optarg==endp)||(*endp!='\0')) {
          cerr << "Error: could not parse frequency offset" << endl;
          ABORT(-1);
        }
        break;
      case 's':
        fs=strtod(optarg,&endp);
        fs_set=true;
        if ((optarg==endp)||(*endp!='\0')) {
          cerr << "Error: could not parse sampling frequency" << endl;
          ABORT(-1);
        }
        break;
      case 'i':
        n_id_cell=strtol(optarg,&endp,10);
        n_id_cell_set=true;
        if ((optarg==endp)||(*endp!='\0')) {
          cerr << "Error: could not parse device index" << endl;
          ABORT(-1);
        }
        break;
      case '1':
        global_1=strtod(optarg,&endp);
        if ((optarg==endp)||(*endp!='\0')) {
          cerr << "Error: could not global variable 1" << endl;
          ABORT(-1);
        }
        break;
      case '2':
        global_2=strtod(optarg,&endp);
        if ((optarg==endp)||(*endp!='\0')) {
          cerr << "Error: could not global variable 2" << endl;
          ABORT(-1);
        }
        break;
      case '3':
        global_3=strtod(optarg,&endp);
        if ((optarg==endp)||(*endp!='\0')) {
          cerr << "Error: could not global variable 3" << endl;
          ABORT(-1);
        }
        break;
      case '4':
        global_4=strtod(optarg,&endp);
        if ((optarg==endp)||(*endp!='\0')) {
          cerr << "Error: could not global variable 4" << endl;
          ABORT(-1);
        }
        break;
      case '5':
        global_5=strtod(optarg,&endp);
        if ((optarg==endp)||(*endp!='\0')) {
          cerr << "Error: could not global variable 5" << endl;
          ABORT(-1);
        }
        break;
      case '6':
        global_6=strtod(optarg,&endp);
        if ((optarg==endp)||(*endp!='\0')) {
          cerr << "Error: could not global variable 6" << endl;
          ABORT(-1);
        }
        break;
      case '7':
        global_7=strtod(optarg,&endp);
        if ((optarg==endp)||(*endp!='\0')) {
          cerr << "Error: could not global variable 7" << endl;
          ABORT(-1);
        }
        break;
      case '8':
        global_8=strtod(optarg,&endp);
        if ((optarg==endp)||(*endp!='\0')) {
          cerr << "Error: could not global variable 8" << endl;
          ABORT(-1);
        }
        break;
      case '9':
        global_9=strtod(optarg,&endp);
        if ((optarg==endp)||(*endp!='\0')) {
          cerr << "Error: could not global variable 9" << endl;
          ABORT(-1);
        }
        break;
      case '?':
        /* getopt_long already printed an error message. */
        ABORT(-1);
      default:
        ABORT(-1);
    }
  }

  // Error if extra arguments are found on the command line
  if (optind<argc) {
    cerr << "Error: unknown/extra arguments specified on command line" << endl;
    ABORT(-1);
  }

  // Second order command line checking. Ensure that command line options
  // are consistent.
  if (!filename_set||!fc_set||!fs_set||!n_id_cell_set||!f_off_set) {
    cerr << "Error: required parameter missing (try --help)" << endl;
    ABORT(-1);
  }
  // Frequencies should be on a 100kHz raster.
  if (fc<1e6) {
    cerr << "Error: frequency must be greater than 1MHz" << endl;
    ABORT(-1);
  }
  if (fc/100e3!=itpp::round(fc/100e3)) {
    fc=itpp::round(fc/100e3)*100e3;
    cout << "Warning: frequency has been rounded to the nearest multiple of 100kHz" << endl;
  }
  // Sampling frequency should be 'reasonable'
  if ((fs<0)||(fs>3e6)) {
    cerr << "Error: sampling frequency must be between 0 and 3MHz" << endl;
    ABORT(-1);
  }
  // Frequency offset should be 'reasonable'
  if (abs(f_off)>200e3) {
    cout << "Warning: frequency offset appears to be too large" << endl;
  }

  cout << "rtl_sdr_check v" << MAJOR_VERSION << "." << MINOR_VERSION << "." << PATCH_LEVEL << " (" << BUILD_TYPE << ") beginning" << endl;
  cout << "  data capture filename: " << filename << endl;
  cout << "  center frequency: " << fc/1e6 << " MHz" << endl;
  cout << "  frequency offset: " << f_off << endl;
  cout << "  sampling frequency: " << fs/1e6 << " MHz" << endl;
  cout << "  cell id: " << n_id_cell << endl;
}

// Main routine.
int main(
  const int argc,
  char * const argv[]
) {
  string filename;
  double fc;
  double f_off;
  double fs;
  uint16 n_id_cell;

  // Get options from command line
  parse_commandline(argc,argv,filename,fc,f_off,fs,n_id_cell);
  const uint8 n_id_1=floor(n_id_cell/3);
  const uint8 n_id_2=n_id_cell-3*n_id_1;
  const double k_factor=(fc-f_off)/fc;

  // Read captured data
  cvec cap_data;
  itpp_ext::rtl_sdr_to_cvec(filename,cap_data);
  // Drop first 4 seconds to allow AGC to converge.
  if (length(cap_data)<fs*4) {
    cerr << "Error: not enough captured data" << endl;
    ABORT(-1);
  }
  cap_data=cap_data(fs*4,-1);
  //cout << "Decimated capture data by 2!!!" << endl;
  //cap_data=cap_data(itpp_ext::matlab_range(0,2,length(cap_data)));
  uint32 n_samp=length(cap_data);
  //n_samp=MIN(3e6,n_samp);
  if (length(cap_data)<1e6) {
    cerr << "Error: not enough captured data" << endl;
    ABORT(-1);
  }
  cout << "Will examine " << n_samp << " samples from capture file" << endl;

  // Extract the 128 points of the PSS and SSS.
  PSS_td pss_td;
  SSS_td sss_td;
  cvec pt=pss_td[n_id_2];
  pt=pt(9,-1);
  cvec st=sss_td(n_id_1,n_id_2,0);
  st=st(9,-1);
  //cvec seq_1p92=concat(sss_td(n_id_1,n_id_2,0),pss_td[n_id_2]);
  //MARK;
  //cout << db10(sigpower(pt)) << endl;
  //cout << db10(sigpower(st)) << endl;

  // Interpolate by a large factor and create the interpolated time domain
  // sequence.
#define FACTOR 1024
  cvec pt_interp=interpft(pt,FACTOR*128);
  cvec st_interp=interpft(st,FACTOR*128);
  //MARK;
  //cout << db10(sigpower(pt_interp)) << endl;
  //cout << db10(sigpower(st_interp)) << endl;
  // Add CP's and create concatenated sequence.
  cvec seq_interp=concat(st_interp(119*FACTOR,-1),st_interp,pt_interp(119*FACTOR,-1),pt_interp);
  //vec seq_interp_timestamp=itpp_ext::matlab_range(0,length(seq_interp)-1)*(1/(FS_LTE/16*FACTOR));

  // seq_interp is an ideal sequence sampled at FS_LTE/16*FACTOR. We want
  // a sequence sampled at fs. Choose the samples from seq_interp that are
  // nearest to the desired sampling instant.
  uint32 n_samp_fs=floor((9+128+9+128)*(1/(FS_LTE/16))/(1/(fs*k_factor)));
  vec desired_time=itpp_ext::matlab_range((uint32)0,n_samp_fs-1)*(1/(fs*k_factor));
  ivec best_sample_index=round_i(desired_time/(1/(FS_LTE/16*FACTOR)));
  if (best_sample_index(length(best_sample_index))>=length(seq_interp)) {
    best_sample_index(length(best_sample_index))=length(seq_interp)-1;
  }
  cvec seq=seq_interp(best_sample_index);
  uint32 n_seq=length(seq);

  //MARK;
  //cout << db10(sigpower(seq_1p92(itpp_ext::matlab_range(0,2,length(seq_1p92)-1))-seq)) << endl;
  //exit(-1);

  // Apply frequency shift
  seq=fshift(seq,f_off,fs*k_factor);
  seq=conj(seq)/((double)n_seq);
  //cout << n_seq << endl;

  // Correlate
  vec xc(n_samp-n_seq+1);
  xc=NAN;
  uint32 t;
  cout << "Correlating" << endl;
#ifdef _OPENMP
#pragma omp parallel for shared(seq,cap_data,xc) private(t)
#endif
  for (t=0;t<=n_samp-n_seq;t++) {
    //if (mod(t,1000000)==0) {
    //  cout << t << endl;
    //}
    xc(t)=sqr(elem_mult_sum(cap_data(t,t+n_seq-1),seq));
  }

  //cout << xc << endl;

  // Search xc for peak
  double peak=-INFINITY;
  for (t=0;t<=n_samp-n_seq;t++) {
    if (xc(t)>peak) {
      peak=xc(t);
    }
  }
  cout << "Maximum correlation found: " << db10(peak) << " dB" << endl;

  // Search for local maxima whose magnitude is near the peak value.
  // Expected number of samples in a frame.
  double expected_period=fs*.010*k_factor;
  cout << "Expected correlation period: " << expected_period << endl;
  int32 prev_peak=-1;
  for (t=1;t<=n_samp-n_seq-1;t++) {
    if ((xc(t)>peak*udb10(-4.0))&&(xc(t)>xc(t-1))&&(xc(t)>xc(t+1))) {
      if (prev_peak==-1) {
        prev_peak=t;
        continue;
      }
      // Have we missed detection of any previous sync signals?
      uint32 n_frames_skipped=MAX(0,round_i((t-prev_peak)/expected_period)-1);
      for (uint32 k=0;k<n_frames_skipped;k++) {
        cout << "Peak loc: " << round_i(prev_peak+(k+1)*expected_period) << " missing" << endl;
      }
      prev_peak=prev_peak+round_i(n_frames_skipped*expected_period);

      int32 n_dropped=round_i(expected_period-(t-prev_peak));
      cout << "Peak loc: " << t << " diff w/ prev: " << t-prev_peak << " n_dropped: " << n_dropped;
      if (abs(n_dropped)>100) {
        cout << " ***" << endl;
      } else if (abs(n_dropped)>10) {
        cout << " **" << endl;
      } else if (abs(n_dropped)>2) {
        cout << " *" << endl;
      } else {
        cout << endl;
      }
      prev_peak=t;
    }
  }

  return 0;
}

