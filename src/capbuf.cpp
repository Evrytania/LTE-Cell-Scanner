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

#include <itpp/itbase.h>
#include <iomanip>
#include <sstream>
#include <queue>
#include "rtl-sdr.h"
#include "common.h"
#include "capbuf.h"
#include "macros.h"

using namespace itpp;
using namespace std;

// Number of complex samples to capture.
#define CAPLENGTH 153600

// This function produces a vector of captured data. The data can either
// come from live data received by the RTLSDR, or from a file containing
// previously captured data.
// Also, optionally, this function can save each set of captured data
// to a file.
void capture_data(
  // Inputs
  const double & fc,
  const double & correction,
  const bool & save_cap,
  const bool & use_recorded_data,
  const string & data_dir,
  // Output
  cvec & capbuf
) {
  // Filename used for recording or loading captured data.
  static uint32 capture_number=0;
  stringstream filename;
  filename << data_dir << "/capbuf_" << setw(4) << setfill('0') << capture_number << ".it";

  if (use_recorded_data) {
    // Read data from a file. Do not use live data.
    if (verbosity>=2) {
      cout << "Reading captured data from file: " << filename.str() << endl;
    }

    it_ifile itf(filename.str());
    itf.seek("capbuf");
    itf>>capbuf;
    itf.seek("fc");
    ivec fc_v;
    itf>>fc_v;
    if (fc!=fc_v(0)) {
      cout << "Warning: while reading capture buffer " << capture_number << ", the read" << endl;
      cout << "center frequency did not match the expected center frequency." << endl;
    }
    itf.close();

  } else {
    if (verbosity>=2) {
      cout << "Capturing live data" << endl;
    }
#define BLOCK_SIZE 16*16384
    uint8 * buffer=(uint8 *)malloc(BLOCK_SIZE*sizeof(uint8));

    // Center frequency
    if (rtlsdr_set_center_freq(dev,itpp::round(fc*correction))<0) {
      cerr << "Error: unable to set center frequency" << endl;
      exit(-1);
    }

    // Reset the buffer
    if (rtlsdr_reset_buffer(dev)<0) {
      cerr << "Error: unable to reset RTLSDR buffer" << endl;
      exit(-1);
    }

    // Read and store the data
    uint32 n_read=0;
    int n_read_current=0;
    uint32 n_saved=0;
    //capbuf.set_size(CAPLENGTH);
#define FIR_LENGTH 10
    cvec capbuf_0p5x(CAPLENGTH/2+FIR_LENGTH+10);
#ifndef NDEBUG
    //capbuf=NAN;
    capbuf_0p5x=NAN;
#endif
    while (true) {
      // Read some data
      if (rtlsdr_read_sync(dev,buffer,BLOCK_SIZE,&n_read_current)<0) {
        cerr << "Error: synchronous read failed" << endl;
        exit(-1);
      }
      if (n_read_current<BLOCK_SIZE) {
        cerr << "Error: short read; samples lost" << endl;
        exit(-1);
      }

      for (uint32 t=0;t<BLOCK_SIZE;t+=2) {
        n_read+=2;
        // Ignore first 20ms... Hopefully PLL will lock by then...
        if (n_read<19200) {
          continue;
        }
        capbuf_0p5x(n_saved++)=complex<double>((buffer[t]-127.0)/128.0,(buffer[t+1]-127.0)/128.0);
        if ((signed)n_saved==length(capbuf_0p5x)) {
          goto cbuf_full;
        }
      }
    }

    cbuf_full:
    free(buffer);
    if ((signed)n_saved!=length(capbuf_0p5x)) {
      cerr << "Error: unable to fill capture buffer..." << endl;
      exit(-1);
    }

    // Interpolate.
    // FIXME: Do proper interpolation.
    capbuf.set_size(CAPLENGTH);
    capbuf=NAN;
    for (uint32 t=0;t<CAPLENGTH;t+=2) {
      capbuf(t)=capbuf_0p5x(t>>1);
      capbuf(t+1)=(capbuf_0p5x(t>>1)+capbuf_0p5x((t>>1)+1))/2;
    }
  }

  // Save the capture data, if requested.
  if (save_cap) {
    if (verbosity>=2) {
      cout << "Saving captured data to file: " << filename.str() << endl;
    }
    it_file itf(filename.str(),true);
    itf << Name("capbuf") << capbuf;
    ivec fc_v(1);
    fc_v(0)=fc;
    itf << Name("fc") << fc_v;
    itf.close();
  }

  capture_number++;
}

