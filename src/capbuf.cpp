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

static void capbuf_rtlsdr_callback(
  unsigned char * buf,
  uint32_t len,
  void * ctx
) {
  vector <char> & capbuf_raw = *((vector <char> *)ctx);
  //cout << capbuf_raw.size() << endl;

  if (len==0) {
    cerr << "Error: received no samples from USB device..." << endl;
    exit(-1);
  }

  for (uint32 t=0;t<len;t++) {
    //cout << capbuf_raw.size() << endl;
    if (capbuf_raw.size()<CAPLENGTH*2) {
      capbuf_raw.push_back(buf[t]);
    }
    if (capbuf_raw.size()==CAPLENGTH*2) {
      rtlsdr_cancel_async(dev);
      break;
    }
  }
  //cout << capbuf_raw.size() << endl;
}

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

    // Read and store the data.
    // This will block until the call to rtlsdr_cancel_async().
    vector <char> capbuf_raw;
    capbuf_raw.reserve(CAPLENGTH*2);
    rtlsdr_read_async(dev,capbuf_rtlsdr_callback,(void *)&capbuf_raw,0,0);
    if (capbuf_raw.size()!=CAPLENGTH*2) {
      cerr << "Error: unable to read sufficient data from USB device" << endl;
      exit(-1);
    }

    // Convert to complex
    capbuf.set_size(CAPLENGTH);
#ifndef NDEBUG
    capbuf=NAN;
#endif
    for (uint32 t=0;t<CAPLENGTH;t++) {
      capbuf(t)=complex<double>((capbuf_raw[(t<<1)]-127.0)/128.0,(capbuf_raw[(t<<1)+1]-127.0)/128.0);
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

