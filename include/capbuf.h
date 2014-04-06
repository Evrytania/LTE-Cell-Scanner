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


#ifdef HAVE_HACKRF
#include "hackrf.h"
#else
typedef struct hackrf_device hackrf_device;
#endif

#ifndef HAVE_CAPBUF_H
#define HAVE_CAPBUF_H

// Number of complex samples to capture.
#define CAPLENGTH 153600


typedef struct {
  std::vector <unsigned char> * buf;
  rtlsdr_dev_t * dev;
} callback_package_t;

#ifdef HAVE_HACKRF

typedef struct {
  std::vector <unsigned char> * buf;
  hackrf_device * dev;
} callback_hackrf_package_t;

#endif

double calculate_fc_programmed_in_context(
  // Inputs
  const double & fc_requested,
  const bool & use_recorded_data,
  const char * load_bin_filename,
  rtlsdr_dev_t * & dev
);

int read_header_from_bin(
  // input
  const char *bin_filename,
  // output, NAN represents invalid header info
  double & fc_requested,
  double & fc_programmed,
  double & fs_requested,
  double & fs_programmed
);

// Returns a capture buffer either from a file or from live data read
// from the dongle.
int capture_data(
  // Inputs
  const double & fc_requested,
  const double & correction,
  const bool & save_cap,
  const char * record_bin_filename,
  const bool & use_recorded_data,
  const char * load_bin_filename,
  const std::string & str,
  rtlsdr_dev_t * & dev,
  hackrf_device * & hackrf_dev,
  const dev_type_t::dev_type_t & dev_use,
  // Output
  itpp::cvec & capbuf,
  double & fc_programmed,
  double & fs_programmed,
  const bool & read_all_in_bin
);

#endif

