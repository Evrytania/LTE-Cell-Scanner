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
#include <itpp/signal/transforms.h>
#include <complex>
#include <cmath>
#include <boost/math/special_functions/gamma.hpp>
#include "rtl-sdr.h"
#include "common.h"
#include "macros.h"
#include "itpp_ext.h"
#include "dsp.h"

using namespace itpp;
using namespace std;

/*
// Shift vector seq up by f Hz assuming that seq was sampled at fs Hz.
cvec fshift(const cvec &seq,const double f,const double fs) {
  complex <double> k=complex<double>(0,pi*f/(fs/2));
  uint32 len=length(seq);
  cvec r(len);
  for (uint32 t=0;t<len;t++) {
    r(t)=seq(t)*exp(k*((double)t));
  }
  return r;
}
// Shift vector seq up by f Hz assuming that seq was sampled at 2 Hz.
cvec fshift(const cvec &seq,const double f) {
  return fshift(seq,f,2);
}
*/

// Perform FFT based interpolation. Assuming input signal is cyclically
// repeating signal sampled at M points, return a cyclically repeating
// signal that is sampled at N points.
cvec interpft(
  const cvec & x,
  const uint32 & n_y_pre
) {
  const uint32 n_x=length(x);

  // If decimation is requested, first interpolate by an integer factor and
  // then decimate.
  uint32 dec_factor;
  uint32 n_y;
  if (n_y_pre<n_x) {
    dec_factor=floor(n_x/n_y_pre)+1;
    n_y=n_y_pre*dec_factor;
  } else {
    dec_factor=1;
    n_y=n_y_pre;
  }

  // Interpolate
  cvec x_fd=fft(x);
  uint32 nyqst=floor(n_x/2);
  cvec y_fd_interp=concat(x_fd(0,nyqst),zeros_c(n_y-n_x),x_fd(nyqst+1,n_x-1));

  // Treat the nyquist sample specially when x is of even length.
  if (mod(n_x,2)==0) {
    y_fd_interp(nyqst)=y_fd_interp(nyqst)/2;
    y_fd_interp(nyqst+n_y-n_x)=y_fd_interp(nyqst);
  }
  cvec y_interp=ifft(y_fd_interp);

  // Decimate and scale
  cvec y(n_y);
  uint32 idx=0;
  for (uint32 t=0;t<n_y;t++) {
    y(t)=(((double)n_y)/n_x)*y_interp(idx);
    idx+=dec_factor;
  }

  return y;
}

