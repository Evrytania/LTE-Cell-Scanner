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

