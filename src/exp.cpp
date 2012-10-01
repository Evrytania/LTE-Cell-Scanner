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

// Main routine.
int main(
  const int argc,
  char * const argv[]
) {
  mat F(62,62);
  F=(double)0;
  for (uint8 t=0;t<62;t++) {
    uint8 lt=MAX(0,t-6);
    uint8 rt=MIN(61,t+6);
    uint8 n_av=rt-lt+1;
    cout << t;
    for (uint8 k=lt;k<=rt;k++) {
      F(t,k)=1.0/n_av;
      cout << k;
    }
    cout << endl;
  }
  mat FmI=F-eye(62);
  cout << F(1,1) << endl;
  cout << FmI(1,1) << endl;
  mat FF=FmI*FmI;
  double factor=0;
  for (uint8 t=0;t<62;t++) {
    factor+=FF(t,t);
  }
  factor=factor/62;
  cout << factor << endl;
  cout << 12.0/13.0 << endl;
  cout << (factor/(12.0/13.0)) << endl;
  cout << db10(factor/(12.0/13.0)) << endl;
  return 0;

  cvec seq_e=blnoise(100);
  cvec seq_e_int=interpft(seq_e,1000);
  cvec seq_e_int_dec=seq_e_int(itpp_ext::matlab_range(0,10,length(seq_e_int)-1));
  cout << db10(sigpower(seq_e-seq_e_int_dec)) << endl;

  cvec seq_o=blnoise(101);
  cvec seq_o_int=interpft(seq_o,1010);
  cvec seq_o_int_dec=seq_o_int(itpp_ext::matlab_range(0,10,length(seq_o_int)-1));
  cout << db10(sigpower(seq_e-seq_e_int_dec)) << endl;
  ABORT(-1);

  const uint16 n_id_cell=271;
  const uint8 n_id_1=floor(n_id_cell/3);
  const uint8 n_id_2=n_id_cell-3*n_id_1;

  cout << "Hello World!" << endl;

  cvec cap_data;
  itpp_ext::rtl_sdr_to_cvec(argv[1],cap_data);
  // Drop first 4 s for AGC to converge.
  cap_data=cap_data(FS_LTE/16*4,-1);
  uint32 n_samp=length(cap_data);
  n_samp=MIN(10e6,n_samp);
  cout << n_samp << endl;

  SSS_td sss_td;
  PSS_td pss_td;

  // Create the sequence we are looking for.
  const double f_off=39914.7;
  const double fc=739e6;
  const double k_factor=(fc-f_off)/fc;
  cvec seq=concat(sss_td(n_id_1,n_id_2,0),pss_td[n_id_2]);
  seq=fshift(seq,f_off,FS_LTE/16*k_factor);
  uint32 n_seq=length(seq);
  seq=conj(seq)/((double)n_seq);
  cout << n_seq << endl;

  vec xc(n_samp-n_seq+1);
  xc=NAN;
  uint32 t;
#ifdef _OPENMP
#pragma omp parallel for shared(seq,cap_data,xc) private(t)
#endif
  for (t=0;t<=n_samp-n_seq;t++) {
    if (mod(t,1000000)==0) {
      cout << t << endl;
    }
    xc(t)=sqr(elem_mult_sum(cap_data(t,t+n_seq-1),seq));
  }

  //cout << xc << endl;

  // Search xc
  double peak=-INFINITY;
  for (t=0;t<=n_samp-n_seq;t++) {
    if (xc(t)>peak) {
      peak=xc(t);
    }
  }
  cout << db10(peak) << " dB" << endl;

  // Search for local maxima near the peak value
  uint32 prev_peak=0;
  for (t=1;t<=n_samp-n_seq-1;t++) {
    if ((xc(t)>peak*udb10(-3.0))&&(xc(t)>xc(t-1))&&(xc(t)>xc(t+1))) {
      double factor=(t-prev_peak)/(19200.0*k_factor);
      cout << t << " " << t-prev_peak << " " << factor;
      if (abs(round_i(factor)-factor)>.0001) {
        cout << " ***" << endl;
      } else {
        cout << endl;
      }
      prev_peak=t;
    }
  }

  return 0;
}

