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
#include <itpp/comm/convcode.h>
#include <itpp/comm/crc.h>
#include <itpp/signal/transforms.h>
#include <list>
#include <complex>
#include <boost/math/special_functions/gamma.hpp>
#include "rtl-sdr.h"
#include "common.h"
#include "lte_lib.h"
#include "constants.h"
#include "macros.h"
#include "searcher.h"
#include "itpp_ext.h"
#include "dsp.h"

using namespace itpp;
using namespace std;

// Implementation of the LTE PN generator.
// A faster implementation is possible using uint32 to store the state.
// This function was coded using vectors and matrices so as to match
// the Matlab implementation as directly as possible.
bvec lte_pn(
  const uint32 & c_init,
  const uint32 & len
) {
  // Initialize the state machines
  bvec x1(31);
  bvec x2(31);
  uint32 c_init_temp(c_init);
  for (uint8 t=0;t<31;t++) {
    x1(t)=0;
    x2(t)=(c_init_temp&1);
    c_init_temp=c_init_temp>>1;
  }
  x1(0)=1;

  // These matrices advance the state machines forward by 1600 steps.
  const static bin x1_p1600_cmat[] = {
    0,1,0,1,1,1,1,0,0,1,0,0,1,0,0,0,0,1,0,1,1,0,0,0,0,1,0,0,0,0,0,
    0,0,1,0,1,1,1,1,0,0,1,0,0,1,0,0,0,0,1,0,1,1,0,0,0,0,1,0,0,0,0,
    0,0,0,1,0,1,1,1,1,0,0,1,0,0,1,0,0,0,0,1,0,1,1,0,0,0,0,1,0,0,0,
    0,0,0,0,1,0,1,1,1,1,0,0,1,0,0,1,0,0,0,0,1,0,1,1,0,0,0,0,1,0,0,
    0,0,0,0,0,1,0,1,1,1,1,0,0,1,0,0,1,0,0,0,0,1,0,1,1,0,0,0,0,1,0,
    0,0,0,0,0,0,1,0,1,1,1,1,0,0,1,0,0,1,0,0,0,0,1,0,1,1,0,0,0,0,1,
    1,0,0,1,0,0,0,1,0,1,1,1,1,0,0,1,0,0,1,0,0,0,0,1,0,1,1,0,0,0,0,
    0,1,0,0,1,0,0,0,1,0,1,1,1,1,0,0,1,0,0,1,0,0,0,0,1,0,1,1,0,0,0,
    0,0,1,0,0,1,0,0,0,1,0,1,1,1,1,0,0,1,0,0,1,0,0,0,0,1,0,1,1,0,0,
    0,0,0,1,0,0,1,0,0,0,1,0,1,1,1,1,0,0,1,0,0,1,0,0,0,0,1,0,1,1,0,
    0,0,0,0,1,0,0,1,0,0,0,1,0,1,1,1,1,0,0,1,0,0,1,0,0,0,0,1,0,1,1,
    1,0,0,1,0,1,0,0,1,0,0,0,1,0,1,1,1,1,0,0,1,0,0,1,0,0,0,0,1,0,1,
    1,1,0,1,1,0,1,0,0,1,0,0,0,1,0,1,1,1,1,0,0,1,0,0,1,0,0,0,0,1,0,
    0,1,1,0,1,1,0,1,0,0,1,0,0,0,1,0,1,1,1,1,0,0,1,0,0,1,0,0,0,0,1,
    1,0,1,0,0,1,1,0,1,0,0,1,0,0,0,1,0,1,1,1,1,0,0,1,0,0,1,0,0,0,0,
    0,1,0,1,0,0,1,1,0,1,0,0,1,0,0,0,1,0,1,1,1,1,0,0,1,0,0,1,0,0,0,
    0,0,1,0,1,0,0,1,1,0,1,0,0,1,0,0,0,1,0,1,1,1,1,0,0,1,0,0,1,0,0,
    0,0,0,1,0,1,0,0,1,1,0,1,0,0,1,0,0,0,1,0,1,1,1,1,0,0,1,0,0,1,0,
    0,0,0,0,1,0,1,0,0,1,1,0,1,0,0,1,0,0,0,1,0,1,1,1,1,0,0,1,0,0,1,
    1,0,0,1,0,1,0,1,0,0,1,1,0,1,0,0,1,0,0,0,1,0,1,1,1,1,0,0,1,0,0,
    0,1,0,0,1,0,1,0,1,0,0,1,1,0,1,0,0,1,0,0,0,1,0,1,1,1,1,0,0,1,0,
    0,0,1,0,0,1,0,1,0,1,0,0,1,1,0,1,0,0,1,0,0,0,1,0,1,1,1,1,0,0,1,
    1,0,0,0,0,0,1,0,1,0,1,0,0,1,1,0,1,0,0,1,0,0,0,1,0,1,1,1,1,0,0,
    0,1,0,0,0,0,0,1,0,1,0,1,0,0,1,1,0,1,0,0,1,0,0,0,1,0,1,1,1,1,0,
    0,0,1,0,0,0,0,0,1,0,1,0,1,0,0,1,1,0,1,0,0,1,0,0,0,1,0,1,1,1,1,
    1,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,1,1,0,1,0,0,1,0,0,0,1,0,1,1,1,
    1,1,0,1,0,0,0,0,0,0,1,0,1,0,1,0,0,1,1,0,1,0,0,1,0,0,0,1,0,1,1,
    1,1,1,1,1,0,0,0,0,0,0,1,0,1,0,1,0,0,1,1,0,1,0,0,1,0,0,0,1,0,1,
    1,1,1,0,1,1,0,0,0,0,0,0,1,0,1,0,1,0,0,1,1,0,1,0,0,1,0,0,0,1,0,
    0,1,1,1,0,1,1,0,0,0,0,0,0,1,0,1,0,1,0,0,1,1,0,1,0,0,1,0,0,0,1,
    1,0,1,0,1,0,1,1,0,0,0,0,0,0,1,0,1,0,1,0,0,1,1,0,1,0,0,1,0,0,0
  };
  const static bmat x1_p1600_mat(x1_p1600_cmat,31,31);
  const static bin x2_p1600_cmat [] = {
    0,1,1,1,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,1,1,0,0,1,0,0,0,0,0,0,0,
    0,0,1,1,1,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,1,1,0,0,1,0,0,0,0,0,0,
    0,0,0,1,1,1,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,1,1,0,0,1,0,0,0,0,0,
    0,0,0,0,1,1,1,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,1,1,0,0,1,0,0,0,0,
    0,0,0,0,0,1,1,1,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,1,1,0,0,1,0,0,0,
    0,0,0,0,0,0,1,1,1,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,1,1,0,0,1,0,0,
    0,0,0,0,0,0,0,1,1,1,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,1,1,0,0,1,0,
    0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,1,1,0,0,1,
    1,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,1,1,0,0,
    0,1,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,1,1,0,
    0,0,1,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,1,1,
    1,1,1,0,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,1,
    1,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,
    0,1,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,1,0,0,0,1,0,0,0,1,0,
    0,0,1,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,1,0,0,0,1,0,0,0,1,
    1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,1,0,0,0,1,0,0,0,
    0,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,1,0,0,0,1,0,0,
    0,0,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,1,0,0,0,1,0,
    0,0,0,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,1,0,0,0,1,
    1,1,1,1,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,1,0,0,0,
    0,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,1,0,0,
    0,0,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,1,0,
    0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,1,
    1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,
    0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1,0,0,0,
    0,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1,0,0,
    0,0,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1,0,
    0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1,
    1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,
    1,0,0,0,1,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,
    1,0,1,1,0,1,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0
  };
  const static bmat x2_p1600_mat(x2_p1600_cmat,31,31);

  // Advance the state machines by 1600 clocks.
  x1=x1_p1600_mat*x1;
  x2=x2_p1600_mat*x2;

  // Create output vector
  bvec rv(len);
  for (uint32 t=0;t<len;t++) {
    // Store output
    rv(t)=x1(0)+x2(0);
    // Advance state machines
    bool x1_new=x1(0)+x1(3);
    bool x2_new=x2(0)+x2(1)+x2(2)+x2(3);
    for (uint8 k=0;k<30;k++) {
      x1(k)=x1(k+1);
      x2(k)=x2(k+1);
    }
    x1(30)=x1_new;
    x2(30)=x2_new;
  }

  return rv;
}

// Instantiate the static members
//vector <cvec> PSS_fd::table(3);
//vector <cvec> PSS_td::table(3);

// PSS
// Return the PSS in the frequency domain.
cvec pss_fd_calc(const uint8 & t) {
  const int zc_map_donotuse[3]={25,29,34};
  const vector <int> zc_map(zc_map_donotuse,zc_map_donotuse+3);
  cvec r=exp((complex<double>(0,-1)*pi*zc_map[t]/63)*elem_mult(ivec("0:62"),ivec("1:63")));
  r.del(31);
  return r;
}

// Initialize table
PSS_fd::PSS_fd(void) {
  table.resize(3);
  for (uint8 t=0;t<3;t++) {
    table[t]=pss_fd_calc(t);
  }
}

// Return the requested PSS in the frequency domain.
const cvec & PSS_fd::operator[](const uint8 & idx) const {
  return table[idx];
}

// Initialize PSS time domain table
PSS_td::PSS_td(void) {
  cvec fd;
  cvec td;
  cvec idft_in;
  table.resize(3);
  for (uint8 t=0;t<3;t++) {
    fd=pss_fd_calc(t);
    idft_in=concat(zeros_c(1),fd(31,61),zeros_c(65),fd(0,30));
    td=idft(idft_in)*sqrt(128.0/62.0);
    table[t]=concat(td(119,127),td);
  }
}

// Return the requested PSS in the time domain.
const cvec & PSS_td::operator[](const uint8 & idx) const {
  return table[idx];
}

// Instantiate the static members
//vector < vector < vector <ivec> > > SSS_fd::table;

// Calculate the SSS in the frequency domain.
ivec sss_fd_calc(
  const uint8 n_id_1,
  const uint8 n_id_2,
  const uint8 slot_num
) {
  // Calculate m0 and m1
  const uint32 qp=floor_i(n_id_1/30);
  const uint32 q=floor_i((n_id_1+qp*(qp+1)/2)/30);
  const uint32 mp=n_id_1+q*(q+1)/2;
  const uint32 m0=itpp_ext::matlab_mod(mp,31);
  const uint32 m1=itpp_ext::matlab_mod(m0+floor_i(mp/31)+1,31);

  //s_td=[0 0 0 0 1];
  //for t=1:26
  //  s_td=[s_td mod(s_td(end-2)+s_td(end-4),2)];
  //end
  ivec s_td("0 0 0 0 1 0 0 1 0 1 1 0 0 1 1 1 1 1 0 0 0 1 1 0 1 1 1 0 1 0 1");
  s_td=1-2*s_td;

  //c_td=[0 0 0 0 1];
  //for t=1:26
  //  c_td=[c_td mod(c_td(end-1)+c_td(end-4),2)];
  //end
  ivec c_td("0 0 0 0 1 0 1 0 1 1 1 0 1 1 0 0 0 1 1 1 1 1 0 0 1 1 0 1 0 0 1");
  c_td=1-2*c_td;

  //z_td=[0 0 0 0 1];
  //for t=1:26
  //  z_td=[z_td mod(z_td(end)+z_td(end-2)+z_td(end-3)+z_td(end-4),2)];
  //end
  ivec z_td("0 0 0 0 1 1 1 0 0 1 1 0 1 1 1 1 1 0 1 0 0 0 1 0 0 1 0 1 0 1 1");
  z_td=1-2*z_td;

  ivec s0_m0=s_td(itpp_ext::matlab_mod(itpp_ext::matlab_range(m0,30+m0),31));
  ivec s1_m1=s_td(itpp_ext::matlab_mod(itpp_ext::matlab_range(m1,30+m1),31));

  ivec c0=c_td(itpp_ext::matlab_mod(itpp_ext::matlab_range(n_id_2,30+n_id_2),31));
  ivec c1=c_td(itpp_ext::matlab_mod(itpp_ext::matlab_range(n_id_2+3,30+n_id_2+3),31));

  ivec z1_m0=z_td(itpp_ext::matlab_mod(itpp_ext::matlab_range(0,30)+itpp_ext::matlab_mod(m0,8),31));
  ivec z1_m1=z_td(itpp_ext::matlab_mod(itpp_ext::matlab_range(0,30)+itpp_ext::matlab_mod(m1,8),31));

  ivec ssc1;
  ivec ssc2;
  if (slot_num==0) {
    ssc2=elem_mult(s1_m1,c1,z1_m0);
    ssc1=elem_mult(s0_m0,c0);
  } else {
    ssc2=elem_mult(s0_m0,c1,z1_m1);
    ssc1=elem_mult(s1_m1,c0);
  }

  // Interleave the SSC
  imat sss(2,31);
  sss.set_row(0,ssc1);
  sss.set_row(1,ssc2);

  return cvectorize(sss);
}

// Initialize table
SSS_fd::SSS_fd(void) {
  table=vector < vector < vector < ivec > > > (168,vector< vector < ivec > >(3, vector < ivec > (2)));
  for (uint8 n_id_1=0;n_id_1<168;n_id_1++) {
    for (uint8 n_id_2=0;n_id_2<3;n_id_2++) {
      for (uint8 n_slot=0;n_slot<2;n_slot++) {
        table[n_id_1][n_id_2][n_slot]=sss_fd_calc(n_id_1,n_id_2,n_slot*10);
      }
    }
  }
}

// Return the requested SSS in the frequency domain.
const ivec & SSS_fd::operator()(const uint8 & n_id_1,const uint8 & n_id_2,const uint8 & n_slot) const {
  return table[n_id_1][n_id_2][n_slot!=0];
}

// Instantiate the static members
vector < vector < vector <cvec> > > SSS_td::table;

// Initialize table
SSS_td::SSS_td(void) {
  table=vector < vector < vector < cvec > > > (168,vector< vector < cvec > >(3, vector < cvec > (2)));
  ivec fd;
  cvec td;
  cvec idft_in;
  for (uint8 n_id_1=0;n_id_1<168;n_id_1++) {
    for (uint8 n_id_2=0;n_id_2<3;n_id_2++) {
      for (uint8 n_slot=0;n_slot<2;n_slot++) {
        fd=sss_fd_calc(n_id_1,n_id_2,n_slot*10);
        idft_in=concat(zeros_c(1),to_cvec(fd(31,61)),zeros_c(65),to_cvec(fd(0,30)));
        td=idft(idft_in)*sqrt(128.0/62.0);
        table[n_id_1][n_id_2][n_slot]=concat(td(119,127),td);
      }
    }
  }
}

// Return the requested SSS in the time domain.
const cvec & SSS_td::operator()(const uint8 & n_id_1,const uint8 & n_id_2,const uint8 & n_slot) const {
  return table[n_id_1][n_id_2][n_slot!=0];
}

// Function to calculate the reference symbols.
// This function executes rather slowly and hence it is usually better
// to use the RS_DL class which precomputes all the needed DL RS at once.
cvec rs_dl_calc(
  const uint32 & slot_num,
  const uint32 & sym_num,
  const uint32 & n_id_cell,
  const uint32 & n_rb_dl,
  const cp_type_t::cp_type_t & cp_type
) {
  // Create c
  const uint32 n_cp=(cp_type==cp_type_t::NORMAL);
  const uint32 c_init=(1<<10)*(7*(slot_num+1)+sym_num+1)*(2*n_id_cell+1)+2*n_id_cell+n_cp;
  const ivec c=to_ivec(lte_pn(c_init,4*N_RB_MAXDL));

  // Create symbols to be used if n_rb_maxdl RB's were used on the DL.
  cvec r_l_ns=(1/pow(2,0.5))*((1-2*c(itpp_ext::matlab_range(0,2,4*N_RB_MAXDL-1)))+J*(1-2*c(itpp_ext::matlab_range(1,2,4*N_RB_MAXDL-1))));

  // Select the appropriate subset of r_l_ns.
  cvec r=r_l_ns(N_RB_MAXDL-n_rb_dl,2*n_rb_dl+N_RB_MAXDL-n_rb_dl-1);

  return r;
}

// Calculate the shift a certain port will use for a certain ofdm symbol.
double rs_dl_shift_calc(
  const uint8 & slot_num,
  const uint8 & sym_num,
  const uint8 & port_num,
  const cp_type_t::cp_type_t & cp_type,
  const uint16 & n_id_cell
) {
  uint8 n_symb_dl=(cp_type==cp_type_t::NORMAL)?7:6;

  double v=NAN;
  if ((port_num==0)&&(sym_num==0))
    v=0;
  else if ((port_num==0)&&(sym_num==n_symb_dl-3))
    v=3;
  else if ((port_num==1)&&(sym_num==0))
    v=3;
  else if ((port_num==1)&&(sym_num==n_symb_dl-3))
    v=0;
  else if ((port_num==2)&&(sym_num==1))
    v=3*(slot_num&1);
  else if ((port_num==3)&&(sym_num==1))
    v=3+3*(slot_num&1);

  return mod(v+n_id_cell,6);
}

// Constructor
RS_DL::RS_DL(
  const uint16 & n_id_cell,
  const uint8 & n_rb_dl,
  const cp_type_t::cp_type_t & cp_type
) {
  // Derive some values
  n_symb_dl=7;
  if (cp_type==cp_type_t::EXTENDED)
    n_symb_dl=6;

  // Fill in the tables
  table.resize(20*n_symb_dl);
  shift_table=mat(20*n_symb_dl,4);
  // Table must be initialized with NAN since NAN is used
  // to signal that there are no RS in this OFDM symbol.
  shift_table=NAN;
  for (uint8 slot_num=0;slot_num<20;slot_num++) {
    for (uint8 t=0;t<3;t++) {
      uint8 sym_num=(t==2)?(n_symb_dl-3):t;
      table[slot_num*n_symb_dl+sym_num]=rs_dl_calc(slot_num,sym_num,n_id_cell,n_rb_dl,cp_type);
      if ((t==0)||(t==2)) {
        shift_table(slot_num*n_symb_dl+sym_num,0)=rs_dl_shift_calc(slot_num,sym_num,0,cp_type,n_id_cell);
        shift_table(slot_num*n_symb_dl+sym_num,1)=rs_dl_shift_calc(slot_num,sym_num,1,cp_type,n_id_cell);
      } else {
        shift_table(slot_num*n_symb_dl+sym_num,2)=rs_dl_shift_calc(slot_num,sym_num,2,cp_type,n_id_cell);
        shift_table(slot_num*n_symb_dl+sym_num,3)=rs_dl_shift_calc(slot_num,sym_num,3,cp_type,n_id_cell);
      }
    }
  }
}

// Return the requested precomputed RS.
const cvec & RS_DL::get_rs (
  const uint8 & slot_num,
  const uint8 & sym_num
) const {
  ASSERT(slot_num<20);
  ASSERT(sym_num<n_symb_dl);
  return table[slot_num*n_symb_dl+sym_num];
}

// Return the requested shift value.
double RS_DL::get_shift (
  const uint8 & slot_num,
  const uint8 & sym_num,
  const uint8 & port_num
) const {
  ASSERT(slot_num<20);
  ASSERT(sym_num<n_symb_dl);
  ASSERT(port_num<4);
  return shift_table(slot_num*n_symb_dl+sym_num,port_num);
}

// Perform LTE ratematching for convolutionally encoded bits d to fit
// an e vector of length n_e.
cvec lte_conv_ratematch(
  const cmat & d,
  const uint32 & n_e
) {
  const uint8 n_c=32;
  const uint32 n_r=ceil((double)d.cols()/n_c);
  ASSERT(d.rows()==3);

  // Subblock interleaving
  const static int perm_pattern_c[]={1,17,9,25,5,21,13,29,3,19,11,27,7,23,15,31,0,16,8,24,4,20,12,28,2,18,10,26,6,22,14,30};
  const ivec perm_pattern(perm_pattern_c,32);
  cmat v(3,n_r*n_c);
#ifndef NDEBUG
  v=NAN;
#endif
  for (uint8 t=0;t<3;t++) {
    cvec temp_row=d.get_row(t);
    cvec temp_nan(n_r*n_c-d.cols());
    temp_nan=NAN;
    temp_row=concat(temp_nan,temp_row);
    cmat y=transpose(reshape(temp_row,n_c,n_r));

    // Permute the columns
    cmat y_perm(n_r,n_c);
#ifndef NDEBUG
    y_perm=NAN;
#endif
    for (uint8 k=0;k<32;k++) {
      y_perm.set_col(k,y.get_col(perm_pattern(k)));
    }

    // Assign
    v.set_row(t,cvectorize(y_perm));
  }

  // Bit collection
  cvec w=cvectorize(transpose(v));

  // Selection
  cvec e(n_e);
#ifndef NDEBUG
  e=NAN;
#endif
  uint32 k=0;
  uint32 j=0;
  while (k<n_e) {
    if (isfinite(w(j).real())) {
      e(k)=w(j);
      k++;
    }
    j=mod(j+1,3*n_r*n_c);
  }

  return e;
}

// Note that this function assumes that e_est is ln(P(e_est==0|r)/P(e_est==1|r))
// This is different from the matlab function which assumes e_est
// is simply P(e_est==0|r).
// This function also returns ln(P(d_est==0|r)/P(d_est==1|r)).
mat lte_conv_deratematch(
  const vec & e_est,
  const uint32 & n_c
) {
  // Probe lte_conv_ratematch to find out which bit came from where.
  // There are better ways to do this...
  mat d_probe_real=repmat(itpp_ext::matlab_range(0.0,2.0),1,n_c,false);
  mat d_probe_imag=repmat(itpp_ext::matlab_range(0.0,n_c-1.0),3,1,true);
  cmat d_probe=to_cmat(d_probe_real,d_probe_imag);
  cvec e_probe=lte_conv_ratematch(d_probe,length(e_est));

  // Convert from probabilities to values along the x axis, assuming symbols
  // are +1/-1 and the received noise power is 1.
  // This is no longer needed since the input to this function is
  // ln(P(e_est==0|r)/P(e_est==1|r)).
  // If a BPSK symbol with value r is received in noise with a power of 2,
  // ln(P(e_est==0|r)/P(e_est==1|r)) is simply r. Thus, if multiple observations
  // of the same transmitted symbol are available and we know
  // ln(P(e_est==0|r)/P(e_est==1|r)) for each observation, we can simply average
  // all of the ln(P(e_est==0|r)/P(e_est==1|r)) values to obtained an combined
  // estimate.
  const vec e_x=e_est;

  // Combine all the observations of the same coded bit.
  mat d_x(3,n_c);
  d_x=0.0;
  imat d_x_count(3,n_c);
  d_x_count=0;
  for (int32 t=0;t<length(e_est);t++) {
    uint8 r=real(e_probe(t));
    uint32 c=imag(e_probe(t));
    d_x(r,c)+=e_x(t);
    d_x_count(r,c)++;
  }
  for (uint8 r=0;r<3;r++) {
    for (uint32 c=0;c<n_c;c++) {
      if (d_x_count(r,c)>1) {
        d_x(r,c)=d_x(r,c)/d_x_count(r,c);
      }
    }
  }

  // No longer needed
  // Convert back to probabilities.
  //mat p1=exp(-sqr(d_x-1)/2);
  //mat m1=exp(-sqr(d_x+1)/2);
  //mat d_est=elem_div(p1,p1+m1);

  return d_x;
}

bmat lte_conv_encode(
  const bvec & c
) {
  Convolutional_Code coder;
  ivec generator(3);
  generator(0)=0133;
  generator(1)=0171;
  generator(2)=0165;
  coder.set_generator_polynomials(generator,7);

  bvec d_vec=coder.encode_tailbite(c);
  bmat d=reshape(d_vec,3,length(c));
  return d;
}

// Note that this function assumes that d_est is ln(P(d_est==0|r)/P(d_est==1|r))
// This is different from the matlab function which assumes d_est
// is simply P(d_est==0|r).
bvec lte_conv_decode(
  const mat & d_est
) {
  Convolutional_Code coder;
  ivec generator(3);
  generator(0)=0133;
  generator(1)=0171;
  generator(2)=0165;
  coder.set_generator_polynomials(generator,7);

  vec d_est_v=cvectorize(d_est);
  bvec c_est=coder.decode_tailbite(d_est_v);
  return c_est;
}

// Instantiate static member
//Array <cvec> Mod_map::table(3);

// Constructor
// This was implemented as a class so that the (slow) initialization only
// needs to be done once, when the program starts.
Mod_map::Mod_map(void) {
  ivec map_qam_real("1 1 -1 -1");
  ivec map_qam_imag("1 -1 1 -1");
  ivec map_qam16_real("1 1 3 3 1 1 3 3 -1 -1 -3 -3 -1 -1 -3 -3");
  ivec map_qam16_imag("1 3 1 3 -1 -3 -1 -3 1 3 1 3 -1 -3 -1 -3");
  ivec map_qam64_real("3 3 1 1 3 3 1 1 5 5 7 7 5 5 7 7 3 3 1 1 3 3 1 1 5 5 7 7 5 5 7 7 -3 -3 -1 -1 -3 -3 -1 -1 -5 -5 -7 -7 -5 -5 -7 -7 -3 -3 -1 -1 -3 -3 -1 -1 -5 -5 -7 -7 -5 -5 -7 -7");
  ivec map_qam64_imag("3 1 3 1 5 7 5 7 3 1 3 1 5 7 5 7 -3 -1 -3 -1 -5 -7 -5 -7 -3 -1 -3 -1 -5 -7 -5 -7 3 1 3 1 5 7 5 7 3 1 3 1 5 7 5 7 -3 -1 -3 -1 -5 -7 -5 -7 -3 -1 -3 -1 -5 -7 -5 -7");
  table.set_size(3);
  table(0)=to_cvec(map_qam_real,map_qam_imag)/sqrt(2);
  table(1)=to_cvec(map_qam16_real,map_qam16_imag)/sqrt(10);
  table(2)=to_cvec(map_qam64_real,map_qam64_imag)/sqrt(42);
}

const itpp::cvec & Mod_map::operator()(
  const modulation_t::modulation_t & mod
) const {
  if (mod==modulation_t::QAM) {
    return table(0);
  } else if (mod==modulation_t::QAM16) {
    return table(1);
  } else if (mod==modulation_t::QAM64) {
    return table(2);
  } else {
    throw("Check code!");
  }
};

// Modulate a vector of bits to QAM/ QAM16/ QAM64 according to LTE specs.
cvec lte_modulate(
  const bvec & bits,
  const modulation_t::modulation_t modulation
) {
  uint8 bps;
  if (modulation==modulation_t::QAM) {
    bps=2;
  } else if (modulation==modulation_t::QAM16) {
    bps=4;
  } else if (modulation==modulation_t::QAM64) {
    bps=6;
  } else {
    throw("Check code!!!");
  }
  ASSERT(mod(length(bits),bps)==0);

  Modulator <complex <double> > modulator(ROM_TABLES.mod_map(modulation),itpp_ext::matlab_range(0,(1<<bps)-1));
  cvec retval=modulator.modulate_bits(bits);

  return retval;
}

// Assumes that the channel has already been removed from the received
// signal and the np is the amount of noise present in each sample.
// This function returns ln(P(b==0|syms)/P(b==1|syms)).
vec lte_demodulate(
  const cvec & syms,
  const vec & np,
  const modulation_t::modulation_t & modulation
) {
  uint8 bps;
  if (modulation==modulation_t::QAM) {
    bps=2;
  } else if (modulation==modulation_t::QAM16) {
    bps=4;
  } else if (modulation==modulation_t::QAM64) {
    bps=6;
  } else {
    throw("Check code!!!");
  }

  Modulator <complex <double> > modulator(ROM_TABLES.mod_map(modulation),itpp_ext::matlab_range(0,(1<<bps)-1));

  cvec gain=1.0/to_cvec(sqrt(np));
  vec retval=modulator.demodulate_soft_bits(elem_mult(syms,gain),gain,1);

  return retval;
}

// Calculate one of the LTE CRC's.
bvec lte_calc_crc(
  const bvec & a,
  const crc_t crc
) {
  bvec poly;
  if (crc==CRC8) {
    poly=bvec("1 1 0 0 1 1 0 1 1");
  } else if (crc==CRC16) {
    poly=bvec("1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1");
  } else if (crc==CRC24A) {
    poly=bvec("1 1 0 0 0 0 1 1 0 0 1 0 0 1 1 0 0 1 1 1 1 1 0 1 1");
  } else if (crc==CRC24B) {
    poly=bvec("1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1 1");
  } else {
    throw("Check code...");
  }

  // Set up the generator
  CRC_Code crc_calc;
  crc_calc.set_generator(poly);

  // Calculate parity
  bvec p;
  crc_calc.parity(a,p);

  return p;
}

