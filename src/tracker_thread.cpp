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
#include <curses.h>
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
  double tp;
  double sp;
  double sp_raw;
  double np;
  cvec ce_filt;
} ce_filt_fifo_pdu_t;
typedef struct {
  uint8 slot_num;
  uint8 sym_num;
  double tp;
  double sp;
  double sp_raw;
  double np;
  cvec ce_interp;
} ce_interp_fifo_pdu_t;
typedef struct {
  cvec syms;
  cmat ce;
  vec sp;
  vec np;
} mib_fifo_pdu_t;

// Pop 128 time domain samples from the fifo, convert to the frequency
// domain, and extract the 72 desired subcarriers.
void get_fd(
  tracked_cell_t & tracked_cell,
  const double & fc_requested,
  const double & fc_programmed,
  const double & fs_programmed,
  const uint8 & slot_num,
  const uint8 & sym_num,
  //const ivec & cn,
  double & bulk_phase_offset,
  cvec & syms,
  double & frequency_offset,
  double & frame_timing
) {
  td_fifo_pdu_t pdu;
  {
    // Lock the fifo until we obtain one element.
    boost::mutex::scoped_lock lock(tracked_cell.fifo_mutex);
    if (tracked_cell.fifo.empty()) {
      tracked_cell.fifo_condition.wait(lock);
    }
    ASSERT(!tracked_cell.fifo.empty());
#ifndef NDEBUG
    if ((tracked_cell.fifo.front().slot_num!=slot_num)||(tracked_cell.fifo.front().sym_num!=sym_num)) {
      // We should never get here...
      cerr << "Error: cell tracker synchronization error! Check code!" << endl;
      ABORT(-1);
    }
#endif
    // Copy locally
    pdu=tracked_cell.fifo.front();
    // POP the fifo
    tracked_cell.fifo.pop();
  }

  // Convert to frequency domain and extract 6 center RB's.
  // Also perform FOC to remove ICI
  frequency_offset=pdu.frequency_offset;
  frame_timing=pdu.frame_timing;
  double k_factor=(fc_requested-frequency_offset)/fc_programmed;

  // Directly manipulate data on the fifo to minimize vector copies.
  // Remove ICI
  fshift_inplace(pdu.data,-frequency_offset,fs_programmed*k_factor);
  // Remove the 2 sample delay
  cvec dft_in(128);
  //dft_in=concat(dft_in(2,-1),dft_in(0,1));
  for (uint8 t=0;t<126;t++) {
    dft_in(t)=pdu.data.get(t+2);
  }
  dft_in(126)=pdu.data.get(0);
  dft_in(127)=pdu.data.get(1);
  cvec dft_out=dft(dft_in);
  //syms=concat(dft_out.right(36),dft_out.mid(1,36));
  syms.set_size(72);
  for (uint8 t=0;t<36;t++) {
    syms(t+36)=dft_out(t+1);
    syms(t)=dft_out(92+t);
  }
  // Compensate for the fact that the DFT was located improperly and also
  // for the bulk phase offset due to frequency errors.
  uint8 n_samp_elapsed;
  // How many time samples have passed since the previous DFT?
  if (tracked_cell.cp_type==cp_type_t::EXTENDED) {
    n_samp_elapsed=128+32;
  } else {
    n_samp_elapsed=(sym_num==0)?128+10:128+9;
  }
  //syms=exp(J*bulk_phase_offset)*elem_mult(syms,exp((-J*2*pi*tracked_cell.fifo.front().late/128)*cn));
  complex <double> coeff;
  double phase;
  const double k=2*pi*pdu.late/128;
  bulk_phase_offset=WRAP(bulk_phase_offset+2*pi*n_samp_elapsed*(1/(FS_LTE/16))*-frequency_offset,-pi,pi);
  const complex <double> bpo_coeff=complex<double>(cos(bulk_phase_offset),sin(bulk_phase_offset));
  for (uint8 t=1;t<=36;t++) {
    phase=-k*t;
    coeff.real()=cos(phase);
    coeff.imag()=sin(phase);
    syms(35+t)*=bpo_coeff*coeff;
    coeff.imag()=-coeff.imag();
    syms(36-t)*=bpo_coeff*coeff;
  }
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
  //if (rs_prev.frame_timing!=rs_next.frame_timing)
  //  return;

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
  double k_factor=(global_thread_data.fc_requested-frequency_offset)/global_thread_data.fc_programmed;
  double residual_f=arg(foe_comb)/(2*pi)/(0.0005+WRAP(rs_next.frame_timing-rs_prev.frame_timing,-19200.0/2,19200.0/2)*(1/(global_thread_data.fs_programmed*k_factor)));
  double residual_f_np=MAX(foe_comb_np/2,.001);
  //residual_f+=1000;
  //cout << residual_f << " f " << db10(residual_f_np) << endl;
  //residual_f=0;
  // Since multiple tracker threads will be executing this code, it's
  // possible that between the read and the write, a different thread will
  // perform a write. This isn't a problem because the worst that will happen
  // is that we will lose one of many (millions?) of updates.
  global_thread_data.frequency_offset((
    global_thread_data.frequency_offset()*(1/.000001)+
    (frequency_offset+residual_f)*(1/residual_f_np)
  )/(1/.000001+1/residual_f_np));
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
  // This is the only thread that can update the frame timing. Reads and
  // writes to frame_timing are automatically locked by the class.
  double diff=WRAP((rs_curr.frame_timing+delay)-tracked_cell.frame_timing(),-19200.0/2,19200.0/2);
  diff=(0*(1/.0001)+diff*(1/delay_np))/(1/.0001+1/delay_np);
  tracked_cell.frame_timing(itpp_ext::matlab_mod(tracked_cell.frame_timing()+diff,19200.0));
  //cout << "TO: " << setprecision(15) << tracked_cell.frame_timing << endl;
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
  // This is the only thread that can update the frame timing. Reads and
  // writes to frame_timing are automatically locked by the class.
  double diff=WRAP((rs_curr.frame_timing+delay)-tracked_cell.frame_timing(),-19200.0/2,19200.0/2);
  diff=(0*(1/.0001)+diff*(1/delay_np))/(1/.0001+1/delay_np);
  tracked_cell.frame_timing(itpp_ext::matlab_mod(tracked_cell.frame_timing()+diff,19200.0));
  //cout << "TO: " << setprecision(15) << tracked_cell.frame_timing << endl;
}

// Update the frequency domain autocorrelation function.
void do_ac_fd(
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
    boost::mutex::scoped_lock lock(tracked_cell.meas_mutex);
    tracked_cell.ac_fd=elem_div(tracked_cell.ac_fd*(1/.00001)+elem_mult(ac_fd,to_cvec(1.0/ac_fd_np)),to_cvec(1/.00001+1.0/ac_fd_np));
  }
}

// Estimate the time domain autocorrelation function.
void do_ac_td(
  tracked_cell_t & tracked_cell,
  const ce_raw_fifo_pdu_t & rs_curr,
  const double & rs_curr_sp,
  deque <cvec> & ce_history
) {
  // Fill historical fifo
  ce_history.push_back(rs_curr.ce);
  if (ce_history.size()>72) {
    ce_history.pop_front();
  }

  // Calculate autocorrelation if fifo is full.
  if (ce_history.size()==72) {
    cvec this_xc(72);
#ifndef NDEBUG
    this_xc=NAN;
#endif
    for (uint8 t=0;t<72;t++) {
      this_xc(t)=elem_mult_sum(conj(ce_history[71]),ce_history[71-t])/12;
    }
    this_xc=this_xc/rs_curr_sp;

    // Update average
    boost::mutex::scoped_lock lock(tracked_cell.meas_mutex);
    tracked_cell.ac_td=(tracked_cell.ac_td*(1/.00001)+this_xc*1/1)/(1/.00001+1);
  }
}

void interp72(
  const ce_filt_fifo_pdu_t & rs,
  cvec & interp
) {
  interp.set_size(72);
  uint8 l_x=rs.shift;
  complex <double> l_y=rs.ce_filt(0);
  uint8 r_x=rs.shift+6;
  complex <double> r_y=rs.ce_filt(1);
  uint8 ptr=1;
  for (uint8 t=0;t<72;t++) {
    // Advance points if necessary
    if ((t>r_x)&&(ptr<11)) {
      l_x=r_x;
      l_y=r_y;
      r_x+=6;
      ptr++;
      r_y=rs.ce_filt(ptr);
    }
    interp(t)=(r_y-l_y)/(r_x-l_x)*(t-l_x)+l_y;
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
  cvec rs_prev_interp;
  interp72(rs_prev,rs_prev_interp);
  cvec rs_curr_interp;
  interp72(rs_curr,rs_curr_interp);

  // Interpolate in the time domain and push onto FIFO
  uint8 slot_num=rs_prev.slot_num;
  uint8 sym_num=rs_prev.sym_num;
  // Time difference between the current and previous channel estimates.
  double time_diff;
  if (port_num>2) {
    time_diff=0.0005;
  } else {
    if (tracked_cell.cp_type==cp_type_t::EXTENDED) {
      time_diff=3*(128+32)*(1/(FS_LTE/16));
    } else {
      if (rs_prev.sym_num==0) {
        time_diff=4*(128+9)*(1/(FS_LTE/16));
      } else {
        time_diff=(2*(128+9)+(128+10))*(1/(FS_LTE/16));
      }
    }
    //time_diff=time_diff*(1/(FS_LTE/16));
  }

  double time_offset=0;
  while ((slot_num!=rs_curr.slot_num)||(sym_num!=rs_curr.sym_num)) {
    // Interpolate in the time domain.
    cvec rs_mid=rs_prev_interp+(rs_curr_interp-rs_prev_interp)*(time_offset/time_diff);
    double rs_mid_tp=rs_prev.tp+(rs_curr.tp-rs_prev.tp)*(time_offset/time_diff);
    double rs_mid_sp=rs_prev.sp+(rs_curr.sp-rs_prev.sp)*(time_offset/time_diff);
    double rs_mid_sp_raw=rs_prev.sp_raw+(rs_curr.sp_raw-rs_prev.sp_raw)*(time_offset/time_diff);
    double rs_mid_np=rs_prev.np+(rs_curr.np-rs_prev.np)*(time_offset/time_diff);

    // Push onto the interpolated CE fifo.
    ce_interp_fifo_pdu_t pdu;
    pdu.ce_interp=rs_mid;
    pdu.tp=rs_mid_tp;
    pdu.sp=rs_mid_sp;
    pdu.sp_raw=rs_mid_sp_raw;
    pdu.np=rs_mid_np;
    if (!ce_interp_fifo_initialized) {
      // Repeat the very first channel estimates so as to provide CE for
      // slot 0 sym 0.
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
  const bvec & scr,
  deque <mib_fifo_pdu_t> & mib_fifo,
  bool & mib_fifo_synchronized
) {
  //static int mib_successes=0;

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
    //ABORT(-1);
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
      //ABORT(-1);
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
    //ABORT(-1);
    //cout << "syms_mib: " << syms_mib << endl;

    // Extract the bits from the complex modulated symbols.
    vec e_est=lte_demodulate(syms_mib,np_mib,modulation_t::QAM);
    // Unscramble
    //bvec scr=lte_pn(tracked_cell.n_id_cell,length(e_est));
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

    // Unpack MIB information bits that are used to determine whether
    // we are still locked on or not. This reduces the probability of
    // falsely detecting noise as an MIB.
    ivec c_est_ivec=to_ivec(c_est);
    const uint8 bw_packed=c_est_ivec(0)*4+c_est_ivec(1)*2+c_est_ivec(2);
    uint8 n_rb_dl_est=0;
    switch (bw_packed) {
      case 0:
        n_rb_dl_est=6;
        break;
      case 1:
        n_rb_dl_est=15;
        break;
      case 2:
        n_rb_dl_est=25;
        break;
      case 3:
        n_rb_dl_est=50;
        break;
      case 4:
        n_rb_dl_est=75;
        break;
      case 5:
        n_rb_dl_est=100;
        break;
    }

    // PHICH duration
    phich_duration_t::phich_duration_t phich_duration_est=c_est_ivec(3)?phich_duration_t::EXTENDED:phich_duration_t::NORMAL;
    // PHICH resources
    uint8 phich_res=c_est_ivec(4)*2+c_est_ivec(5);
    phich_resource_t::phich_resource_t phich_resource_est;
    switch (phich_res) {
      case 0:
        phich_resource_est=phich_resource_t::oneSixth;
        break;
      case 1:
        phich_resource_est=phich_resource_t::half;
        break;
      case 2:
        phich_resource_est=phich_resource_t::one;
        break;
      case 3:
      default:
        phich_resource_est=phich_resource_t::two;
        break;
    }

    // Did we find it?
    if (
      (crc_est==c_est(24,-1)) &&
      (n_rb_dl_est==tracked_cell.n_rb_dl) &&
      (phich_duration_est==tracked_cell.phich_duration) &&
      (phich_resource_est==tracked_cell.phich_resource)
    ) {
      mib_fifo_synchronized=1;
      {
        boost::mutex::scoped_lock lock(tracked_cell.meas_mutex);
        tracked_cell.mib_decode_failures=0;
      }
      for (uint8 t=0;t<16;t++) {
        mib_fifo.pop_front();
      }
    } else {
      if (mib_fifo_synchronized) {
        {
          boost::mutex::scoped_lock lock(tracked_cell.meas_mutex);
          tracked_cell.mib_decode_failures++;
        }
        for (uint8 t=0;t<16;t++) {
          mib_fifo.pop_front();
        }
      } else {
        {
          boost::mutex::scoped_lock lock(tracked_cell.meas_mutex);
          tracked_cell.mib_decode_failures+=0.25;
        }
        for (uint8 t=0;t<4;t++) {
          mib_fifo.pop_front();
        }
      }
    }

    // 10ms of time increases mib_fifo_decode_failures by 0.25.
    // After several seconds of MIB decoding failures, drop the cell.
    if (tracked_cell.mib_decode_failures>=CELL_DROP_THRESHOLD) {
      //cout << "Dropped a cell!" << endl;
      //boost::mutex::scoped_lock lock(tracked_cell.mutex);
      tracked_cell.kill_me=true;
      return -1;
    }
    //cout << "MIB decode finished" << endl;
  }

  return 0;
}

// Measure signal and noise power using the SSS and PSS. More accurate
// than using the CRS.
// Also perform channel estimation and store results.
void do_pss_sss_sigpower_ce(
  tracked_cell_t & tracked_cell,
  const cvec & syms,
  const uint8 & slot_num,
  const uint8 & sym_num,
  cvec & sss_sym
) {
  // Only examine PSS and SSS symbols
  if (((slot_num!=0)&&(slot_num!=10))||((sym_num!=tracked_cell.n_symb_dl()-2)&&(sym_num!=tracked_cell.n_symb_dl()-1))) {
    return;
  }

  // Store the SSS symbol for when the PSS symbol arrives.
  if (sym_num==tracked_cell.n_symb_dl()-2) {
    sss_sym=syms;
    return;
  }

  const cvec & pss_sym=syms;

  // Measure noise power on 'unoccupied' subcarriers.
  const double np_blank=(sigpower(sss_sym(0,4))+sigpower(sss_sym(67,71))+sigpower(pss_sym(0,4))+sigpower(pss_sym(67,71)))/4;

  // Rotate back.
  cvec ce_sss_raw=elem_mult(sss_sym(5,66),to_cvec(ROM_TABLES.sss_fd(tracked_cell.n_id_1,tracked_cell.n_id_2,(slot_num==0)?0:1)));
  cvec ce_pss_raw=elem_mult(pss_sym(5,66),conj(ROM_TABLES.pss_fd[tracked_cell.n_id_2]));

  // Smoothing
  cvec ce_smooth(62);
  for (uint8 t=0;t<62;t++) {
    uint8 lt=MAX(0,t-6);
    uint8 rt=MIN(t+6,61);
    ce_smooth(t)=(sum(ce_sss_raw(lt,rt))+sum(ce_pss_raw(lt,rt)))/(2*(rt-lt+1));
    //ce_smooth(t)=(sum(ce_sss_raw(lt,rt)))/(1*(rt-lt+1));
  }

  // Measure SP and NP
  // Note correction for estimation bias.
  double np=(sigpower(ce_smooth-ce_sss_raw)*13/12+sigpower(ce_smooth-ce_pss_raw)*13/12)/2;
  //double np=(sigpower(concat(ce_smooth-ce_sss_raw,ce_smooth-ce_pss_raw))*13/12)/1;
  //double np=(sigpower(ce_smooth-ce_sss_raw)*13/12)/1;
  double tp=sigpower(ce_smooth);
  // Note that this value can be negative! The expected value of this
  // measurement, however, is positive.
  double sp=tp-np/13;

  // Store results
  {
    boost::mutex::scoped_lock lock(tracked_cell.meas_mutex);
    tracked_cell.sync_tp=tp;
    tracked_cell.sync_sp=sp;
    tracked_cell.sync_np=np;
    tracked_cell.sync_np_blank=np_blank;
    tracked_cell.sync_ce=concat(zeros_c(5),ce_smooth,zeros_c(5));
    if (isnan(tracked_cell.sync_sp_av)) {
      tracked_cell.sync_tp_av=tp;
      tracked_cell.sync_sp_av=sp;
      tracked_cell.sync_np_av=np;
      tracked_cell.sync_np_blank_av=np_blank;
    } else {
      tracked_cell.sync_tp_av=0.999*tracked_cell.sync_tp_av+.001*tp;
      tracked_cell.sync_sp_av=0.999*tracked_cell.sync_sp_av+.001*sp;
      tracked_cell.sync_np_av=0.999*tracked_cell.sync_np_av+.001*np;
      tracked_cell.sync_np_blank_av=0.999*tracked_cell.sync_np_blank_av+.001*np_blank;
    }
  }
}

// Process that tracks a cell that has been found by the searcher.
void tracker_thread(
  tracked_cell_t & tracked_cell,
  global_thread_data_t & global_thread_data
) {
  // Pre-compute some information.
  //ivec cn=concat(itpp_ext::matlab_range(-36,-1),itpp_ext::matlab_range(1,36));
  // Reference symbols
  RS_DL rs_dl(tracked_cell.n_id_cell,6,tracked_cell.cp_type);
  // MIB scrambling sequence.
  const bvec scr=lte_pn(tracked_cell.n_id_cell,(tracked_cell.cp_type==cp_type_t::NORMAL)?1920:1728);

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
  cvec sss_sym;
  //double mib_fifo_decode_failures=0;
  // Store the channel estimates so that the time domain channel
  // autocorrelation function can be estimated.
  vector <deque <cvec> > ce_history(tracked_cell.n_ports);
  // Now that everything has been initialized, indicate to the producer
  // thread that we are ready for data.
  // Data flow starts at the beginning of the next frame.
  tracked_cell.tracker_thread_ready=true;
  // Each iteration of this loop processes one OFDM symbol.
  while (true) {
    // If there is more than 1.5s worth of data in the fifo, dump
    // data to allow the tracker threads to catch up.
    {
      boost::mutex::scoped_lock lock(tracked_cell.fifo_mutex);
      uint16 n_ofdm_1s=(tracked_cell.cp_type==cp_type_t::NORMAL)?(7*2*1000):(6*2*1000);
      while (tracked_cell.fifo.size()>n_ofdm_1s*1.5) {
        for (uint32 t=0;t<n_ofdm_1s;t++) {
          tracked_cell.fifo.pop();
        }
        global_thread_data.cell_seconds_dropped_inc();
      }
    }

    // Get the next frequency domain sample from the fifo.
    cvec syms;
    double frequency_offset;
    double frame_timing;
    //get_fd(tracked_cell,global_thread_data.fc,slot_num,sym_num,cn,bulk_phase_offset,syms,frequency_offset,frame_timing);
    get_fd(tracked_cell,global_thread_data.fc_requested,global_thread_data.fc_programmed,global_thread_data.fs_programmed,slot_num,sym_num,bulk_phase_offset,syms,frequency_offset,frame_timing);

    // Save this information into the data fifo for further processing
    // once channel estimates are ready. Channel estimates for this OFDM
    // symbol may not be ready until several more OFDM symbols have been
    // received.
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
      // Note correction for the estimation bias.
      const double rs_curr_np=sigpower(rs_curr.ce-rs_curr_filt)*7/6;
      //const double rs_curr_np=sigpower(rs_curr.ce-rs_curr_filt)*(1.0/(pow(2.0/7.0,2.0)*2.0+pow(3.0/7.0,2.0)*2.0/3.0));
      // Note that this value can be negative!
      //const double rs_curr_sp=MAX(sigpower(rs_curr_filt)-rs_curr_np/7,0);
      const double rs_curr_tp=sigpower(rs_curr_filt);
      const double rs_curr_sp_raw=rs_curr_tp-rs_curr_np/7;
      const double rs_curr_sp=MAX(.00001,rs_curr_sp_raw);
      //cout << "SP RSC " << db10(sigpower(rs_curr.ce)) << endl;
      //cout << "SP/NP " << db10(rs_curr_sp) << " / " << db10(rs_curr_np) << endl;
      // Store filtered channel estimates.
      ce_filt_fifo_pdu_t pdu;
      pdu.shift=rs_curr.shift;
      pdu.slot_num=rs_curr.slot_num;
      pdu.sym_num=rs_curr.sym_num;
      pdu.tp=rs_curr_tp;
      pdu.sp=rs_curr_sp;
      pdu.sp_raw=rs_curr_sp_raw;
      pdu.np=rs_curr_np;
      pdu.ce_filt=rs_curr_filt;
      ce_filt_fifo[port_num].push_back(pdu);

      // FOE
      do_foe(global_thread_data,rs_prev,rs_next,rs_curr_np,rs_curr_filt);

      // TOE
      //do_toe(tracked_cell,rs_curr,rs_curr_filt,rs_curr_np);
      do_toe_v2(tracked_cell,rs_prev,rs_curr,rs_curr_sp,rs_curr_np);

      // Estimate frequency domain autocorrelations.
      do_ac_fd(tracked_cell,rs_curr,rs_curr_sp,rs_curr_np);

      // Estimate the time domain autocorrelation function.
      do_ac_td(tracked_cell,rs_curr,rs_curr_sp,ce_history[port_num]);

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
#ifndef NDEBUG
      // Synchronization check.
      for (uint8 t=0;t<tracked_cell.n_ports;t++) {
        if (
          (data_fifo.front().slot_num!=ce_interp_fifo[t].front().slot_num)||
          (data_fifo.front().sym_num!=ce_interp_fifo[t].front().sym_num)
        ) {
          cerr << "Error: synchronization error! Check code!" << endl;
          ABORT(-1);
        }
      }
#endif

      // For this OFDM symbol, extract the symbols, the channel estimates,
      // signal power, etc.
      cvec & syms=data_fifo.front().syms;
      cmat ce(tracked_cell.n_ports,72);
      vec tp(tracked_cell.n_ports);
      vec sp(tracked_cell.n_ports);
      vec sp_raw(tracked_cell.n_ports);
      vec np(tracked_cell.n_ports);
      uint8 data_slot_num=data_fifo.front().slot_num;
      uint8 data_sym_num=data_fifo.front().sym_num;
      for (uint8 t=0;t<tracked_cell.n_ports;t++) {
        ce.set_row(t,ce_interp_fifo[t].front().ce_interp);
        tp(t)=ce_interp_fifo[t].front().tp;
        sp(t)=ce_interp_fifo[t].front().sp;
        sp_raw(t)=ce_interp_fifo[t].front().sp_raw;
        np(t)=ce_interp_fifo[t].front().np;
      }

      // Store channel estimates
      {
        boost::mutex::scoped_lock lock(tracked_cell.meas_mutex);
        tracked_cell.ce=ce;
      }

      // Store signal power measurements.
      {
        boost::mutex::scoped_lock lock(tracked_cell.meas_mutex);
        tracked_cell.crs_sp_raw=sp_raw;
        tracked_cell.crs_np=np;
        if (isnan(tracked_cell.crs_sp_raw_av(0))) {
          tracked_cell.crs_tp_av=tp;
          tracked_cell.crs_sp_raw_av=sp_raw;
          tracked_cell.crs_np_av=np;
        } else {
          if (0) {
            tracked_cell.crs_tp_av=0.99999*tracked_cell.crs_tp_av+.00001*tp;
            tracked_cell.crs_sp_raw_av=0.99999*tracked_cell.crs_sp_raw_av+.00001*sp_raw;
            tracked_cell.crs_np_av=0.99999*tracked_cell.crs_np_av+.00001*np;
          } else {
            // This code only averages the measurements for PSS and SSS ofdm
            // symbols.
            if (((data_slot_num==0)||(data_slot_num==10))&&((data_sym_num==5)||(data_sym_num==6))) {
              tracked_cell.crs_tp_av=0.999*tracked_cell.crs_tp_av+.001*tp;
              tracked_cell.crs_sp_raw_av=0.999*tracked_cell.crs_sp_raw_av+.001*sp_raw;
              tracked_cell.crs_np_av=0.999*tracked_cell.crs_np_av+.001*np;
            }
          }
        }
      }

      // Measure signal power and noise power on PSS/SSS (more accurate)
      do_pss_sss_sigpower_ce(tracked_cell,syms,data_slot_num,data_sym_num,sss_sym);

      // Perform MIB decoding
      if (do_mib_decode(tracked_cell,syms,ce,sp,np,data_slot_num,data_sym_num,scr,mib_fifo,mib_fifo_synchronized)==-1) {
        // We have failed to detect an MIB for a long time. Exit this
        // thread.
        //cout << "Tracker thread exiting..." << endl;
        return;
      }

      // Done processing data. Pop data vector and CE vectors.
      data_fifo.pop_front();
      for (uint8 t=0;t<tracked_cell.n_ports;t++) {
        ce_interp_fifo[t].pop_front();
      }
    }

    // Increase the local counter.
    slot_sym_inc(tracked_cell.n_symb_dl(),slot_num,sym_num);
  }
}

