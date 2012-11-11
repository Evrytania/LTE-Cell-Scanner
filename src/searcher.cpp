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

// Relationships between the various frequencies, correction factors, etc.
//
// Fixed relationships:
// xtal_spec=28.8MHz (usually...)
// k_factor=xtal_true/xtal_spec
// xtal_true=xtal_spec*k_factor
// fs_true=k_s*xtal_true
// fs_prog=k_s*xtal_spec
// fc_true=k_c*xtal_true
// fc_prog=k_c*xtal_spec
// freq_offset=fc_req-fc_true
//
// k_factor as a function of known parameters:
// fc_true=fc_req-freq_offset
// xtal_true=fc_true/k_c=(fc_req-freq_offset)/k_c
// k_factor=(fc_req-freq_offset)/(k_c*xtal_spec)
// k_factor=(fc_req-freq_offset)/(fc_prog/xtal_spec*xtal_spec)
// k_factor=(fc_req-freq_offset)/fc_prog
//
// fs_prog*k_factor=k_s*xtal_spec*k_factor=k_s*xtal_true=fs_true
//
// N LTE samples sampled at a rate of FS_LTE/16 is equivalent to this
// many samples sampled at fs_true.
// N*(1/(FS_LTE/16))/(1/(fs_prog*k_factor))
// N*(16/FS_LTE)/(1/(fs_prog*k_factor))
// N*16/FS_LTE*fs_prog*k_factor

#include <itpp/itbase.h>
#include <itpp/signal/transforms.h>
#include <itpp/stat/misc_stat.h>
#include <math.h>
#include <list>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <boost/math/special_functions/gamma.hpp>
#include "rtl-sdr.h"
#include "common.h"
#include "lte_lib.h"
#include "constants.h"
#include "macros.h"
#include "itpp_ext.h"
#include "dsp.h"
#include "searcher.h"

// This LTE cell search algorithm was designed under the following guidelines:
//
// 1) Low performance hardware
// This algorithm was designed to work with the RTL-SDR dongle which has a
// noise figure of 20dB. As such, attempts were made to maximize signal
// processing performance on the SDR side so as to make up for some of the
// peroformance lost in the hardware. This was done at the cost of algorithm
// complexity and execution speed.
//
// 2) Capture-then-process
// The algorithm was designed so that data would first be captured and then
// analysis would begin. It is not possible for the analysis algorithm to
// request 'more' data part-way through the analysis step.
//
// For maximal signal processing performance, the entire MIB must be captured.
// The MIB spans 40ms but the location of the start of the MIB is not knwon
// in advance. Thus, to ensure that an entire MIB is captured, approximately
// 80ms of data must be captured.
//
// 3) Use all available data
// For example, PSS detection only needs 5ms of data. However, since 80ms
// of data will be available, all 80ms of captured data will be used to
// improve PSS detection performance.
//
// 4) Handle all cell types
// This algorithm can handle synchronous or asynchronous networks and can
// also be used in situations where the receiver is receiving a signal
// from both normal and extended CP cells. All found cells will be reported.
//
// 5) Work for all cell loads
// In a highly loaded cell, the transmitted (and received) power is more
// or less constant in time and hence the entire received signal can be used
// to estimate the total received power. In a lightly loaded cell, the TX
// power varies wildly between OFDM symbols containing PSS/SSS, OFDM symbols
// containing RS, and OFDM symbols with no data being transmitted.
//
// Final performance:
// ==================
// Simulations indicate that this algorithm can reliably detect the PSS/SSS
// and therefore the Cell ID down to SNR's (AWGN) of about -12dB. MIB decoding
// is only reliable down to about -10dB and thus overal performance is limited
// by MIB decoding.

using namespace itpp;
using namespace std;

// Correlate the received data against various frequency shifted versions
// of the three PSS sequences.
// This is likely to be the slowest routine since it needs to process so
// much data.
void xc_correlate(
  // Inputs
  const cvec & capbuf,
  const vec & f_search_set,
  const double & fc_requested,
  const double & fc_programmed,
  const double & fs_programmed,
  // Outputs
  vcf3d & xc
) {
  const uint32 n_cap=length(capbuf);
  const uint16 n_f=length(f_search_set);

  // Set aside space for the vector and initialize with NAN's.
#ifndef NDEBUG
  xc=vector < vector < vector < complex < float > > > > (3,vector< vector < complex < float > > >(n_cap-136, vector < complex < float > > (n_f,NAN)));
#else
  xc=vector < vector < vector < complex < float > > > > (3,vector< vector < complex < float > > >(n_cap-136, vector < complex < float > > (n_f)));
#endif

  // Local variables declared outside of the loop.
  double f_off;
  cvec temp;
  complex <double> acc;
  uint16 foi;
  uint8 t;
  uint32 k;
  uint8 m;

  // Loop and perform correlations.
  //Real_Timer tt;
  //tt.tic();
  for (foi=0;foi<n_f;foi++) {
    f_off=f_search_set(foi);
    double k_factor=(fc_requested-f_off)/fc_programmed;
    for (t=0;t<3;t++) {
      temp=ROM_TABLES.pss_td[t];
      temp=fshift(temp,f_off,fs_programmed*k_factor);
      temp=conj(temp)/137;
#ifdef _OPENMP
#pragma omp parallel for shared(temp,capbuf,xc) private(k,acc,m)
#endif
      for (k=0;k<n_cap-136;k++) {
        acc=0;
        for (m=0;m<137;m++) {
          // Correlations are performed at the 2x rate which effectively
          // performs filtering and correlating at the same time. Thus,
          // this algorithm can handle huge frequency offsets limited only
          // by the bandwidth of the capture device.
          // Correlations can also be done at the 1x rate if filtering is
          // peformed first, but this will limit the set of frequency offsets
          // that this algorithm can detect. 1x rate correlations will,
          // however, be nearly twice as fast as the 2x correlations
          // performed here.
          acc+=temp(m)*capbuf(k+m);
        }
        xc[t][k][foi]=acc;
      }
    }
  }
  //tt.toc_print();
}

// Estimate the received signal power within 2 OFDM symbols of a particular
// sample.
//
// In the 6 center RB's, the transmitted power is the same for all PSS and
// SSS OFDM symbols regardless of the cell load.
//
// This function is slightly inaccurate because it estimates the received
// power in all the RB's (approximately 12) instead of the signal power
// only in the center 6 RB's.
void sp_est(
  // Inputs
  const cvec & capbuf,
  // Outputs
  vec & sp,
  vec & sp_incoherent,
  uint16 & n_comb_sp
) {
  const uint32 n_cap=length(capbuf);
  n_comb_sp=floor_i((n_cap-136-137)/9600);
  const uint32 n_sp=n_comb_sp*9600;

  // Set aside space for the vector and initialize with NAN's.
  sp=vec(n_sp);
#ifndef NDEBUG
  sp=NAN;
#endif
  sp[0]=0;
  // Estimate power for first time offset
  for (uint16 t=0;t<274;t++) {
    sp[0]+=pow(capbuf[t].real(),2)+pow(capbuf[t].imag(),2);
  }
  sp[0]=sp[0]/274;
  // Estimate RX power for remaining time offsets.
  for (uint32 t=1;t<n_sp;t++) {
    sp[t]=sp[t-1]+(-pow(capbuf[t-1].real(),2)-pow(capbuf[t-1].imag(),2)+pow(capbuf[t+274-1].real(),2)+pow(capbuf[t+274-1].imag(),2))/274;
  }

  // Combine incoherently
  sp_incoherent=sp.left(9600);
  for (uint16 t=1;t<n_comb_sp;t++) {
    sp_incoherent+=sp.mid(t*9600,9600);
  }
  sp_incoherent=sp_incoherent/n_comb_sp;
  // Shift to the right by 137 samples to align with the correlation peaks.
  tshift(sp_incoherent,137);
}

// Perform incoherent combining on the correlations.
//
// There is no guarantee that PSS/SSS pairs will always be transmitted from
// the same antenna. Thus, PSS/SSS pairs separated by 5ms can only be combined
// incoherently.
//
// Because the size of the capture buffer is very large (80ms) and because this
// algorithm must be able to handle very large frequency offsets (50kHz+),
// care must be taken to combine the correct samples incoherently.
//
// With no frequency error, the start of the next frame will be 19200 samples
// after the start of the current frame. With a +50kHz frequency offset and
// a 740MHz center frequency, the start of the next frame is actually 19198.7
// samples after the start of this frame. This function makes sure that
// the correct samples are combined incoherently.
//
// A side benefit of capturing multiple frames is that the correct downlink
// center frequency can be inferred. For an extreme example, suppose that
// the true downlink center frequency was 740MHz and that the frequency error
// of the oscillator could be up to 100kHz, but currently happens, by chance,
// to be zero. (We do not know in advance that the frequency error is zero.)
//
// For each programmed center frequency, the search will search frequency
// offsets from -100kHz to +100kHz. The searcher will program a center
// frequency of 739.9 MHz and search from -100kHz to +100kHz and it will
// then program a frequency of 740MHz and search from -100kHz to +100kHz,
// and it will also program a frequency of 740.1MHz and search from -100kHz
// to +100kHz.
//
// The correlations performed at 739.9MHz+100kHz, 740MHz+0kHz, and 740.1-100kHz
// will all be nearly the same and it will not be possible for the receiver
// to determine whether the true center frequency was 739.9, 740, or 740.1
// MHz.
//
// Although the correlation peaks for all 3 of the above scenarios will
// have the same magnitude, the spacing between the peaks will vary if
// multiple frames are captured. Only in the 740Mhz+0kHz case will the peaks
// be aligned with a spacing of 19200 samples and thus it is possible
// to determine boht that the true downlink center frequency is 740MHz and
// that the frequency error of the local oscillator is 0Hz.
void xc_combine(
  // Inputs
  const cvec & capbuf,
  const vcf3d & xc,
  const double & fc_requested,
  const double & fc_programmed,
  const double & fs_programmed,
  const vec & f_search_set,
  // Outputs
  vf3d & xc_incoherent_single,
  uint16 & n_comb_xc
) {
  const uint16 n_f=f_search_set.length();
  n_comb_xc=floor_i((xc[0].size()-100)/9600);

  // Create space for some arrays
#ifndef NDEBUG
  xc_incoherent_single=vector < vector < vector < float > > > (3,vector< vector < float > >(9600, vector < float > (n_f,NAN)));
#else
  xc_incoherent_single=vector < vector < vector < float > > > (3,vector< vector < float > >(9600, vector < float > (n_f)));
#endif
  for (uint16 foi=0;foi<n_f;foi++) {
    // Combine incoherently
    const double f_off=f_search_set[foi];
    const double k_factor=(fc_requested-f_off)/fc_programmed;
    for (uint8 t=0;t<3;t++) {
      for (uint16 idx=0;idx<9600;idx++) {
        xc_incoherent_single[t][idx][foi]=0;
      }
      for (uint16 m=0;m<n_comb_xc;m++) {
        // Because of the large supported frequency offsets and the large
        // amount of time represented by the capture buffer, the length
        // in samples, of a frame varies by the frequency offset.
        //double actual_time_offset=m*.005*k_factor;
        //double actual_start_index=itpp::round_i(actual_time_offset*FS_LTE/16);
        double actual_start_index=itpp::round_i(m*.005*k_factor*fs_programmed);
        for (uint16 idx=0;idx<9600;idx++) {
          xc_incoherent_single[t][idx][foi]+=sqr(xc[t][idx+actual_start_index][foi]);
        }
      }
      for (uint16 idx=0;idx<9600;idx++) {
        xc_incoherent_single[t][idx][foi]=xc_incoherent_single[t][idx][foi]/n_comb_xc;
      }
    }
  }
}

// Combine adjacent taps that likely come from the same channel.
// Simply: xc_incoherent(t,idx,foi)=mean(xc_incoherent_single(t,idx-ds_comb_arm:idx+ds_comb_arm,foi);
void xc_delay_spread(
  // Inputs
  const vf3d & xc_incoherent_single,
  const uint8 & ds_comb_arm,
  // Outputs
  vf3d & xc_incoherent
) {
  const int n_f=xc_incoherent_single[0][0].size();

  // Create space for some arrays
#ifndef NDEBUG
  xc_incoherent=vector < vector < vector < float > > > (3,vector< vector < float > >(9600, vector < float > (n_f,NAN)));
#else
  xc_incoherent=vector < vector < vector < float > > > (3,vector< vector < float > >(9600, vector < float > (n_f)));
#endif
  for (uint16 foi=0;foi<n_f;foi++) {
    for (uint8 t=0;t<3;t++) {
      for (uint16 idx=0;idx<9600;idx++) {
        xc_incoherent[t][idx][foi]=xc_incoherent_single[t][idx][foi];
      }
    }
    for (uint8 t=1;t<=ds_comb_arm;t++) {
      for (uint8 k=0;k<3;k++) {
        for (uint16 idx=0;idx<9600;idx++) {
          xc_incoherent[k][idx][foi]+=xc_incoherent_single[k][itpp_ext::matlab_mod(idx-t,9600)][foi]+xc_incoherent_single[k][itpp_ext::matlab_mod(idx+t,9600)][foi];
        }
      }
    }
    // Normalize
    for (uint8 t=0;t<3;t++) {
      for (uint16 idx=0;idx<9600;idx++) {
        xc_incoherent[t][idx][foi]=xc_incoherent[t][idx][foi]/(2*ds_comb_arm+1);
      }
    }
  }
}

// Search for the peak correlation among all frequency offsets.
// For each time offset and each PSS index, examine the correlations
// for all of the frequency offsets and only keep the correlation with
// the largest magnitude.
void xc_peak_freq(
  // Inputs
  const vf3d & xc_incoherent,
  // Outputs
  mat & xc_incoherent_collapsed_pow,
  imat & xc_incoherent_collapsed_frq
) {
  const int n_f=xc_incoherent[0][0].size();

  xc_incoherent_collapsed_pow=mat(3,9600);
  xc_incoherent_collapsed_frq=imat(3,9600);
#ifndef NDEBUG
  xc_incoherent_collapsed_pow=NAN;
  xc_incoherent_collapsed_frq=-1;
#endif

  for (uint8 t=0;t<3;t++) {
    for (uint16 k=0;k<9600;k++) {
      double best_pow=xc_incoherent[t][k][0];
      uint16 best_idx=0;
      for (uint16 foi=1;foi<n_f;foi++) {
        if (xc_incoherent[t][k][foi]>best_pow) {
          best_pow=xc_incoherent[t][k][foi];
          best_idx=foi;
        }
      }
      xc_incoherent_collapsed_pow(t,k)=best_pow;
      xc_incoherent_collapsed_frq(t,k)=best_idx;
    }
  }
}

// Correlate the received signal against all possible PSS and all possible
// frequency offsets.
// This is the main function that calls all of the previously declared
// subfunctions.
void xcorr_pss(
  // Inputs
  const cvec & capbuf,
  const vec & f_search_set,
  const uint8 & ds_comb_arm,
  const double & fc_requested,
  const double & fc_programmed,
  const double & fs_programmed,
  // Outputs
  mat & xc_incoherent_collapsed_pow,
  imat & xc_incoherent_collapsed_frq,
  // Following used only for debugging...
  vf3d & xc_incoherent_single,
  vf3d & xc_incoherent,
  vec & sp_incoherent,
  vcf3d & xc,
  vec & sp,
  uint16 & n_comb_xc,
  uint16 & n_comb_sp
) {
  // Perform correlations
  xc_correlate(capbuf,f_search_set,fc_requested,fc_programmed,fs_programmed,xc);
  // Incoherently combine correlations
  xc_combine(capbuf,xc,fc_requested,fc_programmed,fs_programmed,f_search_set,xc_incoherent_single,n_comb_xc);
  // Combine according to delay spread
  xc_delay_spread(xc_incoherent_single,ds_comb_arm,xc_incoherent);
  // Estimate received signal power
  sp_est(capbuf,sp,sp_incoherent,n_comb_sp);
  // Search for peaks among all the frequency offsets.
  xc_peak_freq(xc_incoherent,xc_incoherent_collapsed_pow,xc_incoherent_collapsed_frq);
}

// Search through all the correlations and determine if any PSS were found.
void peak_search(
  // Inputs
  const mat & xc_incoherent_collapsed_pow,
  const imat & xc_incoherent_collapsed_frq,
  const vec & Z_th1,
  const vec & f_search_set,
  const double & fc_requested,
  const double & fc_programmed,
  const vf3d & xc_incoherent_single,
  const uint8 & ds_comb_arm,
  // Outputs
  list <Cell> & cells
) {
  // Create local copy we can write to and destroy.
  mat xc_incoherent_working=xc_incoherent_collapsed_pow;

  for (;;) {
    // Search for the largest peak. (Not the largest peak relative to
    // the detection threshold Z_th1.)
    ivec peak_ind_v;
    vec peak_pow_v=max(transpose(xc_incoherent_working),peak_ind_v);
    int32 peak_n_id_2;
    double peak_pow=max(peak_pow_v,peak_n_id_2);
    int32 peak_ind=peak_ind_v(peak_n_id_2);
    if (peak_pow<Z_th1(peak_ind)) {
      // This peak has too low of a received power. There are no more
      // interesting peaks. Break!
      break;
    }

    // A peak was found at location peak_ind and has frequency index
    // xc_incoherent_collapsed_frq(peak_n_id_2,peak_ind). This peak
    // is the sum of the energy within ds_comb_arm samples around this
    // peak location. From the samples within ds_comb_arm samples
    // around peak_ind, find the index with the highest power.
    double best_pow=-INFINITY;
    int16 best_ind=-1;
    for (uint16 t=peak_ind-ds_comb_arm;t<=peak_ind+ds_comb_arm;t++) {
      uint16 t_wrap=mod(t,9600);
      if (xc_incoherent_single[peak_n_id_2][t_wrap][xc_incoherent_collapsed_frq(peak_n_id_2,peak_ind)]>best_pow) {
        best_pow=xc_incoherent_single[peak_n_id_2][t_wrap][xc_incoherent_collapsed_frq(peak_n_id_2,peak_ind)];
        best_ind=t_wrap;
      }
    }

    // Record this peak for further processing
    Cell cell;
    cell.fc_requested=fc_requested;
    cell.fc_programmed=fc_programmed;
    cell.pss_pow=peak_pow;
    //cell.ind=peak_ind;
    cell.ind=best_ind;
    cell.freq=f_search_set(xc_incoherent_collapsed_frq(peak_n_id_2,peak_ind));
    cell.n_id_2=peak_n_id_2;
    cells.push_back(cell);

    // Cancel out the false peaks around this one.
    // No peaks with the same pss sequence are allowed within 274 samples of
    // this one.
    for (int16 t=-274;t<=274;t++) {
      //cout <<mod(peak_ind+t,9600)<<endl;
      xc_incoherent_working(peak_n_id_2,itpp_ext::matlab_mod(peak_ind+t,9600))=0;
    }
    // Cancel out other PSS sequences whose power is within 8dB of the current
    // sequence.
    double thresh=peak_pow*udb10(-8.0);
    for (uint8 n=0;n<=3;n++) {
      if (n==peak_n_id_2) {
        continue;
      }
      for (int16 t=-274;t<=274;t++) {
        if (xc_incoherent_working(peak_n_id_2,itpp_ext::matlab_mod(peak_ind+t,9600))<thresh) {
          xc_incoherent_working(peak_n_id_2,itpp_ext::matlab_mod(peak_ind+t,9600))=0;
        }
      }
    }
    // Because of the repetitive nature of the CRS, a PSS at offset I with power
    // P will produce correlation peaks during all CRS OFDM symbols with power
    // approximately P-14dB. Cancel them out.
    thresh=peak_pow*udb10(-12.0);
    for (uint8 r=0;r<3;r++) {
      for (uint16 c=0;c<9600;c++) {
        if (xc_incoherent_working(r,c)<thresh) {
          xc_incoherent_working(r,c)=0;
        }
      }
    }
  }
}

// Simple helper function to perform FOC and return only the subcarriers
// occupied by the PSS or SSS.
//
// Called by more than one function!
inline cvec extract_psss(
  const cvec td_samps,
  const double foc_freq,
  const double & k_factor,
  const double & fs_programmed
) {
  // Frequency shift
  cvec dft_in=fshift(td_samps,foc_freq,fs_programmed*k_factor);
  // Remove the 2 sample time offset
  dft_in=concat(dft_in(2,-1),dft_in.left(2));
  // DFT
  cvec dft_out=dft(dft_in);
  // Extract interesting samples.
  return concat(dft_out.right(31),dft_out.mid(1,31));
}

// Perform channel estimation and extract the SSS subcarriers.
void sss_detect_getce_sss(
  // Inputs
  const Cell & cell,
  const cvec & capbuf,
  const double & fc_requested,
  const double & fc_programmed,
  const double & fs_programmed,
  // Outputs
  vec & sss_h1_np_est,
  vec & sss_h2_np_est,
  cvec & sss_h1_nrm_est,
  cvec & sss_h2_nrm_est,
  cvec & sss_h1_ext_est,
  cvec & sss_h2_ext_est
) {
  // Local copies
  double peak_loc=cell.ind;
  const double peak_freq=cell.freq;
  const uint8 n_id_2_est=cell.n_id_2;

  const double k_factor=(fc_requested-peak_freq)/fc_programmed;

  // Skip to the right by 5 subframes if there is no room here to detect
  // the SSS.
  if (peak_loc+9<162) {
    peak_loc+=9600*k_factor;
  }
  // The location of all PSS's in the capture buffer where we also have
  // access to an SSS.
  vec pss_loc_set=itpp_ext::matlab_range(peak_loc,k_factor*9600,(double)capbuf.length()-125-9);
  uint16 n_pss=length(pss_loc_set);
  vec pss_np(n_pss);
  cmat h_raw(n_pss,62);
  cmat h_sm(n_pss,62);
  cmat sss_nrm_raw(n_pss,62);
  cmat sss_ext_raw(n_pss,62);
#ifndef NDEBUG
  pss_np=NAN;
  h_raw=NAN;
  h_sm=NAN;
  sss_nrm_raw=NAN;
  sss_ext_raw=NAN;
#endif

  for (uint16 k=0;k<n_pss;k++) {
    uint32 pss_loc=itpp::round_i(pss_loc_set(k));
    uint32 pss_dft_location=pss_loc+9-2;

    // Calculate channel response
    h_raw.set_row(k,elem_mult(extract_psss(capbuf.mid(pss_dft_location,128),-peak_freq,k_factor,fs_programmed),conj(ROM_TABLES.pss_fd[n_id_2_est])));
    // Basic smoothing. Average nearest 6 subcarriers.
    for (uint8 t=0;t<62;t++) {
      uint8 lt=MAX(0,t-6);
      uint8 rt=MIN(61,t+6);
      h_sm(k,t)=mean(h_raw.get_row(k).mid(lt,rt-lt+1));
    }

    // Estimate noise power
    pss_np(k)=sigpower(h_sm.get_row(k)-h_raw.get_row(k));

    // Calculate SSS in the frequency domain for extended and normal CP
    uint32 sss_dft_location=pss_dft_location-128-32;
    sss_ext_raw.set_row(k,extract_psss(capbuf.mid(sss_dft_location,128),-peak_freq,k_factor,fs_programmed));
    sss_dft_location=pss_dft_location-128-9;
    sss_nrm_raw.set_row(k,extract_psss(capbuf.mid(sss_dft_location,128),-peak_freq,k_factor,fs_programmed));
  }

  // Combine results from different slots
  // Pre-allocate some vectors that are often used
  vec pss_np_inv_h1=1.0/pss_np(itpp_ext::matlab_range(0,2,n_pss-1));
  vec pss_np_inv_h2=1.0/pss_np(itpp_ext::matlab_range(1,2,n_pss-1));
  sss_h1_np_est.set_size(62);
  sss_h2_np_est.set_size(62);
  sss_h1_nrm_est.set_size(62);
  sss_h2_nrm_est.set_size(62);
  sss_h1_ext_est.set_size(62);
  sss_h2_ext_est.set_size(62);
#ifndef NDEBUG
  sss_h1_np_est=NAN;
  sss_h2_np_est=NAN;
  sss_h1_nrm_est=NAN;
  sss_h2_nrm_est=NAN;
  sss_h1_ext_est=NAN;
  sss_h2_ext_est=NAN;
#endif
  for (uint8 t=0;t<62;t++) {
    // First half (h1) and second half (h2) channel estimates.
    cvec h_sm_h1=h_sm.get_col(t).get(itpp_ext::matlab_range(0,2,n_pss-1));
    cvec h_sm_h2=h_sm.get_col(t).get(itpp_ext::matlab_range(1,2,n_pss-1));
    // Estimate noise power in each subcarrier
    sss_h1_np_est(t)=1/(1+sum(elem_mult(sqr(h_sm_h1),pss_np_inv_h1)));
    sss_h2_np_est(t)=1/(1+sum(elem_mult(sqr(h_sm_h2),pss_np_inv_h2)));
    // Estimate SSS assuming normal CP
    sss_h1_nrm_est(t)=sss_h1_np_est(t)*sum(elem_mult(conj(h_sm_h1),to_cvec(pss_np_inv_h1),sss_nrm_raw.get_col(t).get(itpp_ext::matlab_range(0,2,n_pss-1))));
    sss_h2_nrm_est(t)=sss_h2_np_est(t)*sum(elem_mult(conj(h_sm_h2),to_cvec(pss_np_inv_h2),sss_nrm_raw.get_col(t).get(itpp_ext::matlab_range(1,2,n_pss-1))));
    // Estimate SSS assuming extended CP
    sss_h1_ext_est(t)=sss_h1_np_est(t)*sum(elem_mult(conj(h_sm_h1),to_cvec(pss_np_inv_h1),sss_ext_raw.get_col(t).get(itpp_ext::matlab_range(0,2,n_pss-1))));
    sss_h2_ext_est(t)=sss_h2_np_est(t)*sum(elem_mult(conj(h_sm_h2),to_cvec(pss_np_inv_h2),sss_ext_raw.get_col(t).get(itpp_ext::matlab_range(1,2,n_pss-1))));
  }
}

// Helper function that compares the received SSS against one of the
// known transmitted SSS sequences.
double sss_detect_ml_helper(
  const vec & sss_h12_np_est,
  const cvec & sss_h12_est,
  const ivec & sss_h12_try_orig
) {
  cvec sss_h12_try(to_cvec(sss_h12_try_orig));

  // Compensate for phase errors between the est and try sequences
  double ang=arg(sum(elem_mult(conj(sss_h12_est),sss_h12_try)));
  sss_h12_try*=exp(J*-ang);

  // Calculate the log likelihood
  cvec diff=sss_h12_try-sss_h12_est;
  double log_lik=-sum(elem_div(elem_mult(real(diff),real(diff)),sss_h12_np_est))-sum(elem_div(elem_mult(imag(diff),imag(diff)),sss_h12_np_est));

  return log_lik;
}

// Perform maximum likelihood detection on the combined SSS signals.
void sss_detect_ml(
  // Inputs
  const Cell & cell,
  const vec & sss_h1_np_est,
  const vec & sss_h2_np_est,
  const cvec & sss_h1_nrm_est,
  const cvec & sss_h2_nrm_est,
  const cvec & sss_h1_ext_est,
  const cvec & sss_h2_ext_est,
  // Outputs
  mat & log_lik_nrm,
  mat & log_lik_ext
) {
  log_lik_nrm.set_size(168,2);
  log_lik_ext.set_size(168,2);
#ifndef NDEBUG
  log_lik_nrm=NAN;
  log_lik_ext=NAN;
#endif

  vec sss_h12_np_est=concat(sss_h1_np_est,sss_h2_np_est);
  cvec sss_h12_nrm_est=concat(sss_h1_nrm_est,sss_h2_nrm_est);
  cvec sss_h12_ext_est=concat(sss_h1_ext_est,sss_h2_ext_est);
  for (uint8 t=0;t<168;t++) {
    // Construct the SSS sequence that will be compared against the
    // received sequence.
    ivec sss_h1_try=ROM_TABLES.sss_fd(t,cell.n_id_2,0);
    ivec sss_h2_try=ROM_TABLES.sss_fd(t,cell.n_id_2,10);
    ivec sss_h12_try=concat(sss_h1_try,sss_h2_try);
    ivec sss_h21_try=concat(sss_h2_try,sss_h1_try);

    // Calculate log likelihood for normal/extended and 12/21 ordering
    // of SSS sequence.
    log_lik_nrm(t,0)=sss_detect_ml_helper(sss_h12_np_est,sss_h12_nrm_est,sss_h12_try);
    log_lik_nrm(t,1)=sss_detect_ml_helper(sss_h12_np_est,sss_h12_nrm_est,sss_h21_try);
    log_lik_ext(t,0)=sss_detect_ml_helper(sss_h12_np_est,sss_h12_ext_est,sss_h12_try);
    log_lik_ext(t,1)=sss_detect_ml_helper(sss_h12_np_est,sss_h12_ext_est,sss_h21_try);
  }
}

// Detect the SSS, if present
Cell sss_detect(
  // Inputs
  const Cell & cell,
  const cvec & capbuf,
  const double & thresh2_n_sigma,
  const double & fc_requested,
  const double & fc_programmed,
  const double & fs_programmed,
  // Only used for testing...
  vec & sss_h1_np_est,
  vec & sss_h2_np_est,
  cvec & sss_h1_nrm_est,
  cvec & sss_h2_nrm_est,
  cvec & sss_h1_ext_est,
  cvec & sss_h2_ext_est,
  mat & log_lik_nrm,
  mat & log_lik_ext
) {
  // Get the channel estimates and extract the raw SSS subcarriers
  sss_detect_getce_sss(cell,capbuf,fc_requested,fc_programmed,fs_programmed,sss_h1_np_est,sss_h2_np_est,sss_h1_nrm_est,sss_h2_nrm_est,sss_h1_ext_est,sss_h2_ext_est);
  // Perform maximum likelihood detection
  sss_detect_ml(cell,sss_h1_np_est,sss_h2_np_est,sss_h1_nrm_est,sss_h2_nrm_est,sss_h1_ext_est,sss_h2_ext_est,log_lik_nrm,log_lik_ext);

  // Determine normal/ extended CP
  mat log_lik;
  cp_type_t::cp_type_t cp_type;
  if (max(max(log_lik_nrm))>max(max(log_lik_ext))) {
    log_lik=log_lik_nrm;
    cp_type=cp_type_t::NORMAL;
  } else {
    log_lik=log_lik_ext;
    cp_type=cp_type_t::EXTENDED;
  }

  // Locate the 'frame start' defined as the start of the CP of the frame.
  // The first DFT should be located at frame_start + cp_length.
  // It is expected (not guaranteed!) that a DFT performed at this
  // location will have a measured time offset of 2 samples.
  const double k_factor=(fc_requested-cell.freq)/fc_programmed;
  double frame_start=cell.ind+(128+9-960-2)*16/FS_LTE*fs_programmed*k_factor;
  vec ll;
  if (max(log_lik.get_col(0))>max(log_lik.get_col(1))) {
    ll=log_lik.get_col(0);
  } else {
    ll=log_lik.get_col(1);
    frame_start=frame_start+9600*k_factor*16/FS_LTE*fs_programmed*k_factor;
  }
  frame_start=WRAP(frame_start,-0.5,(2*9600.0-0.5)*16/FS_LTE*fs_programmed*k_factor);

  // Estimate n_id_1.
  int32 n_id_1_est;
  double lik_final=max(ll,n_id_1_est);

  // Second threshold check to weed out some weak signals.
  Cell cell_out(cell);
  vec L=concat(cvectorize(log_lik_nrm),cvectorize(log_lik_ext));
  double lik_mean=mean(L);
  double lik_var=variance(L);
  if (lik_final>=lik_mean+pow(lik_var,0.5)*thresh2_n_sigma) {
    cell_out.n_id_1=n_id_1_est;
    cell_out.cp_type=cp_type;
    cell_out.frame_start=frame_start;
  }

  return cell_out;
}

// Perform FOE using only the PSS and SSS.
// The PSS correlation peak gives us the frequency offset within 2.5kHz.
// The PSS/SSS can be used to estimate the frequency offset within a
// much finer resolution.
Cell pss_sss_foe(
  const Cell & cell_in,
  const cvec & capbuf,
  const double & fc_requested,
  const double & fc_programmed,
  const double & fs_programmed
) {
  const double k_factor=(fc_requested-cell_in.freq)/fc_programmed;

  // Determine where we can find both PSS and SSS
  uint16 pss_sss_dist;
  double first_sss_dft_location;
  if (cell_in.cp_type==cp_type_t::NORMAL) {
    pss_sss_dist=itpp::round_i((128+9)*16/FS_LTE*fs_programmed*k_factor);
    first_sss_dft_location=cell_in.frame_start+(960-128-9-128)*16/FS_LTE*fs_programmed*k_factor;
  } else if (cell_in.cp_type==cp_type_t::EXTENDED) {
    pss_sss_dist=round_i((128+32)*k_factor);
    first_sss_dft_location=cell_in.frame_start+(960-128-32-128)*16/FS_LTE*fs_programmed*k_factor;
  } else {
    throw("Error... check code...");
  }
  uint8 sn;
  first_sss_dft_location=WRAP(first_sss_dft_location,-0.5,9600*2-0.5);
  if (first_sss_dft_location-9600*k_factor>-0.5) {
    first_sss_dft_location-=9600*k_factor;
    sn=10;
  } else {
    sn=0;
  }
  vec sss_dft_loc_set=itpp_ext::matlab_range(first_sss_dft_location,9600*16/FS_LTE*fs_programmed*k_factor,(double)(length(capbuf)-127-pss_sss_dist-100));
  uint16 n_sss=length(sss_dft_loc_set);

  // Loop around for each PSS/SSS pair
  sn=(1-(sn/10))*10;
  complex <double> M(0,0);
  cmat h_raw_fo_pss(n_sss,62);
  cmat h_sm(n_sss,62);
  cmat sss_raw_fo(n_sss,62);
  vec pss_np(n_sss);
#ifndef NDEBUG
  h_raw_fo_pss=NAN;
  h_sm=NAN;
  sss_raw_fo=NAN;
  pss_np=NAN;
#endif
  for (uint16 k=0;k<n_sss;k++) {
    sn=(1-(sn/10))*10;
    uint32 sss_dft_location=round_i(sss_dft_loc_set(k));

    // Find the PSS and use it to estimate the channel.
    uint32 pss_dft_location=sss_dft_location+pss_sss_dist;
    h_raw_fo_pss.set_row(k,extract_psss(capbuf.mid(pss_dft_location,128),-cell_in.freq,k_factor,fs_programmed));
    h_raw_fo_pss.set_row(k,elem_mult(h_raw_fo_pss.get_row(k),conj(ROM_TABLES.pss_fd[cell_in.n_id_2])));

    // Smoothing... Basic...
    for (uint8 t=0;t<62;t++) {
      uint8 lt=MAX(0,t-6);
      uint8 rt=MIN(61,t+6);
      h_sm(k,t)=mean(h_raw_fo_pss.get_row(k).mid(lt,rt-lt+1));
    }

    // Estimate noise power.
    pss_np(k)=sigpower(h_sm.get_row(k)-h_raw_fo_pss.get_row(k));

    // Calculate the SSS in the frequency domain
    sss_raw_fo.set_row(k,extract_psss(capbuf.mid(sss_dft_location,128),-cell_in.freq,k_factor,fs_programmed)*exp(J*pi*-cell_in.freq/(FS_LTE/16/2)*-pss_sss_dist));
    sss_raw_fo.set_row(k,elem_mult(sss_raw_fo.get_row(k),to_cvec(ROM_TABLES.sss_fd(cell_in.n_id_1,cell_in.n_id_2,sn))));

    // Compare PSS to SSS. With no frequency offset, arg(M) is zero.
    M=M+sum(elem_mult(
      conj(sss_raw_fo.get_row(k)),
      h_raw_fo_pss.get_row(k),
      to_cvec(elem_mult(
        sqr(h_sm.get_row(k)),
        1.0/(2*sqr(h_sm.get_row(k))*pss_np(k)+sqr(pss_np(k)))
      ))
    ));
  }

  // Store results.
  Cell cell_out(cell_in);
  cell_out.freq_fine=cell_in.freq+arg(M)/(2*pi)/(1/(fs_programmed*k_factor)*pss_sss_dist);
  return cell_out;
}

// Extract the time/ frequency grid.
//
// Note that this function is inefficient in that it returns the time/
// frequency grid for nearly all samples in the capture buffer whereas
// in reality, we are only interested in the OFDM symbols containing the MIB.
void extract_tfg(
  // Inputs
  const Cell & cell,
  const cvec & capbuf_raw,
  const double & fc_requested,
  const double & fc_programmed,
  const double & fs_programmed,
  // Outputs
  cmat & tfg,
  vec & tfg_timestamp
) {
  // Local shortcuts
  const double frame_start=cell.frame_start;
  const cp_type_t::cp_type_t cp_type=cell.cp_type;
  const double freq_fine=cell.freq_fine;

  // Derive some values
  // fc*k_factor is the receiver's actual RX center frequency.
  const double k_factor=(fc_requested-cell.freq_fine)/fc_programmed;
  const int8 n_symb_dl=cell.n_symb_dl();
  double dft_location;
  if (cp_type==cp_type_t::NORMAL) {
    dft_location=frame_start+10*16/FS_LTE*fs_programmed*k_factor;
  } else if (cp_type==cp_type_t::EXTENDED) {
    dft_location=frame_start+32*16/FS_LTE*fs_programmed*k_factor;
  } else {
    throw("Check code...");
  }

  // See if we can advance the frame start
  if (dft_location-.01*fs_programmed*k_factor>-0.5) {
    dft_location=dft_location-.01*fs_programmed*k_factor;
  }

  // Perform FOC
  cvec capbuf=fshift(capbuf_raw,-freq_fine,fs_programmed*k_factor);

  // Extract 6 frames + 2 slots worth of data
  uint16 n_ofdm_sym=6*10*2*n_symb_dl+2*n_symb_dl;
  tfg=cmat(n_ofdm_sym,72);
  tfg_timestamp=vec(n_ofdm_sym);
#ifndef NDEBUG
  tfg=NAN;
  tfg_timestamp=NAN;
#endif
  uint16 sym_num=0;
  for (uint16 t=0;t<n_ofdm_sym;t++) {
    cvec dft_out=dft(capbuf.mid(round_i(dft_location),128));
    tfg.set_row(t,concat(dft_out.right(36),dft_out.mid(1,36)));
    // Record the time offset where the DFT _should_ have been taken.
    // It was actually taken at the nearest sample boundary.
    tfg_timestamp(t)=dft_location;
    // Calculate location of next DFT
    if (n_symb_dl==6) {
      dft_location+=(128+32)*16/FS_LTE*fs_programmed*k_factor;
    } else {
      if (sym_num==6) {
        dft_location+=(128+10)*16/FS_LTE*fs_programmed*k_factor;
      } else {
        dft_location+=(128+9)*16/FS_LTE*fs_programmed*k_factor;
      }
      sym_num=mod(sym_num+1,7);
    }
  }

  // Compensate for the residual time offset.
  ivec cn=concat(itpp_ext::matlab_range(-36,-1),itpp_ext::matlab_range(1,36));
  for (uint16 t=0;t<n_ofdm_sym;t++) {
    double ideal_offset=tfg_timestamp(t);
    double actual_offset=round_i(ideal_offset);
    // How late were we in locating the DFT
    double late=actual_offset-ideal_offset;
    // Compensate for the improper location of the DFT
    tfg.set_row(t,elem_mult(tfg.get_row(t),exp((-J*2*pi*late/128)*cn)));
  }
  // At this point, tfg(t,:) contains the results of a DFT that was performed
  // at time offset tfg_timestamp(t). Note that tfg_timestamp(t) is not an
  // integer!
}

// Perform 'superfine' TOE/FOE/TOC/FOC.
//
// First, the residual frequency offset is measured using all the samples
// in the TFG. This frequency offset will be less noisy than the one produced
// by comparing the SSS and PSS, and will have significantly better performance
// in low SNR situations. In high SNR situations, the PSS/SSS measured
// frequency offset is sufficient.
//
// Then, FOC is performed. At this point, if the UE had transmitted the same
// symbol on all subcarriers of all OFDM symbols:
//   arg(E[conj(tfg(t,f))*tfg(t+1,f)])=0;
//
// Next, TOE is performed followed by TOC. At this point, if the UE had
// transmitted the same symbol on all subcarriers of all OFDM symbols:
//   arg(E[conj(tfg(t,f))*tfg(t,f+1)])=0;
Cell tfoec(
  // Inputs
  const Cell & cell,
  const cmat & tfg,
  const vec & tfg_timestamp,
  const double & fc_requested,
  const double & fc_programmed,
  const RS_DL & rs_dl,
  // Outputs
  cmat & tfg_comp,
  vec & tfg_comp_timestamp
) {
  // Local shortcuts
  const int8 n_symb_dl=cell.n_symb_dl();
  uint16 n_ofdm=tfg.rows();
  uint16 n_slot=floor(((double)n_ofdm)/n_symb_dl);

  // Perform super-fine FOE
  complex <double> foe;
  for (uint8 sym_num=0;sym_num<=n_symb_dl-3;sym_num+=n_symb_dl-3) {
    // Extract all the RS and compensate for the known transmitted symbols.
    cmat rs_extracted(n_slot,12);
#ifndef NDEBUG
    rs_extracted=NAN;
#endif
    for (uint16 t=0;t<n_slot;t++) {
      // Extract RS
      rs_extracted.set_row(t,tfg.get_row(t*n_symb_dl+sym_num).get(itpp_ext::matlab_range((uint32)rs_dl.get_shift(mod(t,20),sym_num,0),(uint32)6,(uint32)71)));
      // Compensate
      rs_extracted.set_row(t,elem_mult(rs_extracted.get_row(t),conj(rs_dl.get_rs(mod(t,20),sym_num))));
    }
    // FOE, subcarrier by subcarrier.
    for (uint16 t=0;t<12;t++) {
      cvec col=rs_extracted.get_col(t);
      foe=foe+sum(elem_mult(conj(col(0,n_slot-2)),col(1,-1)));
    }
  }
  double residual_f=arg(foe)/(2*pi)/0.0005;

  // Perform FOC. Does not fix ICI!
  double k_factor_residual=(fc_requested-residual_f)/fc_programmed;
  tfg_comp=cmat(n_ofdm,72);
#ifndef NDEBUG
  tfg_comp=NAN;
#endif
  tfg_comp_timestamp=k_factor_residual*tfg_timestamp;
  ivec cn=concat(itpp_ext::matlab_range(-36,-1),itpp_ext::matlab_range(1,36));
  for (uint16 t=0;t<n_ofdm;t++) {
    tfg_comp.set_row(t,tfg.get_row(t)*exp(J*2*pi*-residual_f*tfg_comp_timestamp(t)/(FS_LTE/16)));
    // How late were we in locating the DFT
    double late=tfg_timestamp(t)-tfg_comp_timestamp(t);
    // Compensate for the improper location of the DFT
    tfg_comp.set_row(t,elem_mult(tfg_comp.get_row(t),exp((-J*2*pi*late/128)*cn)));
  }

  // Perform TOE.
  // Implemented by comparing subcarrier k of one OFDM symbol with subcarrier
  // k+3 of another OFDM symbol. This is why FOE must be performed first.
  // Slightly less performance but faster execution time could be obtained
  // by comparing subcarrier k with subcarrier k+6 of the same OFDM symbol.
  complex <double> toe=0;
  for (uint16 t=0;t<2*n_slot-1;t++) {
    // Current OFDM symbol containing RS
    uint8 current_sym_num=(t&1)?(n_symb_dl-3):0;
    uint8 current_slot_num=mod((t>>1),20);
    uint16 current_offset=(t>>1)*n_symb_dl+current_sym_num;
    // Since we are using port 0, the shift is the same for all slots.
    uint8 current_shift=rs_dl.get_shift(0,current_sym_num,0);
    // Next OFDM symbol containing RS
    uint8 next_sym_num=((t+1)&1)?(n_symb_dl-3):0;
    uint8 next_slot_num=mod(((t+1)>>1),20);
    uint16 next_offset=((t+1)>>1)*n_symb_dl+next_sym_num;
    // Since we are using port 0, the shift is the same for all slots.
    uint8 next_shift=rs_dl.get_shift(0,next_sym_num,0);

    uint16 r1_offset,r2_offset;
    uint8 r1_shift,r2_shift;
    uint8 r1_sym_num,r2_sym_num;
    uint8 r1_slot_num,r2_slot_num;
    if (current_shift<next_shift) {
      r1_offset=current_offset;
      r1_shift=current_shift;
      r1_sym_num=current_sym_num;
      r1_slot_num=current_slot_num;
      r2_offset=next_offset;
      r2_shift=next_shift;
      r2_sym_num=next_sym_num;
      r2_slot_num=next_slot_num;
    } else {
      r1_offset=next_offset;
      r1_shift=next_shift;
      r1_sym_num=next_sym_num;
      r1_slot_num=next_slot_num;
      r2_offset=current_offset;
      r2_shift=current_shift;
      r2_sym_num=current_sym_num;
      r2_slot_num=current_slot_num;
    }
    cvec r1v=tfg_comp.get_row(r1_offset).get(itpp_ext::matlab_range(r1_shift,6,71));
    r1v=elem_mult(r1v,conj(rs_dl.get_rs(r1_slot_num,r1_sym_num)));
    cvec r2v=tfg_comp.get_row(r2_offset).get(itpp_ext::matlab_range(r2_shift,6,71));
    r2v=elem_mult(r2v,conj(rs_dl.get_rs(r2_slot_num,r2_sym_num)));
    complex<double> toe1=sum(elem_mult(conj(r1v),r2v));
    complex<double> toe2=sum(elem_mult(conj(r2v(0,10)),r1v(1,11)));
    toe+=toe1+toe2;
  }
  double delay=-arg(toe)/3/(2*pi/128);

  // Perform TOC
  cvec comp_vector=exp((J*2*pi/128*delay)*cn);
  for (uint16 t=0;t<n_ofdm;t++) {
    tfg_comp.set_row(t,elem_mult(tfg_comp.get_row(t),comp_vector));
  }

  Cell cell_out(cell);
  cell_out.freq_superfine=cell_out.freq_fine+residual_f;
  return cell_out;
}

// Helper function to remove entries that are out of bounds
void del_oob(
  ivec & v
) {
  int32 t=0;
  while (t<v.length()) {
    if ((v(t)<0)||(v(t)>11)) {
      v.del(t);
    } else {
      t++;
    }
  }
}

// Interpolate between the filtered channel estimates and the full
// time/frequency grid.
// This algorithm interpolates in the frequency domain and then
// in the time domain.
void ce_interp_freq_time(
  // Inputs
  const cmat & ce_filt,
  const ivec & shift,
  const int16 & n_ofdm,
  const int16 & n_rs_ofdm,
  const ivec & rs_set,
  // Outputs
  cmat & ce_tfg
) {
  // Interpolate in the frequency domain
  cmat ce_filt_frq(n_rs_ofdm,72);
#ifndef NDEBUG
  ce_filt_frq=NAN;
#endif
  for (int32 t=0;t<n_rs_ofdm;t++) {
    vec X=itpp_ext::matlab_range((double)shift(t&1),6.0,71.0);
    cvec Y=ce_filt.get_row(t);
    vec x=itpp_ext::matlab_range(0.0,71.0);
    ce_filt_frq.set_row(t,interp1(X,Y,x));
  }

  // Interpolate in the time domain
  ce_tfg=cmat(n_ofdm,72);
  for (uint8 t=0;t<=71;t++) {
    vec X=to_vec(rs_set);
    cvec Y=ce_filt_frq.get_col(t);
    vec x=itpp_ext::matlab_range(0.0,n_ofdm-1.0);
    ce_tfg.set_col(t,interp1(X,Y,x));
  }
}

// Interpolate between the filtered channel estimates and the full
// time/frequency grid.
// This algorithm first creates a uniformly spaced grid from the hexagonally
// spaced RS grid and then interpolates between the uniformly spaced grid.
void ce_interp_2stage(
  // Inputs
  const cmat & ce_filt,
  const ivec & shift,
  const int16 & n_ofdm,
  const int16 & n_rs_ofdm,
  const ivec & rs_set,
  // Outputs
  cmat & ce_tfg
) {
  // Fill in some missing samples to create a uniformly sampled grid.
  cmat ce_filt_exp(n_rs_ofdm,24);
#ifndef NDEBUG
  ce_filt_exp=NAN;
#endif
  bool current_row_leftmost=shift(0)<shift(1);
  for (uint16 t=0;t<n_rs_ofdm;t++) {
    for (int16 k=0;k<24;k++) {
      if ((k&1)==current_row_leftmost) {
        complex <double> total=0;
        uint8 n_total=0;
        // Up one
        if (t-1>=0) {
          total+=ce_filt(t-1,k>>1);
          n_total++;
        }
        // Down one
        if (t+1<n_rs_ofdm) {
          total+=ce_filt(t+1,k>>1);
          n_total++;
        }
        // Left
        if (((k-1)>>1)>=0) {
          total+=ce_filt(t,(k-1)>>1);
          n_total++;
        }
        // Right
        if (((k+1)>>1)<12) {
          total+=ce_filt(t,(k+1)>>1);
          n_total++;
        }
        ce_filt_exp(t,k)=total/n_total;
      } else {
        // There is already a sample here
        ce_filt_exp(t,k)=ce_filt(t,k>>1);
      }
    }
    current_row_leftmost=!current_row_leftmost;
  }
  ivec ce_filt_exp_x=itpp_ext::matlab_range(min(shift),3,71);

  // Interpolate (linearly) the uniformly sampled grid to create channel
  // estimates for every RE.
  ce_tfg=cmat(n_ofdm,72);
#ifndef NDEBUG
  ce_tfg=NAN;
#endif
  // First, interpolate in the frequency dimension.
  for (uint16 t=0;t<n_rs_ofdm;t++) {
    ivec X=ce_filt_exp_x;
    cvec Y=ce_filt_exp.get_row(t);
    ivec x=itpp_ext::matlab_range(0,71);
    ce_tfg.set_row(rs_set(t),interp1(to_vec(X),Y,to_vec(x)));
  }
  // Now, interpolate in the time dimension.
  for (uint16 t=0;t<72;t++) {
    ivec X=rs_set;
    cvec Y=ce_tfg.get_col(t).get(rs_set);
    ivec x=itpp_ext::matlab_range(0,n_ofdm-1);
    ce_tfg.set_col(t,interp1(to_vec(X),Y,to_vec(x)));
  }
}

// Helper function. Linearly interpolates the left/right edge samples so as to
// guarantee a vertex at the left and right subcarrier.
void ce_interp_hex_extend(
  vec & row_x,
  cvec & row_val
) {
  if (row_x(0)!=0) {
    row_val.ins(0,row_val(0)-row_x(0)*(row_val(1)-row_val(0))/(row_x(1)-row_x(0)));
    row_x.ins(0,0);
  }
  if (itpp_ext::last(row_x)!=71) {
    uint16 len=length(row_val);
    row_val.ins(len,row_val(len-1)+(71-itpp_ext::last(row_x))*(row_val(len-1)-row_val(len-2))/(row_x(len-1)-row_x(len-2)));
    row_x.ins(len,71);
  }
}

// Use Delaunay triangles to perform interpolation. (Similar to Matlab's
// griddata function.)
// This struct is only used by this function.
struct triangle_vertex_t {
  uint8 x_sc;
  uint16 y_symnum;
  complex <double> val;
};
void ce_interp_hex(
  // Inputs
  const cmat & ce_filt,
  const ivec & shift,
  const int16 & n_ofdm,
  const int16 & n_rs_ofdm,
  const ivec & rs_set,
  // Outputs
  cmat & ce_tfg
) {
  ce_tfg=cmat(n_ofdm,72);
#ifndef NDEBUG
  ce_tfg=NAN;
#endif
  for (uint16 t=0;t<=n_rs_ofdm-2;t++) {
    // Extract two rows and ensure that there are samples for subcarriers
    // 0 and 71.
    // In general, top_row_* is actually equal to bot_row_* from the
    // previous iteration.
    vec top_row_x=to_vec(itpp_ext::matlab_range((t&1)?shift(1):shift(0),6,71));
    cvec top_row_val=ce_filt.get_row(t);
    ce_interp_hex_extend(top_row_x,top_row_val);
    vec bot_row_x=to_vec(itpp_ext::matlab_range((t&1)?shift(0):shift(1),6,71));
    cvec bot_row_val=ce_filt.get_row(t+1);
    ce_interp_hex_extend(bot_row_x,bot_row_val);

    // First row is not handled inside the main loop.
    if (t==0) {
      ce_tfg.set_row(rs_set(0),interp1(top_row_x,top_row_val,itpp_ext::matlab_range(0.0,71.0)));
    }

    // Create initial triangle
    uint8 top_row_last_used;
    uint8 bot_row_last_used;
    vector <triangle_vertex_t> triangle(3);
    if (top_row_x(1)<bot_row_x(1)) {
      triangle[0].x_sc=top_row_x(0);
      triangle[0].y_symnum=rs_set(t);
      triangle[0].val=top_row_val(0);
      triangle[1].x_sc=bot_row_x(0);
      triangle[1].y_symnum=rs_set(t+1);
      triangle[1].val=bot_row_val(0);
      triangle[2].x_sc=top_row_x(1);
      triangle[2].y_symnum=rs_set(t);
      triangle[2].val=top_row_val(1);
      top_row_last_used=1;
      bot_row_last_used=0;
    } else {
      triangle[0].x_sc=bot_row_x(0);
      triangle[0].y_symnum=rs_set(t+1);
      triangle[0].val=bot_row_val(0);
      triangle[1].x_sc=top_row_x(0);
      triangle[1].y_symnum=rs_set(t);
      triangle[1].val=top_row_val(0);
      triangle[2].x_sc=bot_row_x(1);
      triangle[2].y_symnum=rs_set(t+1);
      triangle[2].val=bot_row_val(1);
      top_row_last_used=0;
      bot_row_last_used=1;
    }

    // This loop succesively creates triangles to cover the space between
    // top_row and bot_row.
    uint8 spacing=rs_set(t+1)-rs_set(t);
    vec x_offset(spacing+1);
    x_offset=0.0;
    while (true) {
      // Calculate the parameters of a plane passing through all points
      // of the triangle.
      // value=a_p*x_sc+b_p*y_symnum+c_p;
      cmat M(3,3);
      M(0,0)=triangle[0].x_sc;
      M(1,0)=triangle[1].x_sc;
      M(2,0)=triangle[2].x_sc;
      M(0,1)=triangle[0].y_symnum;
      M(1,1)=triangle[1].y_symnum;
      M(2,1)=triangle[2].y_symnum;
      M(0,2)=1;
      M(1,2)=1;
      M(2,2)=1;
      cvec V(3);
      V(0)=triangle[0].val;
      V(1)=triangle[1].val;
      V(2)=triangle[2].val;
      // In the future, inv(M) can be calculated directly for speed
      // since the last column of M is all ones.
      cvec abc=inv(M)*V;
      complex <double> a_p=abc(0);
      complex <double> b_p=abc(1);
      complex <double> c_p=abc(2);

      // Calculate the parameters of the line defining the rightmost
      // edge of the triangle.
      // x_sc=a_l*y_symnum+b_l;
      double x1=triangle[1].x_sc;
      double x2=triangle[2].x_sc;
      double y1=triangle[1].y_symnum;
      double y2=triangle[2].y_symnum;
      double a_l=(x1-x2)/(y1-y2);
      double b_l=(y1*x2-y2*x1)/(y1-y2);

      for (uint8 r=1;r<=spacing;r++) {
        while (x_offset(r)<=a_l*(rs_set(t)+r)+b_l) {
          ce_tfg(rs_set(t)+r,x_offset(r))=a_p*x_offset(r)+b_p*(rs_set(t)+r)+c_p;
          x_offset(r)++;
        }
      }

      if ((x_offset(1)==72)&&(itpp_ext::last(x_offset)==72)) {
        break;
      }

      // We are not done yet. Choose the points for the next triangle.
      if (triangle[2].y_symnum==rs_set(t)) {
        triangle[0]=triangle[1];
        triangle[1]=triangle[2];
        bot_row_last_used++;
        triangle[2].x_sc=bot_row_x(bot_row_last_used);
        triangle[2].y_symnum=rs_set(t+1);
        triangle[2].val=bot_row_val(bot_row_last_used);
      } else {
        triangle[0]=triangle[1];
        triangle[1]=triangle[2];
        top_row_last_used++;
        triangle[2].x_sc=top_row_x(top_row_last_used);
        triangle[2].y_symnum=rs_set(t);
        triangle[2].val=top_row_val(top_row_last_used);
      }
    }
  }

  // Rows before the first and after the last OFDM symbol with RS are
  // created simply by copying the nearest OFDM symbol with RS.
  for (uint8 t=0;t<rs_set(0);t++) {
    ce_tfg.set_row(t,ce_tfg.get_row(rs_set(0)));
  }
  for (uint16 t=itpp_ext::last(rs_set)+1;t<n_ofdm;t++) {
    ce_tfg.set_row(t,ce_tfg.get_row(itpp_ext::last(rs_set)));
  }
}

// Examine the time/frequency grid and perform channel estimation
// and filtering to produce channel estimates for every RE.
//
// Each invocation of this function performs CE/ filtering for one particular
// antenna port.
void chan_est(
  // Inputs
  const Cell & cell,
  const RS_DL & rs_dl,
  const cmat & tfg,
  const uint8 & port,
  // Outputs
  cmat & ce_tfg,
  double & np
) {
  const int8 n_symb_dl=cell.n_symb_dl();
  const uint16 n_ofdm=tfg.rows();

  // Set of OFDM symbols containing reference symbols.
  ivec rs_set;
  if (port<=1) {
    // There are better ways to implement this...
    Sort <int> sort;
    rs_set=concat(itpp_ext::matlab_range(0,n_symb_dl,n_ofdm-1),itpp_ext::matlab_range(n_symb_dl-3,n_symb_dl,n_ofdm-1));
    sort.sort(0,rs_set.length()-1,rs_set);
    //rs_set=reverse(rs_set);
  } else {
    rs_set=itpp_ext::matlab_range(1,n_symb_dl,n_ofdm-1);
  }
  const uint16 n_rs_ofdm=length(rs_set);

  // Extract the raw channel estimates. 12 raw channel estimates per OFDM
  // symbol containing RS.
  cmat ce_raw(n_rs_ofdm,12);
#ifndef NDEBUG
  ce_raw=NAN;
#endif
  uint8 slot_num=0;
  ivec shift(2);
  shift=-1000;
  for (uint16 t=0;t<n_rs_ofdm;t++) {
    uint8 sym_num=mod(rs_set(t),n_symb_dl);
    if (t<=1) {
      shift(t)=rs_dl.get_shift(mod(slot_num,20),sym_num,port);
    }

    cvec rs=rs_dl.get_rs(slot_num,sym_num);
    // Extract
    cvec raw_row=tfg.get_row(rs_set(t)).get(itpp_ext::matlab_range((int32)rs_dl.get_shift(mod(slot_num,20),sym_num,port),6,71));
    ce_raw.set_row(t,raw_row);
    // Compensate for known RS
    ce_raw.set_row(t,elem_mult(ce_raw.get_row(t),conj(rs)));
    if (((t&1)==1)||(port>=2)) {
      slot_num=mod(slot_num+1,20);
    }
  }

  // Simple filtering.
  //
  // 1   2   3   4   5   6
  //   7   8   9   A   B
  // C   D   E   F   G   H
  //
  // If the above is a representation of the location of the raw channel
  // estimates in the TFG, the filtered channel estimate for position 8
  // is the mean of the raw channel estimates for positions 2,3,7,8,9,D,
  // and E.
  cmat ce_filt(n_rs_ofdm,12);
  bool current_row_leftmost=shift(0)<shift(1);
  for (uint16 t=0;t<n_rs_ofdm;t++) {
    complex <double> total;
    uint8 n_total;
    for (uint16 k=0;k<12;k++) {
      // Current time offset
      ivec ind;
      ind=itpp_ext::matlab_range(k-1,k+1);
      del_oob(ind);
      total=sum(ce_raw.get_row(t).get(ind));
      n_total=length(ind);

      // Add in the previous and next time offset (if they exist)
      if (shift(0)==shift(1)) {
        ind=itpp_ext::matlab_range(k-1,k+1);
      } else {
        if (current_row_leftmost) {
          ind=itpp_ext::matlab_range(k-1,k);
        } else {
          ind=itpp_ext::matlab_range(k,k+1);
        }
      }
      del_oob(ind);
      // Previous time offset
      if (t!=0) {
        total+=sum(ce_raw.get_row(t-1).get(ind));
        n_total+=length(ind);
      }
      if (t!=n_rs_ofdm-1) {
        total+=sum(ce_raw.get_row(t+1).get(ind));
        n_total+=length(ind);
      }
      ce_filt(t,k)=total/n_total;
    }
    current_row_leftmost=!current_row_leftmost;
  }

  // Estimate noise power
  np=sigpower(cvectorize(ce_filt)-cvectorize(ce_raw));

  // There is no appreciable difference in performance between these
  // algorithms for high SNR values.
  //ce_interp_2stage(ce_filt,shift,n_ofdm,n_rs_ofdm,rs_set,ce_tfg);
  //ce_interp_freq_time(ce_filt,shift,n_ofdm,n_rs_ofdm,rs_set,ce_tfg);
  ce_interp_hex(ce_filt,shift,n_ofdm,n_rs_ofdm,rs_set,ce_tfg);
}

// Examine the time/ frequency grid and extract the RE that belong to the PBCH.
// Also return the channel estimates for that RE from all 4 possible eNodeB
// ports.
void pbch_extract(
  // Inputs
  const Cell & cell,
  const cmat & tfg,
  const Array <cmat> & ce,
  // Outputs
  cvec & pbch_sym,
  cmat & pbch_ce
) {
  // Shortcuts
  const int8 n_symb_dl=cell.n_symb_dl();
  const uint16 m_bit=(cell.cp_type==cp_type_t::NORMAL)?1920:1728;
  const uint8 v_shift_m3=mod(cell.n_id_cell(),3);

  pbch_sym=cvec(m_bit/2);
  // One channel estimate from each of 4 ports for each RE.
  pbch_ce=cmat(4,m_bit/2);
#ifndef NDEBUG
  pbch_sym=NAN;
  pbch_ce=NAN;
#endif
  uint32 idx=0;
  for (uint8 fr=0;fr<=3;fr++) {
    for (uint8 sym=0;sym<=3;sym++) {
      for (uint8 sc=0;sc<=71;sc++) {
        // Skip if there might be an RS occupying this position.
        if ((mod(sc,3)==v_shift_m3)&&((sym==0)||(sym==1)||((sym==3)&&(n_symb_dl==6)))) {
          continue;
        }
        uint16 sym_num=fr*10*2*n_symb_dl+n_symb_dl+sym;
        pbch_sym(idx)=tfg(sym_num,sc);
        pbch_ce(0,idx)=ce(0).get(sym_num,sc);
        pbch_ce(1,idx)=ce(1).get(sym_num,sc);
        pbch_ce(2,idx)=ce(2).get(sym_num,sc);
        pbch_ce(3,idx)=ce(3).get(sym_num,sc);
        idx++;
      }
    }
  }
  ASSERT(idx==m_bit/2);
}

// Blindly try various frame alignments and numbers of antennas to try
// to find a valid MIB.
Cell decode_mib(
  const Cell & cell,
  const cmat & tfg,
  const RS_DL & rs_dl
) {
  // Local shortcuts
  const int8 n_symb_dl=cell.n_symb_dl();

  Cell cell_out=cell;

  // Channel estimation. This is automatically performed for four antennas
  // and for every RE, not only the RE's that contain an MIB!!!
  Array <cmat> ce_tfg(4);
  vec np_v(4);
  chan_est(cell,rs_dl,tfg,0,ce_tfg(0),np_v(0));
  chan_est(cell,rs_dl,tfg,1,ce_tfg(1),np_v(1));
  chan_est(cell,rs_dl,tfg,2,ce_tfg(2),np_v(2));
  chan_est(cell,rs_dl,tfg,3,ce_tfg(3),np_v(3));

  // Try various frame offsets and number of TX antennas.
  bvec c_est;
  for (uint8 frame_timing_guess=0;frame_timing_guess<=3;frame_timing_guess++) {
    const uint16 ofdm_sym_set_start=frame_timing_guess*10*2*n_symb_dl;
    ivec ofdm_sym_set=itpp_ext::matlab_range(ofdm_sym_set_start,ofdm_sym_set_start+3*10*2*n_symb_dl+2*n_symb_dl-1);

    // Extract only the portion of the TFG containing the four frames
    // we are interested in.
    cmat tfg_try=tfg.get_rows(ofdm_sym_set);
    Array <cmat> ce_try(4);
    for (uint8 t=0;t<4;t++) {
      ce_try(t)=ce_tfg(t).get_rows(ofdm_sym_set);
    }

    // Extract symbols and channel estimates for the PBCH
    cvec pbch_sym;
    cmat pbch_ce;
    pbch_extract(cell,tfg_try,ce_try,pbch_sym,pbch_ce);

    // Try 1, 2, and 4 ports.
    vec np;
    cvec syms;
    for (uint8 n_ports_pre=1;n_ports_pre<=3;n_ports_pre++) {
      const uint8 n_ports=(n_ports_pre==3)?4:n_ports_pre;
      // Perform channel compensation and also estimate noise power in each
      // symbol.
      if (n_ports==1) {
        cvec gain=conj(elem_div(pbch_ce.get_row(0),to_cvec(sqr(pbch_ce.get_row(0)))));
        syms=elem_mult(pbch_sym,gain);
        np=np_v(0)*sqr(gain);
      } else {
        syms.set_size(length(pbch_sym));
        np.set_size(length(pbch_sym));
#ifndef NDEBUG
        syms=NAN;
        np=NAN;
#endif
        for (int32 t=0;t<length(syms);t+=2) {
          // Simple zero-forcing
          // http://en.wikipedia.org/wiki/Space-time_block_coding_based_transmit_diversity
          complex <double> h1,h2;
          double np_temp;
          if (n_ports==2) {
            h1=(pbch_ce(0,t)+pbch_ce(0,t+1))/2;
            h2=(pbch_ce(1,t)+pbch_ce(1,t+1))/2;
            np_temp=mean(np_v(0,1));
          } else {
            if (mod(t,4)==0) {
              h1=(pbch_ce(0,t)+pbch_ce(0,t+1))/2;
              h2=(pbch_ce(2,t)+pbch_ce(2,t+1))/2;
              np_temp=(np_v(0)+np_v(2))/2;
            } else {
              h1=(pbch_ce(1,t)+pbch_ce(1,t+1))/2;
              h2=(pbch_ce(3,t)+pbch_ce(3,t+1))/2;
              np_temp=(np_v(1)+np_v(3))/2;
            }
          }
          complex <double> x1=pbch_sym(t);
          complex <double> x2=pbch_sym(t+1);
          double scale=pow(h1.real(),2)+pow(h1.imag(),2)+pow(h2.real(),2)+pow(h2.imag(),2);
          syms(t)=(conj(h1)*x1+h2*conj(x2))/scale;
          syms(t+1)=conj((-conj(h2)*x1+h1*conj(x2))/scale);
          np(t)=(pow(abs(h1)/scale,2)+pow(abs(h2)/scale,2))*np_temp;
          np(t+1)=np(t);
        }
        // 3dB factor comes from precoding for transmit diversity
        syms=syms*pow(2,0.5);
      }

      // Extract the bits from the complex modulated symbols.
      vec e_est=lte_demodulate(syms,np,modulation_t::QAM);
      // Unscramble
      bvec scr=lte_pn(cell.n_id_cell(),length(e_est));
      for (int32 t=0;t<length(e_est);t++) {
        if (scr(t)) e_est(t)=-e_est(t);
      }
      // Undo ratematching
      mat d_est=lte_conv_deratematch(e_est,40);
      // Decode
      c_est=lte_conv_decode(d_est);
      // Calculate received CRC
      bvec crc_est=lte_calc_crc(c_est(0,23),CRC16);
      // Apply CRC mask
      if (n_ports==2) {
        for (uint8 t=0;t<16;t++) {
          crc_est(t)=1-((int)crc_est(t));
        }
      } else if (n_ports==4) {
        for (uint8 t=1;t<length(crc_est);t+=2) {
          crc_est(t)=1-((int)crc_est(t));
        }
      }
      // Did we find it?
      if (crc_est==c_est(24,-1)) {
        // YES!
        cell_out.n_ports=n_ports;
        // Unpack the MIB
        ivec c_est_ivec=to_ivec(c_est);
        // DL bandwidth
        const uint8 bw_packed=c_est_ivec(0)*4+c_est_ivec(1)*2+c_est_ivec(2);
        switch (bw_packed) {
          case 0:
            cell_out.n_rb_dl=6;
            break;
          case 1:
            cell_out.n_rb_dl=15;
            break;
          case 2:
            cell_out.n_rb_dl=25;
            break;
          case 3:
            cell_out.n_rb_dl=50;
            break;
          case 4:
            cell_out.n_rb_dl=75;
            break;
          case 5:
            cell_out.n_rb_dl=100;
            break;
        }
        // PHICH duration
        cell_out.phich_duration=c_est_ivec(3)?phich_duration_t::EXTENDED:phich_duration_t::NORMAL;
        // PHICH resources
        uint8 phich_res=c_est_ivec(4)*2+c_est_ivec(5);
        switch (phich_res) {
          case 0:
            cell_out.phich_resource=phich_resource_t::oneSixth;
            break;
          case 1:
            cell_out.phich_resource=phich_resource_t::half;
            break;
          case 2:
            cell_out.phich_resource=phich_resource_t::one;
            break;
          case 3:
            cell_out.phich_resource=phich_resource_t::two;
            break;
        }
        // Calculate SFN
        int8 sfn_temp=128*c_est_ivec(6)+64*c_est_ivec(7)+32*c_est_ivec(8)+16*c_est_ivec(9)+8*c_est_ivec(10)+4*c_est_ivec(11)+2*c_est_ivec(12)+c_est_ivec(13);
        cell_out.sfn=itpp_ext::matlab_mod(sfn_temp*4-frame_timing_guess,1024);
        return cell_out;
      }
    }
  }

  return cell_out;
}

