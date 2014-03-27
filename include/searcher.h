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

// Improved by Jiao Xianjun (putaoshu@gmail.com):
// 1. TD-LTE support
// 2. fast pre-search frequencies (external mixer/LNB support)
// 3. multiple tries at one frequency
// 4. .bin file recording and replaying

#ifndef HAVE_SEARCHER_H
#define HAVE_SEARCHER_H

//#define USE_OPENCL // just for debug purpose. It should be removed before formal release

#ifdef USE_OPENCL

#include <CL/cl.h>
#define MAX_NUM_PLATFORM 8
#define MAX_NUM_DEVICE 8
// Class related to OpenCL
class lte_opencl_t {
  public:
    // Initializer
    lte_opencl_t(
      const uint & platform_id,
      const uint & device_id
    );

    // de-Initializer
    virtual ~lte_opencl_t();

    uint platform_id;
    uint device_id;

    cl_uint num_platform;
    cl_platform_id platforms[MAX_NUM_PLATFORM];
    cl_uint num_device;
    cl_device_id devices[MAX_NUM_DEVICE];
    cl_context context;
    cl_command_queue cmdQueue;

    // setup OpenCL environment
    int setup_opencl();

    // for filter_my
    size_t filter_my_length;
    uint filter_my_workitem;
    size_t filter_my_capbuf_length;

    float *filter_my_in_host;
    float *filter_my_out_host;

    cl_mem filter_my_orig;
    cl_mem filter_my_in;
    cl_mem filter_my_mid;
    cl_mem filter_my_out;

    uint filter_my_buf_in_len;
    uint filter_my_buf_mid_len;
    uint filter_my_buf_out_len;

    cl_kernel filter_my_skip2cols;
    cl_kernel filter_my_multi_filter;
    cl_kernel filter_my_result_combine;

    int setup_filter_my(std::string filter_my_kernels_filename, const size_t & capbuf_length_in, const uint & filter_workitem_in);

    int filter_my(itpp::cvec & capbuf);

    // for xcorr_pss
    size_t filter_mchn_length;
    uint filter_mchn_workitem;
    size_t filter_mchn_capbuf_length;
    size_t filter_mchn_num_chn;

    float *filter_mchn_coef_host;
    float *filter_mchn_in_host;
    float *filter_mchn_out_abs2_host;

    cl_mem filter_mchn_coef;
    cl_mem filter_mchn_orig;
    cl_mem filter_mchn_in;
    cl_mem filter_mchn_mid;
    cl_mem filter_mchn_out;
    cl_mem filter_mchn_out_abs2;

    uint filter_mchn_buf_coef_len;
    uint filter_mchn_buf_in_len;
    uint filter_mchn_buf_mid_len;
    uint filter_mchn_buf_out_len;

    cl_kernel filter_mchn_skip2cols;
    cl_kernel filter_mchn_multi_filter;
    cl_kernel filter_mchn_result_combine;

    int setup_filter_mchn(std::string filter_mchn_kernels_filename, const size_t & capbuf_length_in, const size_t & num_filter_in, const size_t & filter_length_in, const uint & xcorr_workitem_in);

    int filter_mchn(const itpp::cvec & capbuf, const itpp::cmat & pss_fo_set, itpp::mat & corr_store);

  private:

};

#else
class lte_opencl_t {
  public:
    // Initializer
    lte_opencl_t(
      const uint & platform_id,
      const uint & device_id);

    uint platform_id;
    uint device_id;
};
#endif

// FIR 6RB filter
void filter_my(
  //Inputs
  const itpp::vec & coef,
  //Outputs
  itpp::cvec & capbuf
);

void pss_fo_set_gen_non_twist(
  // Input
  const itpp::vec & fo_search_set,
  const double & fs_programmed,
  const double & k_factor,
  // Output
  itpp::cmat & pss_fo_set
);

void pss_fo_set_gen_twist(
  // Input
  const itpp::vec & fo_search_set,
  const double & fc_requested,
  const double & fc_programmed,
  const double & fs_programmed,
  // Output
  itpp::cmat & pss_fo_set
);

void pss_fo_set_gen(
  // Input
  const itpp::vec & fo_search_set,
  // Output
  itpp::cmat & pss_fo_set
);

void sampling_ppm_f_search_set_by_pss(
  // Inputs
  lte_opencl_t & lte_ocl,
  const uint16 & num_loop,
  const itpp::cvec & s,
  const itpp::cmat & pss_fo_set,
  const bool & sampling_carrier_twist,
  const uint16 & max_reserve,
  // Inputs&Outputs
  itpp::vec & fo_search_set,
  // Outpus
  itpp::vec & ppm,
  std::vector <itpp::mat> & xc
);

void sampling_ppm_f_search_set_by_pss_old(
  // Inputs
  const itpp::cvec & capbuf,
  const itpp::cmat & pss_fo_set,
  const bool & sampling_carrier_twist,
  // Inputs&Outputs
  itpp::vec & f_search_set,
  // Outpus
  double & ppm,
  vf3d & xc
);

// Correlate the captured data against the PSS.
void xcorr_pss(
  // Inputs
  const itpp::cvec & capbuf,
  const itpp::vec & f_search_set,
  const uint8 & ds_comb_arm,
  const double & fc_requested,
  const double & fc_programmed,
  const double & fs_programmed,
  const std::vector <itpp::mat> & xc,
  // Outputs
  itpp::mat & xc_incoherent_collapsed_pow,
  itpp::imat & xc_incoherent_collapsed_frq,
  // Only used for testing
  std::vector <itpp::mat> & xc_incoherent_single,
  std::vector <itpp::mat> & xc_incoherent,
  itpp::vec & sp_incoherent,
  itpp::vec & sp,
  uint16 & n_comb_xc,
  uint16 & n_comb_sp,
  const bool & sampling_carrier_twist,
  const double k_facotr
);

// Search the correlations for peaks.
void peak_search(
  // Inputs
  const itpp::mat & xc_incoherent_collapsed_pow,
  const itpp::imat & xc_incoherent_collapsed_frq,
  const itpp::vec & Z_th1,
  const itpp::vec & f_search_set,
  const double & fc_requested,
  const double & fc_programmed,
  const std::vector <itpp::mat> & xc_incoherent_single,
  const uint8 & ds_comb_arm,
  const bool & sampling_carrier_twist,
  const double k_factor,
  // Outputs
  std::list <Cell> & cells
);

// For a certain detected PSS, attempt to find the SSS.
Cell sss_detect(
  // Inputs
  const Cell & cell,
  const itpp::cvec & capbuf,
  const double & thresh2_n_sigma,
  const double & fc_requested,
  const double & fc_programmed,
  const double & fs_programmed,
  // Only used for testing
  itpp::vec & sss_h1_np_est,
  itpp::vec & sss_h2_np_est,
  itpp::cvec & sss_h1_nrm_est,
  itpp::cvec & sss_h2_nrm_est,
  itpp::cvec & sss_h1_ext_est,
  itpp::cvec & sss_h2_ext_est,
  itpp::mat & log_lik_nrm,
  itpp::mat & log_lik_ext,
  const bool & sampling_carrier_twist,
  const int & tdd_flag
);

// Perform FOE based only on the PSS and SSS
Cell pss_sss_foe(
  const Cell & cell_in,
  const itpp::cvec & capbuf,
  const double & fc_requested,
  const double & fc_programmed,
  const double & fs_programmed,
  const bool & sampling_carrier_twist,
  const int & tdd_flag
);

// Extract the time and frequency grid.
void extract_tfg(
  // Inputs
  const Cell & cell,
  const itpp::cvec & capbuf_raw,
  const double & fc_requested,
  const double & fc_programmed,
  const double & fs_programmed,
  // Outputs
  itpp::cmat & tfg,
  itpp::vec & tfg_timestamp,
  const bool & sampling_carrier_twist
);

// Perform TOE/FOE/TOC/FOC on the time/ frequency grid.
Cell tfoec(
  // Inputs
  const Cell & cell,
  const itpp::cmat & tfg,
  const itpp::vec & tfg_timestamp,
  const double & fc_requested,
  const double & fc_programmed,
  const RS_DL & rs_dl,
  // Outputs
  itpp::cmat & tfg_comp,
  itpp::vec & tfg_comp_timestamp,
  const bool & sampling_carrier_twist
);

// Attempt to decode the MIB.
Cell decode_mib(
  const Cell & cell,
  const itpp::cmat & tfg,
  const RS_DL & rs_dl
);

// Small helper function that is used by LTE-Tracker
void del_oob(
  itpp::ivec & v
);

#endif

