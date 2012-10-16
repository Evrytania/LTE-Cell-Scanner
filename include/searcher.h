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

#ifndef HAVE_SEARCHER_H
#define HAVE_SEARCHER_H

// Correlate the captured data against the PSS.
void xcorr_pss(
  // Inputs
  const itpp::cvec & capbuf,
  const itpp::vec & f_search_set,
  const uint8 & ds_comb_arm,
  const double & fc_requested,
  const double & fc_programmed,
  const double & fs_programmed,
  // Outputs
  itpp::mat & xc_incoherent_collapsed_pow,
  itpp::imat & xc_incoherent_collapsed_frq,
  // Only used for testing
  vf3d & xc_incoherent_single,
  vf3d & xc_incoherent,
  itpp::vec & sp_incoherent,
  vcf3d & xc,
  itpp::vec & sp,
  uint16 & n_comb_xc,
  uint16 & n_comb_sp
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
  const vf3d & xc_incoherent_single,
  const uint8 & ds_comb_arm,
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
  itpp::mat & log_lik_ext
);

// Perform FOE based only on the PSS and SSS
Cell pss_sss_foe(
  const Cell & cell_in,
  const itpp::cvec & capbuf,
  const double & fc_requested,
  const double & fc_programmed,
  const double & fs_programmed
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
  itpp::vec & tfg_timestamp
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
  itpp::vec & tfg_comp_timestamp
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

