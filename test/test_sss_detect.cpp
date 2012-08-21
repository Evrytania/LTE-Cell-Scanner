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

// This program actually tests both sss_detect() and pss_sss_foe()
#include <itpp/itbase.h>
#include <list>
#include "rtl-sdr.h"
#include "common.h"
#include "macros.h"
#include "lte_lib.h"
#include "searcher.h"
#include "itpp_ext.h"

using namespace std;
using namespace itpp;

#define READ_VAR(NAME,READ_TYPE,LOCAL_TYPE) \
f_in.seek(#NAME); \
LOCAL_TYPE NAME; \
{ \
  READ_TYPE temp; \
  f_in>>temp; \
  NAME=to_ ## LOCAL_TYPE(temp); \
}

int main(
  int argc,
  char *argv[]
) {
  if (argc<2) {
    cout << "Must specify input data filename" << endl;
    exit(-1);
  }

  uint32 failed=0;

  it_ifile f_in(argv[1]);

  // Read input arguments
  READ_VAR(peaks_pow,vec,vec);
  READ_VAR(peaks_ind,ivec,ivec);
  peaks_ind-=1;
  READ_VAR(peaks_freq,ivec,vec);
  READ_VAR(peaks_n_id_2,ivec,ivec);
  READ_VAR(capbuf,cvec,cvec);
  READ_VAR(thresh2_n_sigma,ivec,ivec);
  READ_VAR(fc,ivec,ivec);
  // Read expected outputs
  READ_VAR(peaks_out_n_id_1,vec,vec);
  READ_VAR(peaks_out_frame_start,vec,vec);
  peaks_out_frame_start-=1;
  READ_VAR(peaks_out_cp_type,bvec,ivec);
  READ_VAR(peaks_out_freq_fine,vec,vec);
  READ_VAR(sss_h1_np_est,mat,mat);
  READ_VAR(sss_h2_np_est,mat,mat);
  READ_VAR(sss_h1_nrm_est,cmat,cmat);
  READ_VAR(sss_h2_nrm_est,cmat,cmat);
  READ_VAR(sss_h1_ext_est,cmat,cmat);
  READ_VAR(sss_h2_ext_est,cmat,cmat);

  f_in.close();

  uint16 n_cells=length(peaks_pow);
  for (uint16 t=0;t<n_cells;t++) {
    Cell cell_in;
    cell_in.pss_pow=peaks_pow[t];
    cell_in.ind=peaks_ind[t];
    cell_in.freq=peaks_freq[t];
    cell_in.n_id_2=peaks_n_id_2[t];

    vec sss_h1_np_est_meas;
    vec sss_h2_np_est_meas;
    cvec sss_h1_nrm_est_meas;
    cvec sss_h2_nrm_est_meas;
    cvec sss_h1_ext_est_meas;
    cvec sss_h2_ext_est_meas;
    mat log_lik_nrm;
    mat log_lik_ext;

    Cell cell_out;
    cell_out=sss_detect(cell_in,capbuf,thresh2_n_sigma[0],fc[0],sss_h1_np_est_meas,sss_h2_np_est_meas,sss_h1_nrm_est_meas,sss_h2_nrm_est_meas,sss_h1_ext_est_meas,sss_h2_ext_est_meas,log_lik_nrm,log_lik_ext);

    failed+=
      (max(abs(sss_h1_np_est_meas-sss_h1_np_est.get_row(t)))>1e-12)||
      (max(abs(sss_h2_np_est_meas-sss_h2_np_est.get_row(t)))>1e-12)||
      (max(abs(sss_h1_nrm_est_meas-sss_h1_nrm_est.get_row(t)))>1e-12)||
      (max(abs(sss_h2_nrm_est_meas-sss_h2_nrm_est.get_row(t)))>1e-12)||
      (max(abs(sss_h1_ext_est_meas-sss_h1_ext_est.get_row(t)))>1e-12)||
      (max(abs(sss_h2_ext_est_meas-sss_h2_ext_est.get_row(t)))>1e-12)
    ;
    if (isfinite(peaks_out_n_id_1(t))) {
      failed+=(cell_out.n_id_1!=peaks_out_n_id_1(t));
      failed+=(cell_out.cp_type!=((peaks_out_cp_type(t)==0)?cp_type_t::NORMAL:cp_type_t::EXTENDED));
      failed+=isnan(cell_out.frame_start)||(abs(cell_out.frame_start-peaks_out_frame_start(t))>1e-6);
    } else {
      failed+=(cell_out.n_id_1!=-1);
      failed+=(cell_out.cp_type!=cp_type_t::UNKNOWN);
      failed+=!isnan(cell_out.frame_start);
    }

    if (cell_out.n_id_1>=0) {
      Cell cell_out2;
      cell_out2=pss_sss_foe(cell_out,capbuf,fc[0]);
      failed+=(abs(peaks_out_freq_fine(t)-cell_out2.freq_fine)>1e-8);
    }
  }

  if (failed) {
    cout << "FAILED!!!" << endl;
  } else {
    cout << "passed" << endl;
  }

  return failed;
}

