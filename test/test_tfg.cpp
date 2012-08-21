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

// This program tests all of the key routines starting with the generation of
// the time/frequency grid tfg.
#include <itpp/itbase.h>
#include <list>
#include "rtl-sdr.h"
#include "common.h"
#include "macros.h"
#include "lte_lib.h"
#include "searcher.h"
#include "itpp_ext.h"

using namespace itpp;
using namespace std;

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
  READ_VAR(peaks_in_n_id_1,ivec,ivec);
  READ_VAR(peaks_in_n_id_2,ivec,ivec);
  READ_VAR(peaks_in_cp_type,ivec,ivec);
  READ_VAR(peaks_in_frame_start,vec,vec);
  peaks_in_frame_start-=1;
  READ_VAR(peaks_in_freq_fine,vec,vec);
  READ_VAR(capbuf,cvec,cvec);
  READ_VAR(fc,ivec,vec);
  // Read expected outputs
  READ_VAR(tfg,cmat,cmat);
  READ_VAR(tfg_timestamp,vec,vec);
  tfg_timestamp-=1;
  READ_VAR(tfg_comp,cmat,cmat);
  READ_VAR(tfg_comp_timestamp,vec,vec);
  tfg_comp_timestamp-=1;
  READ_VAR(peaks_out_freq_superfine,vec,vec);
  READ_VAR(peaks_out_n_rb_dl,ivec,ivec);
  READ_VAR(peaks_out_phich_dur,ivec,ivec);
  READ_VAR(peaks_out_phich_res,vec,vec);
  READ_VAR(peaks_out_sfn,ivec,ivec);
  f_in.close();

  Cell cell;
  cell.n_id_1=peaks_in_n_id_1[0];
  cell.n_id_2=peaks_in_n_id_2[0];
  cell.cp_type=peaks_in_cp_type[0]?cp_type_t::EXTENDED:cp_type_t::NORMAL;
  cell.frame_start=peaks_in_frame_start[0];
  cell.freq_fine=peaks_in_freq_fine[0];

  // Call the functions to be tested and check the outputs.
  cmat tfg_actual;
  vec tfg_timestamp_actual;
  extract_tfg(cell,capbuf,fc[0],tfg_actual,tfg_timestamp_actual);
  failed+=(max(max(abs(tfg-tfg_actual)))>1e-10);
  failed+=(max(abs(tfg_timestamp-tfg_timestamp_actual))>1e-10);

  RS_DL rs_dl(cell.n_id_cell(),6,cell.cp_type);

  cmat tfg_comp_actual;
  vec tfg_comp_timestamp_actual;
  Cell cell_out_tfg=tfoec(cell,tfg_actual,tfg_timestamp_actual,fc[0],rs_dl,tfg_comp_actual,tfg_comp_timestamp_actual);
  failed+=(max(max(abs(tfg_comp-tfg_comp_actual)))>1e-10);
  failed+=(max(abs(tfg_comp_timestamp-tfg_comp_timestamp_actual))>1e-10);
  failed+=(abs(cell_out_tfg.freq_superfine-peaks_out_freq_superfine[0])>1e-7);

  Cell cell_out_dmib=decode_mib(cell_out_tfg,tfg_comp_actual,rs_dl);
  failed+=(cell_out_dmib.n_rb_dl!=50);

  if (failed) {
    cout << "FAILED!!!" << endl;
  } else {
    cout << "passed" << endl;
  }

  return failed;
}

