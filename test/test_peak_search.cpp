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
    cerr << "Must specify input data filename" << endl;
    exit(-1);
  }

  uint32 failed=0;

  // Read input arguments
  it_ifile f_in(argv[1]);

  READ_VAR(xc_incoherent_collapsed_pow,mat,mat);
  READ_VAR(xc_incoherent_collapsed_frq,imat,imat);
  // Convert from matlab 1 based indexing to c 0 based indexing.
  xc_incoherent_collapsed_frq-=1;
  READ_VAR(Z_th1,vec,vec);
  READ_VAR(f_search_set,ivec,vec);
  READ_VAR(peaks_pow,vec,vec);
  READ_VAR(peaks_ind,ivec,ivec);
  // Convert from matlab 1 based indexing to c 0 based indexing.
  peaks_ind-=1;
  READ_VAR(peaks_freq,ivec,ivec);
  READ_VAR(peaks_n_id_2,ivec,ivec);

  f_in.close();

  uint8 n_f=length(f_search_set);
  vf3d xc_incoherent_single=vector < vector < vector < float > > > (3,vector< vector < float > >(9600, vector < float > (n_f,NAN)));
  for (uint32 k=0;k<3;k++) {
    for (uint32 m=0;m<9600;m++) {
      for (uint8 t=0;t<n_f;t++) {
        xc_incoherent_single[k][m][t]=xc_incoherent_collapsed_pow(k,m);
      }
    }
  }

  // Call peak_search
  list <Cell> cells;
  peak_search(xc_incoherent_collapsed_pow,xc_incoherent_collapsed_frq,Z_th1,f_search_set,739e6,xc_incoherent_single,0,cells);

  // Check results
  if (cells.size()!=(unsigned)peaks_pow.length()) {
    failed=-1;
  } else {
    uint16 offset=0;
    for (list<Cell>::iterator iterator=cells.begin();iterator!=cells.end();iterator++) {
      failed+=
        ((*iterator).pss_pow-peaks_pow[offset]>1e-6)||
        ((*iterator).ind!=peaks_ind[offset])||
        ((*iterator).freq!=peaks_freq[offset])||
        ((*iterator).n_id_2!=peaks_n_id_2[offset])
      ;
      offset++;
    }
  }

  if (failed) {
    cout << "FAILED!!!" << endl;
  } else {
    cout << "passed" << endl;
  }

  return failed;
}

