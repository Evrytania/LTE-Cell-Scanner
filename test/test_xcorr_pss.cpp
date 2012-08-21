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

using namespace itpp;
using namespace std;

int main(
  int argc,
  char *argv[]
) {
  if (argc<2) {
    cout << "Must specify input data filename" << endl;
    exit(-1);
  }

  uint32 failed=0;

  // Read input arguments
  it_ifile f_in(argv[1]);

  f_in.seek("capbuf");
  cvec capbuf;
  f_in >> capbuf;

  f_in.seek("f_search_set");
  ivec f_search_set_temp;
  f_in >> f_search_set_temp;
  vec f_search_set=to_vec(f_search_set_temp);

  f_in.seek("ds_comb_arm");
  ivec ds_comb_arm_temp;
  f_in >> ds_comb_arm_temp;
  uint8 ds_comb_arm=ds_comb_arm_temp[0];

  ivec fc_temp;
  f_in.seek("fc");
  f_in >> fc_temp;
  double fc=fc_temp[0];

  // Read expected outputs
  f_in.seek("xc_incoherent_collapsed_pow");
  vec xc_incoherent_collapsed_pow_flat_expected;
  f_in >> xc_incoherent_collapsed_pow_flat_expected;

  f_in.seek("xc_incoherent_collapsed_frq");
  ivec xc_incoherent_collapsed_frq_flat_expected;
  f_in >> xc_incoherent_collapsed_frq_flat_expected;

  f_in.seek("xc_incoherent_single");
  vec xc_incoherent_single_flat_expected;
  f_in >> xc_incoherent_single_flat_expected;

  f_in.seek("xc_incoherent");
  vec xc_incoherent_flat_expected;
  f_in >> xc_incoherent_flat_expected;

  f_in.seek("sp_incoherent");
  vec sp_incoherent_expected;
  f_in >> sp_incoherent_expected;

  f_in.seek("xc");
  cvec xc_flat_expected;
  f_in >> xc_flat_expected;

  f_in.seek("sp");
  vec sp_expected;
  f_in >> sp_expected;

  f_in.close();

  // Call xcorr_pss
  mat xc_incoherent_collapsed_pow;
  imat xc_incoherent_collapsed_frq;
  vf3d xc_incoherent_single;
  vf3d xc_incoherent;
  vec sp_incoherent;
  vcf3d xc;
  vec sp;
  uint16 n_comb_xc;
  uint16 n_comb_sp;
  xcorr_pss(capbuf,f_search_set,ds_comb_arm,fc,xc_incoherent_collapsed_pow,xc_incoherent_collapsed_frq,xc_incoherent_single,xc_incoherent,sp_incoherent,xc,sp,n_comb_xc,n_comb_sp);

  cvec xc_flat=itpp_ext::flatten(xc);
  failed+=max(abs(xc_flat-xc_flat_expected))>1e-6;

  failed+=max(abs(sp-sp_expected.left(length(sp))))>1e-14;

  failed+=max(abs(sp_incoherent-sp_incoherent_expected))>1e-15;

  vec xc_incoherent_single_flat=itpp_ext::flatten(xc_incoherent_single);
  failed+=max(abs(xc_incoherent_single_flat-xc_incoherent_single_flat_expected))>1e-7;

  vec xc_incoherent_flat=itpp_ext::flatten(xc_incoherent);
  failed+=max(abs(xc_incoherent_flat-xc_incoherent_flat_expected))>1e-8;

  vec xc_incoherent_collapsed_pow_flat=cvectorize(xc_incoherent_collapsed_pow);
  failed+=max(abs(xc_incoherent_collapsed_pow_flat-xc_incoherent_collapsed_pow_flat_expected))>1e-8;

  ivec xc_incoherent_collapsed_frq_flat=cvectorize(xc_incoherent_collapsed_frq);
  failed+=max(abs(xc_incoherent_collapsed_frq_flat-(xc_incoherent_collapsed_frq_flat_expected-1)))!=0;

  if (failed) {
    cout << "FAILED!!!" << endl;
  } else {
    cout << "passed" << endl;
  }

  return failed;
}

