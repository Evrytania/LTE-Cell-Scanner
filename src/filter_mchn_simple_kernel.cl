// An OpenCL accelerated LTE Cell Scanner
//
// Written by Jiao Xianjun <putaoshu@gmail.com>
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

// kernels for multi channels PSS correlation
// !!!This file should be put into $PATH directory or any other location where program can discover in runtime!!!
// This is a simplified version of filter_mchn_kernels.cl to test speedup.

#define CAPLENGTH 153600  // remember to conform with main code
#define COEF_LEN 137  // remember to conform with main code
#define OUT_LEN (CAPLENGTH-COEF_LEN+1)  // remember to conform with main code

__kernel void multi_filter( __global float2* in,
                            __global float* out,
                            __global float2* coef
                            )
{// one work item per work group; the second dimension for each channels
  const size_t m = get_global_id(0);
  const size_t base_out = m*OUT_LEN;
  const size_t base_coef = m*COEF_LEN;

  float2 acc, in_tmp, coef_tmp;
  size_t i, j;
  for (i=0; i<OUT_LEN; i++){
    acc = (float2)(0.0f, 0.0f);

    for (j=0; j<COEF_LEN; j++) {
      coef_tmp = coef[base_coef+ j];
      in_tmp = in[i+j];
      acc = acc + (float2)( in_tmp.x*coef_tmp.x - in_tmp.y*coef_tmp.y,  in_tmp.x*coef_tmp.y + in_tmp.y*coef_tmp.x );
    }

    out[base_out+i] = acc.x*acc.x + acc.y*acc.y;
  }

}
