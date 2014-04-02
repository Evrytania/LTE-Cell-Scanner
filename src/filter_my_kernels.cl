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

// kernels for FIR 6RB channel filter.
// !!!This file should be put into $PATH directory or any other location where program can discover in runtime!!!

// filter coefficients // original order actuall is the same with reversed order
__constant float chn_6RB_filter_coef[47] = { \
8.193313185354206e-04,     3.535548569572820e-04,    -1.453429245341695e-03,     1.042805860697287e-03,     1.264224526451337e-03, \
  -3.219586065044259e-03,     1.423981657254563e-03,     3.859884310477692e-03,    -6.552708013395765e-03,     8.590509694961493e-04, \
  9.363722386299336e-03,    -1.120357391780316e-02,    -2.423088424232164e-03,     1.927528718829535e-02,    -1.646405738285926e-02, \
  -1.143040384534755e-02,     3.652830082843752e-02,    -2.132986170036144e-02,    -3.396829121834471e-02,     7.273086636811442e-02, \
  -2.476823886110626e-02,    -1.207789042999466e-01,     2.861583432079335e-01,     6.398255789896659e-01,     2.861583432079335e-01, \
  -1.207789042999466e-01,    -2.476823886110626e-02,     7.273086636811442e-02,    -3.396829121834471e-02,    -2.132986170036144e-02, \
  3.652830082843752e-02,    -1.143040384534755e-02,    -1.646405738285926e-02,     1.927528718829535e-02,    -2.423088424232164e-03, \
  -1.120357391780316e-02,     9.363722386299336e-03,     8.590509694961493e-04,    -6.552708013395765e-03,     3.859884310477692e-03, \
  1.423981657254563e-03,    -3.219586065044259e-03,     1.264224526451337e-03,     1.042805860697287e-03,    -1.453429245341695e-03, \
  3.535548569572820e-04,     8.193313185354206e-04
};


__kernel void skip2cols( __global float2* in,
                         __global float2* out,
                         uint len_in)
{// one work item per work group
  const size_t n = get_global_size(0);
  const size_t m = get_global_id(0);
  const size_t sub_len = len_in/n;
  const size_t base_idx = m*sub_len;

  size_t i, new_base_idx;
  for (i=0; i<sub_len; i++) {
    new_base_idx = i*n;
    out[new_base_idx + m] = in[base_idx + i];
  }
}

__kernel void multi_filter( __global float2* in,
                            __global float2* out,
                            uint len_in)
{// one work item per work group
  const size_t n = get_global_size(0);
  const size_t m = get_global_id(0);
  const size_t filter_len = sizeof(chn_6RB_filter_coef)/sizeof(float);
  const size_t sub_len_in = len_in/n;
  const size_t sub_len_out = sub_len_in + filter_len - 1;

  size_t i, j, base_idx, coef_idx;
  if (m==0){
    for (i=(sub_len_out-filter_len+1); i<sub_len_out; i++){
      base_idx = i*n;
      out[base_idx] = (float2)(0.0f, 0.0f);
    }
  }

  float2 acc;
  for (i=0; i<filter_len-1; i++){
    acc = (float2)(0.0f, 0.0f);

    for (j=0; j<i+1; j++) {
      base_idx = j*n;
      coef_idx = (filter_len-1)-i+j;
      acc = acc + in[base_idx + m] * chn_6RB_filter_coef[coef_idx];
    }

    base_idx = i*n;
    out[base_idx+m+1] = acc;
  }

  for (i=filter_len-1; i<=sub_len_in-1; i++){
    acc = (float2)(0.0f, 0.0f);

    for (j=0; j<filter_len; j++) {
      base_idx = (i-(filter_len-1)+j)*n;
      acc = acc + in[base_idx + m] * chn_6RB_filter_coef[j];
    }

    base_idx = i*n;
    out[base_idx+m+1] = acc;
  }

  for (i=sub_len_in; i<sub_len_out; i++){
    acc = (float2)(0.0f, 0.0f);

    for (j=0; j<sub_len_out-i; j++) {
      base_idx = (i-(filter_len-1)+j)*n;
      acc = acc + in[base_idx + m] * chn_6RB_filter_coef[j];
    }

    base_idx = i*n;
    out[base_idx+m+1] = acc;
  }
}

__kernel void result_combine( __global float2* in,
                              __global float2* out,
                              uint len_out)
{
  const size_t n = get_global_size(0);
  const size_t m = get_global_id(0);
  const size_t filter_len = sizeof(chn_6RB_filter_coef)/sizeof(float);
  const size_t sub_len_out = len_out/n;

  size_t base_linear_idx = m*sub_len_out;

  size_t i, base_idx, base_tail_idx;
  for (i=0; i<filter_len-1; i++){
    base_idx = i*n;
    base_tail_idx = (sub_len_out+i)*n;
    out[base_linear_idx+i] = in[base_idx+m+1] + in[base_tail_idx+m];
  }

  for (i=filter_len-1; i<sub_len_out; i++){
    base_idx = i*n;
    out[base_linear_idx+i] = in[base_idx+m+1];
  }
}

