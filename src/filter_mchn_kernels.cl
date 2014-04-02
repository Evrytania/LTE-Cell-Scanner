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
                            __global float2* coef,
                            uint len_in,
                            uint len_coef)
{// one work item per work group; the second dimension for each channels
  const size_t n = get_global_size(0);
  const size_t m = get_global_id(0);
  const size_t m1 = get_global_id(1);
  const size_t filter_len = len_coef;
  const size_t sub_len_in = len_in/n;
  const size_t sub_len_out = sub_len_in + filter_len - 1;
  const size_t base_coef = m1*len_coef;
  const size_t base_out = m1*sub_len_out*n;

  float2 coef_tmp;
  size_t i, j, base_idx, coef_idx;
  if (m==0){
    for (i=(sub_len_out-filter_len+1); i<sub_len_out; i++){
      base_idx = i*n;
      out[base_out+ base_idx] = (float2)(0.0f, 0.0f);
    }
  }

  for (i=0; i<filter_len-1; i++){
    float2 acc = (float2)(0.0f, 0.0f);

    for (j=0; j<i+1; j++) {
      base_idx = j*n;
      coef_idx = (filter_len-1)-i+j;
      coef_tmp = coef[base_coef+ coef_idx];
      acc = acc + (float2)( in[base_idx + m].x*coef_tmp.x - in[base_idx + m].y*coef_tmp.y,  in[base_idx + m].x*coef_tmp.y + in[base_idx + m].y*coef_tmp.x );
    }

    base_idx = i*n;
    out[base_out+ base_idx+m+1] = acc;
  }

  for (i=filter_len-1; i<=sub_len_in-1; i++){
    float2 acc = (float2)(0.0f, 0.0f);

    for (j=0; j<filter_len; j++) {
      base_idx = (i-(filter_len-1)+j)*n;
      coef_idx = j;
      coef_tmp = coef[base_coef+ coef_idx];
      acc = acc + (float2)( in[base_idx + m].x*coef_tmp.x - in[base_idx + m].y*coef_tmp.y,  in[base_idx + m].x*coef_tmp.y + in[base_idx + m].y*coef_tmp.x );
    }

    base_idx = i*n;
    out[base_out+ base_idx+m+1] = acc;
  }

  for (i=sub_len_in; i<sub_len_out; i++){
    float2 acc = (float2)(0.0f, 0.0f);

    for (j=0; j<sub_len_out-i; j++) {
      base_idx = (i-(filter_len-1)+j)*n;
      coef_idx = j;
      coef_tmp = coef[base_coef+ coef_idx];
      acc = acc + (float2)( in[base_idx + m].x*coef_tmp.x - in[base_idx + m].y*coef_tmp.y,  in[base_idx + m].x*coef_tmp.y + in[base_idx + m].y*coef_tmp.x );
    }

    base_idx = i*n;
    out[base_out+ base_idx+m+1] = acc;
  }
}

__kernel void result_combine( __global float2* in,
                              __global float2* out,
                              __global float* out_abs2,
                              uint len_out,
                              uint len_coef)
{
  const size_t n = get_global_size(0);
  const size_t m = get_global_id(0);
  const size_t m1 = get_global_id(1);
  const size_t filter_len = len_coef;
  const size_t sub_len_out = len_out/n;
  const size_t sub_len_in = sub_len_out + filter_len - 1;
  const size_t base_in = m1*sub_len_in*n;
  const size_t base_out = m1*len_out;

  size_t base_linear_idx = m*sub_len_out;

  size_t i, base_idx, base_tail_idx;
  float2 tmp_val;
  for (i=0; i<filter_len-1; i++){
    base_idx = i*n;
    base_tail_idx = (sub_len_out+i)*n;
    tmp_val = in[base_in+ base_idx+m+1] + in[base_in+ base_tail_idx+m];
    out[base_out+ base_linear_idx+i] = tmp_val;
    out_abs2[base_out+ base_linear_idx+i] = tmp_val.x*tmp_val.x + tmp_val.y*tmp_val.y;
  }

  for (i=filter_len-1; i<sub_len_out; i++){
    base_idx = i*n;
    tmp_val = in[base_in+ base_idx+m+1];
    out[base_out+ base_linear_idx+i] = tmp_val;
    out_abs2[base_out+ base_linear_idx+i] = tmp_val.x*tmp_val.x + tmp_val.y*tmp_val.y;
  }
}

