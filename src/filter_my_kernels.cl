constant float chn_6RB_filter_coef[47] = { \
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

  for (size_t i=0; i<sub_len; i++) {
    size_t new_base_idx = i*n;
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

  if (m==0){
    for (size_t i=(sub_len_out-filter_len+1); i<sub_len_out; i++){
      size_t base_idx = i*n;
      out[base_idx] = (float2)(0.0f, 0.0f);
    }
  }

  for (size_t i=0; i<filter_len-1; i++){
    float2 acc = (float2)(0.0f, 0.0f);

    for (size_t j=0; j<i+1; j++) {
      size_t base_idx = j*n;
      acc = acc + in[base_idx + m] * chn_6RB_filter_coef[i-j];
    }

    size_t base_idx = i*n;
    out[base_idx+m+1] = acc;
  }

//  for (size_t i=filter_len-1; i<=sub_len_out-filter_len; i++){
  for (size_t i=filter_len-1; i<=sub_len_in-1; i++){
    float2 acc = (float2)(0.0f, 0.0f);

    for (size_t j=0; j<filter_len; j++) {
      size_t base_idx = (i-(filter_len-1)+j)*n;
      acc = acc + in[base_idx + m] * chn_6RB_filter_coef[(filter_len-1)-j];
    }

    size_t base_idx = i*n;
    out[base_idx+m+1] = acc;
  }

//  for (size_t i=sub_len_out-filter_len+1; i<sub_len_out; i++){
  for (size_t i=sub_len_in; i<sub_len_out; i++){
    float2 acc = (float2)(0.0f, 0.0f);

//    for (size_t j=0; j<(filter_len- (i-(sub_len_out-filter_len))); j++) {
    for (size_t j=0; j<sub_len_out-i; j++) {
      size_t base_idx = (i-(filter_len-1)+j)*n;
      acc = acc + in[base_idx + m] * chn_6RB_filter_coef[(filter_len-1)-j];
    }

    size_t base_idx = i*n;
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

  for (size_t i=0; i<filter_len-1; i++){
    size_t base_idx = i*n;
    size_t base_tail_idx = (sub_len_out+i)*n;
    out[base_linear_idx+i] = in[base_idx+m+1] + in[base_tail_idx+m];
  }

  for (size_t i=filter_len-1; i<sub_len_out; i++){
    size_t base_idx = i*n;
    out[base_linear_idx+i] = in[base_idx+m+1];
  }
}

