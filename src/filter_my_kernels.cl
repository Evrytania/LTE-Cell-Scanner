__kernel void filter1( __global float* in_buf,
                       __global float* out_buf,
                       size_t len_in)
{
  const size_t n = get_global_size(0);
  const size_t i = get_global_id(0);
  out_buf[i] = (float)( n+i );
}

__kernel void filter2( __global float* in_buf,
                       __global float* out_buf,
                       size_t len_in)
{
  const size_t n = get_global_size(0);
  const size_t i = get_global_id(0);
  out_buf[i] = (float)( 2*(n+i) );
}

