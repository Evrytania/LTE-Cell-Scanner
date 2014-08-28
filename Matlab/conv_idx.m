

function out_idx = conv_idx(in_idx, decimation_ratio)
out_idx = (in_idx-1).*decimation_ratio + 1;

