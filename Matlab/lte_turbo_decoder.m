function bits = lte_turbo_decoder(llr, niter)
code_length = (length(llr)-12)/3;

if code_length ~= floor(code_length)
    disp('Abnormal code length!');
    return;
end

[info, parity0, parity1] = demultiplex(llr, code_length);

[end_state_conv0, end_state_conv1] = pre_trace_back(info, parity0, parity1, code_length);
intlv_table = intlvLTE(code_length);

inf_val = 8192;
info1 = info(intlv_table);
llr_ex = zeros(1, code_length);
for iter_idx = 1 : (niter-1)
    [llr_ex, ~] = trellis_proc(llr_ex, [0, -inf_val.*ones(1,7)], end_state_conv0, info, parity0, code_length, 0);
    llr_ex = llr_ex(intlv_table);
    [llr_ex, ~] = trellis_proc(llr_ex, [0, -inf_val.*ones(1,7)], end_state_conv1, info1, parity1, code_length, 0);
    llr_ex(intlv_table) = llr_ex;
end
[llr_ex, ~] = trellis_proc(llr_ex, [0, -inf_val.*ones(1,7)], end_state_conv0, info, parity0, code_length, 0);
llr_ex = llr_ex(intlv_table);
[~, bits] = trellis_proc(llr_ex, [0, -inf_val.*ones(1,7)], end_state_conv1, info1, parity1, code_length, 1);
bits(intlv_table) = bits;

function [llr_ex, bits] = trellis_proc(llr_pre, init_state_forward, init_state_backward, info, parity, code_length, bits_out_flag)

bits = inf(1, code_length);
llr_ex = zeros(1, code_length);

num_state = length(init_state_forward);

% forward transverse
state_metric = init_state_forward;
state_metric_store = zeros(code_length, num_state);
for i = 1 : code_length
    state_metric_store(i,:) = state_metric;
    
    tmp_info = info(i) + llr_pre(i);
    
    path_metric(1) = -tmp_info - parity(i);
    path_metric(2) = -tmp_info + parity(i);
    path_metric(3) = -path_metric(2);
    path_metric(4) = -path_metric(1);
    
    state_metric = state_proc(state_metric, path_metric);
end

% backward transverse
state_metric = init_state_backward;
if bits_out_flag == 0
    for i = code_length : -1 : 1
        tmp_info = info(i) + llr_pre(i);

        path_metric(1) = -tmp_info - parity(i);
        path_metric(2) = -tmp_info + parity(i);
        path_metric(3) = -path_metric(2);
        path_metric(4) = -path_metric(1);

        llr_ex_tmp = bit_llr_proc(state_metric_store(i,:), state_metric, parity(i));
        llr_ex(i) = 0.7*llr_ex_tmp;
        
        tmp_state_in = [state_metric(1), state_metric(5), state_metric(3), state_metric(7), state_metric(2), state_metric(6), state_metric(4), state_metric(8)];
        tmp_state_out = state_proc(tmp_state_in, [path_metric(1) path_metric(3) path_metric(2) path_metric(4)]);
        state_metric = [tmp_state_out(1) tmp_state_out(5) tmp_state_out(3) tmp_state_out(7) tmp_state_out(2) tmp_state_out(6) tmp_state_out(4) tmp_state_out(8)];
    end
else
    for i = code_length : -1 : 1
        tmp_info = info(i) + llr_pre(i);

        path_metric(1) = -tmp_info - parity(i);
        path_metric(2) = -tmp_info + parity(i);
        path_metric(3) = -path_metric(2);
        path_metric(4) = -path_metric(1);

        llr_ex_tmp = bit_llr_proc(state_metric_store(i,:), state_metric, parity(i));
        bits(i) = (1 + sign(tmp_info + llr_ex_tmp))/2;
        
        tmp_state_in = [state_metric(1), state_metric(5), state_metric(3), state_metric(7), state_metric(2), state_metric(6), state_metric(4), state_metric(8)];
        tmp_state_out = state_proc(tmp_state_in, [path_metric(1) path_metric(3) path_metric(2) path_metric(4)]);
        state_metric = [tmp_state_out(1) tmp_state_out(5) tmp_state_out(3) tmp_state_out(7) tmp_state_out(2) tmp_state_out(6) tmp_state_out(4) tmp_state_out(8)];
    end
end

function llr = bit_llr_proc(start_state, end_state, parity_bit)
% % bit 1
tmp1 = start_state(1) + end_state(1);
tmp2 = start_state(2) + end_state(5);
s0 = max(tmp1, tmp2);

tmp1 = start_state(7) + end_state(8);
tmp2 = start_state(8) + end_state(4);
s1 = max(tmp1, tmp2);

tmp1 = start_state(3) + end_state(6);
tmp2 = start_state(4) + end_state(2);
s2 = max(tmp1, tmp2);

tmp1 = start_state(5) + end_state(3);
tmp2 = start_state(6) + end_state(7);
s3 = max(tmp1, tmp2);

tmp1 = max(s0, s1) - parity_bit;
tmp2 = max(s2, s3) + parity_bit;
L_bit1 = max(tmp1, tmp2);

% % bit 0
tmp1 = start_state(1) + end_state(5);
tmp2 = start_state(2) + end_state(1);
s0 = max(tmp1, tmp2);

tmp1 = start_state(7) + end_state(4);
tmp2 = start_state(8) + end_state(8);
s1 = max(tmp1, tmp2);

tmp1 = start_state(3) + end_state(2);
tmp2 = start_state(4) + end_state(6);
s2 = max(tmp1, tmp2);

tmp1 = start_state(5) + end_state(7);
tmp2 = start_state(6) + end_state(3);
s3 = max(tmp1, tmp2);

tmp1 = max(s0, s1) + parity_bit;
tmp2 = max(s2, s3) - parity_bit;
L_bit0 = max(tmp1, tmp2);

llr = (L_bit0 - L_bit1)/2;

function [end_state_conv0, end_state_conv1] = pre_trace_back(info, parity0, parity1, code_length)
inf_val = 8192;
end_state_conv0 = -inf_val.*ones(1, 8);
end_state_conv1 = -inf_val.*ones(1, 8);
end_state_conv0(1) = 0;
end_state_conv1(1) = 0;

path_metric = zeros(1, 4);
for i = code_length+3 : -1 : code_length+1
    % ----------conv0--------------
    path_metric(1) = -info(i) - parity0(i);
    path_metric(2) = -info(i) + parity0(i);
    path_metric(3) = -path_metric(2);
    path_metric(4) = -path_metric(1);
    
    tmp_state_in = [end_state_conv0(1), end_state_conv0(5), end_state_conv0(3), end_state_conv0(7), end_state_conv0(2), end_state_conv0(6), end_state_conv0(4), end_state_conv0(8)];
    tmp_state_out = state_proc(tmp_state_in, [path_metric(1) path_metric(3) path_metric(2) path_metric(4)]);
    end_state_conv0 = [tmp_state_out(1) tmp_state_out(5) tmp_state_out(3) tmp_state_out(7) tmp_state_out(2) tmp_state_out(6) tmp_state_out(4) tmp_state_out(8)];
    
    % ----------conv1--------------
    path_metric(1) = -info(i+3) - parity1(i);
    path_metric(2) = -info(i+3) + parity1(i);
    path_metric(3) = -path_metric(2);
    path_metric(4) = -path_metric(1);
    
    tmp_state_in = [end_state_conv1(1), end_state_conv1(5), end_state_conv1(3), end_state_conv1(7), end_state_conv1(2), end_state_conv1(6), end_state_conv1(4), end_state_conv1(8)];
    tmp_state_out = state_proc(tmp_state_in, [path_metric(1) path_metric(3) path_metric(2) path_metric(4)]);
    end_state_conv1 = [tmp_state_out(1) tmp_state_out(5) tmp_state_out(3) tmp_state_out(7) tmp_state_out(2) tmp_state_out(6) tmp_state_out(4) tmp_state_out(8)];
end

function [info, parity0, parity1] = demultiplex(llr, code_length)
parity0 = zeros(1, code_length+3);
parity1 = zeros(1, code_length+3);
info = zeros(1, code_length+6);

info(1:code_length) = llr(1 : 3 : (3*code_length) );
parity0(1:code_length) = llr(2 : 3 : (3*code_length) );
parity1(1:code_length) = llr(3 : 3 : (3*code_length) );

info(code_length+1 : code_length+3) = llr((3*code_length+1) : 2 : (3*code_length+6) );
parity0(code_length+1 : code_length+3) = llr((3*code_length+2) : 2 : (3*code_length+6) );

info(code_length+4 : code_length+6) = llr((3*code_length+7) : 2 : (3*code_length+12) );
parity1(code_length+1 : code_length+3) = llr((3*code_length+8) : 2 : (3*code_length+12) );

function state_metric_out = state_proc(state_metric_in, path_metric)

tmp1 = state_metric_in(1)  + path_metric(1);
tmp2 = state_metric_in(2)  + path_metric(4);
state_metric_out(1) = max(tmp1, tmp2);

tmp1 = state_metric_in(3)  + path_metric(3);
tmp2 = state_metric_in(4)  + path_metric(2);
state_metric_out(2) = max(tmp1, tmp2);

tmp1 = state_metric_in(5)  + path_metric(2);
tmp2 = state_metric_in(6)  + path_metric(3);
state_metric_out(3) = max(tmp1, tmp2);

tmp1 = state_metric_in(7)  + path_metric(4);
tmp2 = state_metric_in(8)  + path_metric(1);
state_metric_out(4) = max(tmp1, tmp2);

tmp1 = state_metric_in(1)  + path_metric(4);
tmp2 = state_metric_in(2)  + path_metric(1);
state_metric_out(5) = max(tmp1, tmp2);

tmp1 = state_metric_in(3)  + path_metric(2);
tmp2 = state_metric_in(4)  + path_metric(3);
state_metric_out(6) = max(tmp1, tmp2);

tmp1 = state_metric_in(5)  + path_metric(3);
tmp2 = state_metric_in(6)  + path_metric(2);
state_metric_out(7) = max(tmp1, tmp2);

tmp1 = state_metric_in(7)  + path_metric(1);
tmp2 = state_metric_in(8)  + path_metric(4);
state_metric_out(8) = max(tmp1, tmp2);
