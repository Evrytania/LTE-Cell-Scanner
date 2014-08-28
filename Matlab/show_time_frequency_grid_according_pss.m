function show_time_frequency_grid_according_pss(pss_loc, k_factor, s)
clf;

s = s(:).';

peak_loc = find(pss_loc ~= -inf, 1, 'first');
peak_loc = pss_loc(peak_loc);

tdd_flag_set = [0 1 0 1];
cp_type_flag_set = [0 0 1 1];

for i = 1 : 4
    tdd_flag = tdd_flag_set(i);
    cp_type_flag = cp_type_flag_set(i);
    slot_start = get_slot_start(tdd_flag, cp_type_flag, peak_loc, k_factor);
    [tf_grid, sp_set] = get_time_frequency_grid(cp_type_flag, slot_start, k_factor, s);
    subplot(2,2,i); pcolor(tf_grid); shading flat; 
end

function slot_start = get_slot_start(tdd_flag, cp_type_flag, peak_loc, k_factor)
if tdd_flag==1
    if cp_type_flag == 0
        slot_start=peak_loc+(-(2*(128+9)+1)-2)*k_factor; % TDD NORMAL CP
    else
        slot_start=peak_loc+(-(2*(128+32))-2)*k_factor; % TDD EXTENDED CP
    end
else
    if cp_type_flag == 0
        slot_start=peak_loc+(-(6*(128+9)+1)-2)*k_factor; % FDD NORMAL CP
    else
        slot_start=peak_loc+(-(5*(128+32))-2)*k_factor; % FDD EXTENDED CP
    end
end

slot_start=wrap(slot_start, 0.5, 960+0.5);

decimation_ratio = 16;
slot_start = conv_idx(slot_start, decimation_ratio);

function [tf_grid, sp_set] = get_time_frequency_grid(cp_type_flag, slot_start, k_factor, s)

total_len = length(s);

decimation_ratio = 16;
len_symbol_core = 128*decimation_ratio;
if cp_type_flag == 0
    len_symbol = (128+9)*decimation_ratio;
else
    len_symbol = (128+32)*decimation_ratio;
end
len_cp = len_symbol - len_symbol_core;

num_symbol_rough = ceil(total_len/len_symbol);
tf_grid = zeros(1200, num_symbol_rough);
sp_set = zeros(1, num_symbol_rough);

symbol_count = 0;
symbol_start_tmp = slot_start;
symbol_round_idx = 0;
while symbol_start_tmp+len_symbol-1 <= total_len
    if cp_type_flag == 1
        len_cp_tmp = len_cp;
        len_symbol_tmp = len_symbol;
    else
        if symbol_round_idx == 0
            len_cp_tmp = len_cp + 16;
            len_symbol_tmp = len_symbol + 16;
        else
            len_cp_tmp = len_cp;
            len_symbol_tmp = len_symbol;
        end
        symbol_round_idx = mod(symbol_round_idx+1, 7);
    end

    sp = round(symbol_start_tmp + len_cp_tmp*k_factor);
    ep = sp + len_symbol_core - 1;
    symbol_start_tmp = symbol_start_tmp + len_symbol_tmp*k_factor;

    sym = s(sp:ep);
    fft_sym = fft(sym);
    spec = abs([fft_sym(end-600+1 : end)  fft_sym(1:600)]).';
    
    tf_grid(:, symbol_count+1) = spec;
    sp_set(symbol_count+1) = sp;
    
    symbol_count = symbol_count + 1;
end

tf_grid = tf_grid(:, 1:symbol_count);
sp_set = (1:symbol_count);
