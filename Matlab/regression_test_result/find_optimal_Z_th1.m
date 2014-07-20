
clear all; close all;

format_str = {'r', 'b.', 'gs', 'k'};
% load simple_resv1_2_4_pss_period_try1_max_peak_all_improved_record_more/CellSearch_test1to33_twist0_fo-140to140_resv2_numPtry1_Prange160_parTh15_numPth0.5.mat;
% load simple_resv1_2_4_pss_period_try1_max_peak_all_improved_record_more_ex_gain2.4/CellSearch_test1to33_twist0_fo-140to140_resv2_numPtry1_Prange160_parTh15_numPth0.5.mat;
% load simple_resv1_2_4_pss_period_try1_max_peak_all_improved_record_more_ex_gain2/CellSearch_test1to33_twist0_fo-140to140_resv2_numPtry1_Prange160_parTh15_numPth0.5.mat;
load simple_resv1_2_4_pss_period_try1_max_peak_all_improved_record_more_ex_gain2_peak_range-1/CellSearch_test1to33_twist0_fo-140to140_resv2_numPtry1_Prange-1_parTh8_numPth0.5.mat;

idx_raw = [];
z_th1_store_raw = [];
pow_store_raw = [];
z_th1_store = [];
pow_store = [];
top_idx = [];
sub_idx = [];
try_idx = [];
par_valid = [];
par_invalid = [];
for i = 1 : length(cell_info)
    a = cell_info{i};
    for k = 1 : length(a)
        z_th1_store_raw = [z_th1_store_raw a(k).Z_th1];
        pow_store_raw = [pow_store_raw a(k).pow];
        flag = 0;
        if ~isnan(a(k).n_id_cell)
            flag = 1;
            z_th1_store = [z_th1_store a(k).Z_th1];
            try_idx = [try_idx a(k).extra_info.try_idx];
            par_valid = [par_valid a(k).extra_info.par];
            pow_store = [pow_store a(k).pow];
            top_idx = [top_idx i];
            sub_idx = [sub_idx k];
        else
            par_invalid = [par_invalid a(k).extra_info.par];
        end
        idx_raw = [idx_raw flag];
    end
end

figure;
plot(z_th1_store); hold on;
plot(pow_store, 'r');

figure;
plot(pow_store./z_th1_store);

figure;
plot(pow_store_raw./z_th1_store_raw); hold on;

idx = find(idx_raw);
plot(idx, pow_store_raw(idx)./z_th1_store_raw(idx), 'r.'); hold on;

figure; plot(try_idx, 'b.-');

figure; plot(par_valid, 'b'); hold on; plot(par_invalid, 'r');

num_id = zeros(1, 33);
min_par = NaN(1, 33);
for j = 1 : length(cell_info)
    a = cell_info{j};
    for k = 1 : length(a)
        if ~isnan(a(k).n_id_cell)
            num_id(j) = num_id(j) + 1;
            min_par(j) = min(a(1).extra_info.par);
        end
    end
end
figure;
plot(num_id, 'b.-'); hold on;
load simple_resv1_2_4_pss_period_try1_max_peak_all_improved_record_more_ex_gain2/CellSearch_test1to33_twist0_fo-140to140_resv2_numPtry1_Prange160_parTh15_numPth0.5.mat;
num_id = zeros(1, 33);
for j = 1 : length(cell_info)
    a = cell_info{j};
    for k = 1 : length(a)
        if ~isnan(a(k).n_id_cell)
            num_id(j) = num_id(j) + 1;
        end
    end
end
plot(num_id, 'rs-'); hold on;

figure; plot(min_par);
