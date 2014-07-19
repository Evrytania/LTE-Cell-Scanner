
clear all; close all;

format_str = {'r', 'b.', 'gs', 'k'};
% load simple_resv1_2_4_pss_period_try1_max_peak_all_improved_record_more/CellSearch_test1to33_twist0_fo-140to140_resv2_numPtry1_Prange160_parTh15_numPth0.5.mat;
% load simple_resv1_2_4_pss_period_try1_max_peak_all_improved_record_more_ex_gain2.4/CellSearch_test1to33_twist0_fo-140to140_resv2_numPtry1_Prange160_parTh15_numPth0.5.mat;
load simple_resv1_2_4_pss_period_try1_max_peak_all_improved_record_more_ex_gain2/CellSearch_test1to33_twist0_fo-140to140_resv2_numPtry1_Prange160_parTh15_numPth0.5.mat;

idx_raw = [];
z_th1_store_raw = [];
pow_store_raw = [];
z_th1_store = [];
pow_store = [];
top_idx = [];
sub_idx = [];
for i = 1 : length(cell_info)
    a = cell_info{i};
    for k = 1 : length(a)
        z_th1_store_raw = [z_th1_store_raw a(k).Z_th1];
        pow_store_raw = [pow_store_raw a(k).pow];
        flag = 0;
        if ~isnan(a(k).n_id_cell)
            flag = 1;
            z_th1_store = [z_th1_store a(k).Z_th1];
            pow_store = [pow_store a(k).pow];
            top_idx = [top_idx i];
            sub_idx = [sub_idx k];
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
