
clear all; close all;

format_str = {'r', 'b.', 'gs', 'k'};
load par_explore/CellSearch_test1to33_twist0_fo-140to140_resv2_numPtry1_Prange-1_parTh8_numPth0.5.mat;

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
par_combined_max_valid = [];
par_combined_max_invalid = [];
par_max_max_valid = [];
par_max_max_invalid = [];
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
            par_combined_max_valid = [par_combined_max_valid a(k).extra_info.par_combined_max];
            par_max_max_valid = [par_max_max_valid a(k).extra_info.par_max_max];
            pow_store = [pow_store a(k).pow];
            top_idx = [top_idx i];
            sub_idx = [sub_idx k];
        else
            par_invalid = [par_invalid a(k).extra_info.par];
            par_combined_max_invalid = [par_combined_max_invalid a(k).extra_info.par_combined_max];
            par_max_max_invalid = [par_max_max_invalid a(k).extra_info.par_max_max];
        end
        idx_raw = [idx_raw flag];
    end
end

figure; plot(try_idx, 'b.-');

figure; plot(par_valid, 'b'); hold on; plot(par_invalid, 'r');
figure; plot(par_combined_max_valid, 'b'); hold on; plot(par_combined_max_invalid, 'r');
figure; plot(par_max_max_valid, 'b'); hold on; plot(par_max_max_invalid, 'r');
