
clear all; close all;

format_str = {'r', 'b.', 'gs', 'k'};
resv_set = [1 2 4];

subplot(4,1,1);
num_id = zeros(length(resv_set), 33);
for i = 1 : length(resv_set)
    load(['baseline/CellSearch_test1to33_twist0_fo-100to100_resv' num2str(resv_set(i)) '_numPtry3_Prange160.mat']);
    for j = 1 : length(cell_info)
        num_id(i,j) = length(cell_info{j});
    end
    plot(num_id(i,:), format_str{i}); hold on;
end
legend('rsv1', 'rsv2', 'rsv4');
cell_info1 = cell_info;

subplot(4,1,2);
num_id = zeros(length(resv_set), 33);
for i = 1 : length(resv_set)
    load(['baseline_num_peak_th0.5/CellSearch_test1to33_twist0_fo-140to140_resv' num2str(resv_set(i)) '_numPtry3_Prange160_parTh8.5_numPth0.5.mat']);
    for j = 1 : length(cell_info)
        a = cell_info{j};
        for k = 1 : length(a)
            if ~isnan(a(k).n_id_cell)
                num_id(i,j) = num_id(i,j) + 1;
            end
        end
    end
    plot(num_id(i,:), format_str{i}); hold on;
end
legend('rsv1', 'rsv2', 'rsv4');
cell_info2 = cell_info;

subplot(4,1,3);
num_id = zeros(length(resv_set), 33);
for i = 1 : length(resv_set)
    load(['simple_resv1_2_4_pss_period_try1/CellSearch_test1to33_twist0_fo-140to140_resv' num2str(resv_set(i)) '_numPtry1_Prange160_parTh8.5_numPth0.5.mat']);
    for j = 1 : length(cell_info)
        num_id(i,j) = length(cell_info{j});
    end
    plot(num_id(i,:), format_str{i}); hold on;
end
legend('rsv1', 'rsv2', 'rsv4');
cell_info3 = cell_info;

subplot(4,1,4);
num_id = zeros(length(resv_set), 33);
max_par = NaN(length(resv_set), 33);
min_par = NaN(length(resv_set), 33);
for i = 1 : length(resv_set)
    load(['simple_resv1_2_4_pss_period_try1_max_peak_all_improved/CellSearch_test1to33_twist0_fo-140to140_resv' num2str(resv_set(i)) '_numPtry1_Prange160_parTh8.5_numPth0.5.mat']);
    for j = 1 : length(cell_info)
        a = cell_info{j};
        for k = 1 : length(a)
            if ~isnan(a(k).n_id_cell)
                num_id(i,j) = num_id(i,j) + 1;
            end
        end
        if ~isempty(a)
            max_par(i,j) = max(a(1).extra_info.par);
            min_par(i,j) = min(a(1).extra_info.par);
        end
    end
    plot(num_id(i,:), format_str{i}); hold on;
end
legend('rsv1', 'rsv2', 'rsv4');
cell_info4 = cell_info;

figure;
subplot(2,1,1); plot(max_par(1,:),'b.-'); hold on; plot(max_par(2,:),'r.-'); hold on; plot(max_par(3,:),'k.-');
subplot(2,1,2); plot(min_par(1,:),'b.-'); hold on; plot(min_par(2,:),'r.-'); hold on; plot(min_par(3,:),'k.-');
