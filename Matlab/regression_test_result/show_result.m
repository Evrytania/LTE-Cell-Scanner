% function show_result
clear all; close all;

format_str = {'r', 'b.', 'gs', 'k'};
resv_set = [1 2 4];
num_id = zeros(length(resv_set), 33);
figure;
for i = 1 : length(resv_set)
    load(['baseline/CellSearch_test1to33_twist0_fo-100to100_resv' num2str(resv_set(i)) '_numPtry3_Prange160.mat']);
    for j = 1 : length(cell_info)
        num_id(i,j) = length(cell_info{j});
    end
    plot(num_id(i,:), format_str{i}); hold on;
end
cell_info1 = cell_info;

figure;
for i = 1 : length(resv_set)
    load(['simple_resv1_2_4_pss_period_try1/CellSearch_test1to33_twist0_fo-140to140_resv' num2str(resv_set(i)) '_numPtry1_Prange160_parTh8.5_numPth0.5.mat']);
    for j = 1 : length(cell_info)
        num_id(i,j) = length(cell_info{j});
    end
    plot(num_id(i,:), format_str{i}); hold on;
end
cell_info2 = cell_info;
