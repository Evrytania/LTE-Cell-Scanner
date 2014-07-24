% Jiao Xianjun (putaoshu@gmail.com; putaoshu@msn.com)
% CellSearch.m
% Improved LTE-Cell-Scanner (written by James Peroulas: https://github.com/Evrytania/LTE-Cell-Scanner).
% See also README in root directory, ../test, ../../rtl-sdr-LTE/scan-capture/.

% Some scripts are borrowed from:
% https://github.com/JiaoXianjun/rtl-sdr-LTE
% https://github.com/Evrytania/LTE-Cell-Scanner
% https://github.com/Evrytania/Matlab-Library
% https://github.com/JiaoXianjun/multi-rtl-sdr-calibration

% clear all;
% close all;
function cell_info = test_CellSearch(test_sp, test_ep)
test_source_info = regression_test_source('../regression_test_signal_file');

% test_sp = 1;
% test_ep = length(test_source_info);
% test_sp = 10; test_ep = 10;
sampling_carrier_twist = 0;
% f_search_set = -140e3:5e3:140e3;
f_search_set = -140e3:5e3:140e3; % align to C
pss_peak_max_reserve = 2;
num_pss_period_try = 1;
% combined_pss_peak_range = 160;
combined_pss_peak_range = -1; % set it to -1 to use complementary range of peak.
par_th = 8.5;
num_peak_th = 1/2; % originally is 2/3;

filename = ['CellSearch_test' num2str(test_sp) 'to' num2str(test_ep) '_twist' num2str(sampling_carrier_twist) '_fo' num2str(min(f_search_set)/1e3) 'to' num2str(max(f_search_set)/1e3) '_resv' num2str(pss_peak_max_reserve) '_numPtry' num2str(num_pss_period_try) '_Prange' num2str(combined_pss_peak_range) '_parTh' num2str(par_th) '_numPth' num2str(num_peak_th) '.mat'];

cell_info = cell(1, length(test_source_info));
for i = test_sp : test_ep
%     if isempty( strfind(test_source_info(i).filename, 'dimitri') )
%         continue;
%     end
    disp(test_source_info(i).filename);
    coef_pbch = pbch_filter_coef_gen(test_source_info(i).fs);
    
    r_raw = get_signal_from_bin(test_source_info(i).filename, inf, test_source_info(i).dev);
    r_raw = r_raw - mean(r_raw); % remove DC

    r_pbch = filter_wo_tail(r_raw, coef_pbch, (30.72e6/16)/test_source_info(i).fs);
    [~, ~, ~, cell_info{i}] = CellSearch(r_pbch, [], f_search_set, test_source_info(i).fc, sampling_carrier_twist, pss_peak_max_reserve, num_pss_period_try, combined_pss_peak_range, par_th, num_peak_th);
    save(filename, 'test_source_info', 'cell_info', 'test_sp', 'test_ep', 'sampling_carrier_twist', 'f_search_set', 'pss_peak_max_reserve', 'num_pss_period_try', 'combined_pss_peak_range');
end

