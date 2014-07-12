% Jiao Xianjun (putaoshu@gmail.com; putaoshu@msn.com)
% CellSearch.m
% Improved LTE-Cell-Scanner (written by James Peroulas: https://github.com/Evrytania/LTE-Cell-Scanner).
% See also README in root directory, ../test, ../../rtl-sdr-LTE/scan-capture/.

% Some scripts are borrowed from:
% https://github.com/JiaoXianjun/rtl-sdr-LTE
% https://github.com/Evrytania/LTE-Cell-Scanner
% https://github.com/Evrytania/Matlab-Library
% https://github.com/JiaoXianjun/multi-rtl-sdr-calibration

clear all;
close all;

[filename, file_sampling_rate] = regression_test_source(source_dir_idx, varargin);

sampling_carrier_twist = 0; % ATTENTION! If this is 1, make sure fc is aligned with bin file!!!

num_try = 10; % how many times we try for each frequency or file
num_radioframe = 8; % each radio frame length 10ms. MIB period is 4 radio frame

raw_sampling_rate = 19.2e6; % constrained by hackrf board
sampling_rate = 30.72e6;
sampling_rate_pbch = sampling_rate/16; % LTE spec. 30.72MHz/16.

num_subframe_per_radioframe = 10;
len_time_subframe = 1e-3; % 1ms. LTE spec
num_sample_per_radioframe = num_subframe_per_radioframe*len_time_subframe*sampling_rate_pbch;
num_sample_pbch = num_radioframe*num_sample_per_radioframe;

coef_pbch = fir1(254, (0.18e6*6+150e3)/raw_sampling_rate); %freqz(coef_pbch, 1, 1024);
coef_8x_up = fir1(254, 20e6/(raw_sampling_rate*8)); %freqz(coef_8x_up, 1, 1024);

DS_COMB_ARM = 2;
FS_LTE = 30720000;
thresh1_n_nines=12;
rx_cutoff=(6*12*15e3/2+4*15e3)/(FS_LTE/16/2);
THRESH2_N_SIGMA = 3;

f_search_set = 20e3:5e3:30e3; % change it wider if you don't know pre-information

    r_raw = get_signal_from_bin(filename, inf);
    r_raw = r_raw - mean(r_raw); % remove DC

    r_pbch = filter_wo_tail(r_raw, coef_pbch.*5, sampling_rate_pbch/raw_sampling_rate);
    r_20M = filter_wo_tail(r_raw, coef_8x_up.*8, 8);
    r_20M = r_20M(1:5:end);
    
    plot(real(r_raw)); drawnow;
    [cell_info, r_pbch, r_20M] = CellSearch(r_pbch, r_20M, f_search_set, fc);