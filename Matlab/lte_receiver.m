% Jiao Xianjun (putaoshu@gmail.com)
% lte_receiver.m
% Some are from updated https://github.com/JiaoXianjun/rtl-sdr-LTE
clear all; close all;

% filename = '../test/f1860_s15.36_bw10_l32_g34_1s.bin'; % FDD 20MHz
% filename = '../test/f1860_s15.36_bw10_l32_g40_1s.bin'; % FDD 20MHZ
% filename = '../test/f1860_s15.36_bw10_l24_g40_1s.bin'; % FDD 20MHZ
% filename = '../test/f1860_s15.36_bw10_l40_g32_1s.bin'; % FDD 20MHZ
% filename = '../test/f1860_s15.36_bw10_l40_g40_1s.bin'; % FDD 20MHZ

% filename = '../test/f1890_s15.36_bw10_l32_g34_1s.bin'; % TDD 20MHz
% filename = '../test/f1890_s15.36_bw10_l32_g32_1s.bin'; % TDD 20MHz
% filename = '../test/f1890_s15.36_bw10_l32_g40_1s.bin'; % TDD 20MHz
% filename = '../test/f1890_s15.36_bw10_l32_g22_1s.bin'; % TDD 20MHz

filename = '../test/f1860_s15.36_bw10_l32_g36_1s.bin'; % FDD 20MHz
% filename = '../test/f1890_s15.36_bw10_l32_g36_1s.bin'; % TDD 20MHz

sampling_rate = 15.36e6;
coef = fir1(158, ((9e6+400e3))/sampling_rate); %freqz(coef, 1, 1024); %10M channel filter. 

s = get_signal_from_bin(filename, 100e-3*sampling_rate);

s = filter_wo_tail(s, coef, 1);

f_search_set = 20e3:5e3:30e3; % change it wider if you don't know pre-information
[peaks, s] = pss_sss_detect(s, sampling_rate, f_search_set);
