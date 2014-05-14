% Jiao Xianjun (putaoshu@gmail.com; putaoshu@msn.com)
% LTE_DL_receiver.m

% Some scripts are borrowed from:
% https://github.com/JiaoXianjun/rtl-sdr-LTE
% https://github.com/Evrytania/LTE-Cell-Scanner
% https://github.com/Evrytania/Matlab-Library
% https://github.com/JiaoXianjun/multi-rtl-sdr-calibration

% See also README in root directory and ../scan-capture.

clear all;
close all;

% ------------------------------------------------------------------------------------
% % bin file captured by hackrf_transfer  
% filename = '../test/f2585_s19.2_bw20_1s_hackrf_bda.bin'; fc = 2585e6;
% filename = '../test/f2585_s19.2_bw20_1s_hackrf_bda1.bin'; fc = 2585e6;
% filename = '../test/f1860_s19.2_bw20_1s_hackrf_home1.bin'; fc = 1860e6;
% filename = '../test/f1860_s19.2_bw20_1s_hackrf_home.bin'; fc = 1860e6;
% filename = '../test/f1890_s19.2_bw20_1s_hackrf_home.bin'; fc = 1890e6;
% filename = '../test/f1890_s19.2_bw20_1s_hackrf_home1.bin'; fc = 1890e6;
filename = '../test/f2360_s19.2_bw20_1s_hackrf_bda.bin'; fc = 2360e6;

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

if isempty(dir([filename(1:end-4) '.mat']))
    r_raw = get_signal_from_bin(filename, inf);
    r_raw = r_raw - mean(r_raw); % remove DC

    r_pbch = filter_wo_tail(r_raw, coef_pbch.*5, sampling_rate_pbch/raw_sampling_rate);
    r_20M = filter_wo_tail(r_raw, coef_8x_up.*8, 8);
    r_20M = r_20M(1:5:end);
    
    plot(real(r_raw)); drawnow;
    [cell_info, r_pbch, r_20M] = CellSearch(r_pbch, r_20M, f_search_set, fc);
    
    r_pbch = r_pbch.';
    r_20M = r_20M.';
    save([filename(1:end-4) '.mat'], 'r_pbch', 'r_20M', 'cell_info');
else
    load([filename(1:end-4) '.mat']);
end

% cell_info

% for cell_idx = 1 : length(cell_info)
for cell_idx = 1 : 1
    cell_tmp = cell_info(cell_idx);
    [tfg, tfg_timestamp, cell_tmp]=extract_tfg(cell_tmp,r_20M,fc,sampling_carrier_twist, cell_tmp.n_rb_dl);
%     [tfg_comp, tfg_comp_timestamp, cell_tmp]=tfoec(cell_tmp, tfg, tfg_timestamp, fc, sampling_carrier_twist, cell_tmp.n_rb_dl);
%     cell_tmp=decode_mib(cell_tmp,tfg_comp(:, 565:636));
    
    n_symb_per_radioframe = 10*2*cell_tmp.n_symb_dl;
    num_radioframe = floor(size(tfg,1)/n_symb_per_subframe);
    num_subframe = num_radioframe*10;
    pdcch_info = cell(1, num_subframe);
    pcfich_info = zeros(1, num_subframe);
    pcfich_corr = zeros(1, num_subframe);

    % % ----------------following process radio frame by radio frame-------------------
    for radioframe_idx = 1 : num_radioframe
        
        subframe_base_idx = (radioframe_idx-1)*10;
        
        % % decode pcfich and identify uldl_cfg if TDD mode
        for subframe_idx = 1 : 10
            sp = subframe_base_idx + (subframe_idx-1)*n_symb_per_subframe + 1;
            ep = sp + n_symb_per_subframe - 1;

            [tfg_comp, tfg_comp_timestamp, cell_tmp]=tfoec_subframe(cell_tmp, subframe_idx-1, tfg(sp:ep, :), tfg_timestamp(sp:ep), fc, sampling_carrier_twist);
        end
        
        % % decode pdcch
        for subframe_idx = 1 : 10
            sp = (subframe_idx-1)*n_symb_per_subframe + 1;
            ep = sp + n_symb_per_subframe - 1;

            [tfg_comp, tfg_comp_timestamp, cell_tmp]=tfoec_subframe(cell_tmp, subframe_idx-1, tfg(sp:ep, :), tfg_timestamp(sp:ep), fc, sampling_carrier_twist);

            [pdcch_info{subframe_idx}, pcfich_info(subframe_idx), pcfich_corr(subframe_idx)] = decode_pdcch(cell_tmp, subframe_idx-1, tfg_comp);
        end
        
    end
end

disp(num2str(pcfich_corr));
sf_set = find(pcfich_info>0);
val_set = pcfich_info(pcfich_info>0);
disp(['subframe  ' num2str(sf_set)]);
disp(['num pdcch ' num2str(val_set)]);

subplot(3,1,1); plot(pcfich_corr); axis tight;
subplot(3,1,2); plot(sf_set, val_set, 'b.-'); axis tight;
subplot(3,1,3);
a = zeros(1, max(sf_set)); a(sf_set) = 1;
pcolor([a;a]); shading faceted;  axis tight;

