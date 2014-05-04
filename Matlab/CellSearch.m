% Jiao Xianjun (putaoshu@gmail.com; putaoshu@msn.com)
% CellSearch.m
% Improved LTE-Cell-Scanner (written by James Peroulas: https://github.com/Evrytania/LTE-Cell-Scanner).
% Add 1) TD-LTE; 2) external mixer (no assumption on relationship between sampling and carrier error) support
% test HACKRF now.

% Some scripts are borrowed from:
% https://github.com/JiaoXianjun/rtl-sdr-LTE
% https://github.com/Evrytania/LTE-Cell-Scanner
% https://github.com/Evrytania/Matlab-Library
% https://github.com/JiaoXianjun/multi-rtl-sdr-calibration

% See also README in root directory and ../scan-capture.

clear all;
close all;

fc = 1890e6;
% fc = 1860e6;

% ------------------------------------------------------------------------------------
% % bin file captured by hackrf_transfer  
% filename = '../test/f1860_s20_bw20_l40_g40_1s_hackrf_home.bin';
% filename = '../test/f1890_s20_bw20_1s_hackrf_bda.bin';
% filename = '../test/f1890_s20_bw20_l32_g40_1s_hackrf_bda.bin';
% filename = '../test/f1890_s20_bw20_l32_g40_1s_hackrf_home.bin';
% filename = '../test/f2360_s20_bw20_1s_hackrf_bda.bin';
% filename = '../test/f2360_s20_bw20_l32_g40_1s_hackrf_bda.bin';
filename = '../test/f2585_s20_bw20_l40_g40_1s_hackrf_bda.bin';

sampling_carrier_twist = 0; % ATTENTION! If this is 1, make sure fc is aligned with bin file!!!

num_try = 10; % how many times we try for each frequency or file
num_radioframe = 8; % each radio frame length 10ms. MIB period is 4 radio frame
sampling_rate = 1.92e6; % LTE spec. 30.72MHz/16.
num_subframe_per_radioframe = 10;
len_time_subframe = 1e-3; % 1ms. LTE spec
num_sample_per_radioframe = num_subframe_per_radioframe*len_time_subframe*sampling_rate;
num_sample = num_radioframe*num_sample_per_radioframe;

coef_10M = fir1(158, ((9e6+400e3))/15.36e6); %freqz(coef, 1, 1024); %10M channel filter. 
coef_pbch = fir1(158, (0.18e6*6+150e3)/15.36e6); %freqz(coef, 1, 1024);

DS_COMB_ARM = 2;
FS_LTE = 30720000;
thresh1_n_nines=12;
rx_cutoff=(6*12*15e3/2+4*15e3)/(FS_LTE/16/2);
THRESH2_N_SIGMA = 3;

[~, td_pss] = pss_gen;

f_search_set = 20e3:5e3:30e3; % change it wider if you don't know pre-information
pss_fo_set = pss_fo_set_gen(td_pss, f_search_set);

loop_size = 1;
freq_idx = 1;

peaks_store = cell(1,loop_size);
detect_flag_store = cell(1,loop_size);
tdd_flags_store = cell(1,loop_size);

r_all = get_signal_from_bin(filename, inf);

r_all = filter_wo_tail(r_all, coef_10M, 1);
r_all = filter_wo_tail(r_all, coef_pbch, 15.36/1.92);

plot(real(r_all)); drawnow;

for try_idx = 1 : num_try
    disp(['Try idx ' num2str(try_idx)]);

    sp = (try_idx-1)*num_sample + 1;
    ep = sp + num_sample - 1;
    r = r_all(sp:ep);

    disp(['Input averaged abs: ' num2str( mean(abs([real(r); imag(r)])) )]);
    disp(['Processing  at ' filename]);

    capbuf = r.';

    disp('sampling_ppm_f_search_set_by_pss: try ... ... ');
    [period_ppm, dynamic_f_search_set, xc] = sampling_ppm_f_search_set_by_pss(r, f_search_set, pss_fo_set, sampling_carrier_twist);
    if sampling_carrier_twist==0
        if period_ppm == inf
            disp('No valid PSS is found at pre-proc phase! Please try again.');
            peaks = [];
            detect_flag = [];
            tdd_flags = [];
            continue;
        else
            k_factor_set=(1+period_ppm.*1e-6);
        end

        peaks = [];
        for i=1:length(k_factor_set)
            col_idx = i:length(k_factor_set):3*length(k_factor_set);

            [xc_incoherent_collapsed_pow, xc_incoherent_collapsed_frq, n_comb_xc, n_comb_sp, xc_incoherent_single, xc_incoherent, sp_incoherent, sp]= ...
            xcorr_pss(capbuf,dynamic_f_search_set(i),DS_COMB_ARM,fc,sampling_carrier_twist,k_factor_set(i), xc(:,:,i));

            R_th1=chi2inv(1-(10.0^(-thresh1_n_nines)), 2*n_comb_xc*(2*DS_COMB_ARM+1));
            Z_th1=R_th1*sp_incoherent/rx_cutoff/137/2/n_comb_xc/(2*DS_COMB_ARM+1);

            tmp_peak=peak_search(xc_incoherent_collapsed_pow,xc_incoherent_collapsed_frq,Z_th1,dynamic_f_search_set(i),fc,sampling_carrier_twist,k_factor_set(i));
            peaks = [peaks tmp_peak];
        end
    else
        [xc_incoherent_collapsed_pow, xc_incoherent_collapsed_frq, n_comb_xc, n_comb_sp, xc_incoherent_single, xc_incoherent, sp_incoherent, sp]= ...
        xcorr_pss(capbuf,dynamic_f_search_set,DS_COMB_ARM,fc,sampling_carrier_twist,NaN, xc);

        R_th1=chi2inv(1-(10.0^(-thresh1_n_nines)), 2*n_comb_xc*(2*DS_COMB_ARM+1));
        Z_th1=R_th1*sp_incoherent/rx_cutoff/137/2/n_comb_xc/(2*DS_COMB_ARM+1);

        peaks=peak_search(xc_incoherent_collapsed_pow,xc_incoherent_collapsed_frq,Z_th1,dynamic_f_search_set,fc, sampling_carrier_twist,NaN);
    end

    detect_flag = zeros(1, length(peaks));
    tdd_flags = zeros(1, length(peaks));
    for i=1:length(peaks)
        for tdd_flag=0:1
            peak = sss_detect(peaks(i),capbuf,THRESH2_N_SIGMA,fc,sampling_carrier_twist,tdd_flag);
            if ~isnan( peak.n_id_1 )
                break;
            end
        end
        if isnan( peak.n_id_1 )
            continue;
        end
        peak=pss_sss_foe(peak,capbuf,fc,sampling_carrier_twist,tdd_flag);
        [tfg, tfg_timestamp]=extract_tfg(peak,capbuf,fc,sampling_carrier_twist);
        [tfg_comp, tfg_comp_timestamp, peak]=tfoec(peak,tfg,tfg_timestamp,fc,sampling_carrier_twist);
        peak=decode_mib(peak,tfg_comp);
        if isnan( peak.n_rb_dl)
            continue;
        end
        if tdd_flag == 1
            disp('  Detected a TDD cell!');
        else
            disp('  Detected a FDD cell!');
        end
        disp(['  at ' filename]);
        disp(['    cell ID: ' num2str(peak.n_id_cell)]);
        disp(['    PSS  ID: ' num2str(peak.n_id_2+1)]);
        disp(['    RX power level: ' num2str(10*log10(peak.pow))]);
        disp(['    residual frequency offset: ' num2str(peak.freq_superfine)]);
        peaks(i) = peak;
        detect_flag(i) = 1;
        tdd_flags(i) = tdd_flag;
    end
    peaks_store{freq_idx} = peaks;
    detect_flag_store{freq_idx} = detect_flag;
    tdd_flags_store{freq_idx} = tdd_flags;
    if sum(detect_flag)==0
        disp('No LTE cells were found...');
    else
        break;
    end
end


% show all Cells information
disp(' ');
disp('-------------------------------Cells information summary-------------------------------');

disp(['At ' filename]);
peaks = peaks_store{freq_idx};
detect_flag = detect_flag_store{freq_idx};
tdd_flags = tdd_flags_store{freq_idx};
if isempty(detect_flag)
    disp('No valid PSS is found at pre-proc phase! Please try again.');
else
    if sum(detect_flag)
        hit_idx = find(detect_flag);
        for i=1:length(hit_idx);
            peak = peaks(hit_idx(i));
            tdd_flag = tdd_flags(hit_idx(i));
            if tdd_flag == 1
                cell_mode_str = 'TDD';
            else
                cell_mode_str = 'FDD';
            end
            disp(['Cell ' num2str(i) ' information:--------------------------------------------------------']);
            disp(['            Cell mode: ' num2str(cell_mode_str)]);
            disp(['              Cell ID: ' num2str(peak.n_id_cell)]);
            disp(['   Num. eNB Ant ports: ' num2str(peak.n_ports)]);
            disp(['    Carrier frequency: ' num2str(fc/1e6) 'MHz']);
            disp(['Residual freq. offset: ' num2str(peak.freq_superfine/1e3) 'kHz']);
            disp(['       RX power level: ' num2str(10*log10(peak.pow))]);
            disp(['              CP type: ' peak.cp_type]);
            disp(['              Num. RB: ' num2str(peak.n_rb_dl)]);
            disp(['       PHICH duration: ' peak.phich_dur]);
            disp(['  PHICH resource type: ' num2str(peak.phich_res)]);
        end
    else
        disp('No LTE cells were found...  Please try again.');
    end
end
