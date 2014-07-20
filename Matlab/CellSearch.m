% Jiao Xianjun (putaoshu@gmail.com; putaoshu@msn.com)
% CellSearch.m
% Improved LTE-Cell-Scanner (written by James Peroulas: https://github.com/Evrytania/LTE-Cell-Scanner).
% See also README in root directory, ../test, ../../rtl-sdr-LTE/scan-capture/.

% Some scripts are borrowed from:
% https://github.com/JiaoXianjun/rtl-sdr-LTE
% https://github.com/Evrytania/LTE-Cell-Scanner
% https://github.com/Evrytania/Matlab-Library
% https://github.com/JiaoXianjun/multi-rtl-sdr-calibration

function [cell_info, r_pbch_sub, r_20M_sub, cell_info_return] = CellSearch(r_pbch, r_20M, f_search_set, fc, sampling_carrier_twist, pss_peak_max_reserve, num_pss_period_try, combined_pss_peak_range, par_th, num_peak_th)
r_20M_sub = -1;

[~, td_pss] = pss_gen;

% f_search_set = 20e3:5e3:30e3; % change it wider if you don't know pre-information
pss_fo_set = pss_fo_set_gen(td_pss, f_search_set);

% sampling_carrier_twist = 0; % ATTENTION! If this is 1, make sure fc is aligned with bin file!!!

num_radioframe = 8; % each radio frame length 10ms. MIB period is 4 radio frame

pbch_sampling_ratio = 16;
sampling_rate = 30.72e6;
sampling_rate_pbch = sampling_rate/pbch_sampling_ratio; % LTE spec. 30.72MHz/16.

num_subframe_per_radioframe = 10;
len_time_subframe = 1e-3; % 1ms. LTE spec
num_sample_per_radioframe = num_subframe_per_radioframe*len_time_subframe*sampling_rate_pbch;
num_sample_pbch = num_radioframe*num_sample_per_radioframe;

DS_COMB_ARM = 2;
FS_LTE = 30720000;
thresh1_n_nines=12;
rx_cutoff=(6*12*15e3/2+4*15e3)/(FS_LTE/16/2);
THRESH2_N_SIGMA = 3;
ex_gain = 2;

num_try = floor(length(r_pbch)/num_sample_pbch);
tdd_fdd_str = {'TDD', 'FDD'};
for try_idx = 1 : num_try
    disp(['Try idx ' num2str(try_idx)]);

    sp = (try_idx-1)*num_sample_pbch + 1;
    ep = sp + num_sample_pbch - 1;
    capbuf_pbch = r_pbch(sp:ep).';
    
    sp_20M = (sp-1)*pbch_sampling_ratio + 1;
    ep_20M = ep*pbch_sampling_ratio;
    r_pbch_sub = capbuf_pbch;
    
    if ~isempty(r_20M)
        r_20M_sub = r_20M(sp_20M:ep_20M);
    end

    disp(['Input averaged abs: ' num2str( mean(abs([real(capbuf_pbch) imag(capbuf_pbch)])) )]);

    disp('sampling_ppm_f_search_set_by_pss: try ... ... ');
    [period_ppm, dynamic_f_search_set, xc, ~, ~, ~, ~, extra_info] = sampling_ppm_f_search_set_by_pss(capbuf_pbch.', f_search_set, pss_fo_set, sampling_carrier_twist, pss_peak_max_reserve, num_pss_period_try, combined_pss_peak_range, par_th, num_peak_th);
    if sampling_carrier_twist==0
        if period_ppm == inf
            disp('No valid PSS is found at pre-proc phase! Please try again.');
            peaks = [];
            detect_flag = [];
            continue;
        else
            k_factor_set=(1+period_ppm.*1e-6);
        end

        peaks = [];
        for i=1:length(k_factor_set)
            [xc_incoherent_collapsed_pow, xc_incoherent_collapsed_frq, n_comb_xc, ~, ~, ~, sp_incoherent, ~]= ...
            xcorr_pss(capbuf_pbch,dynamic_f_search_set(i),DS_COMB_ARM,fc,sampling_carrier_twist,k_factor_set(i), xc(:,:,i));

            R_th1=chi2inv(1-(10.0^(-thresh1_n_nines)), 2*n_comb_xc*(2*DS_COMB_ARM+1));
            Z_th1=ex_gain.*R_th1*sp_incoherent/rx_cutoff/137/2/n_comb_xc/(2*DS_COMB_ARM+1);

            tmp_peak=peak_search(xc_incoherent_collapsed_pow,xc_incoherent_collapsed_frq,Z_th1,dynamic_f_search_set(i),fc,sampling_carrier_twist,k_factor_set(i));
            for j = 1 : length(tmp_peak)
                extra_info.try_idx = try_idx;
                tmp_peak(j).extra_info = extra_info;
            end
            peaks = [peaks tmp_peak];
        end
    else
        [xc_incoherent_collapsed_pow, xc_incoherent_collapsed_frq, n_comb_xc, ~, ~, ~, sp_incoherent, sp]= ...
        xcorr_pss(capbuf_pbch,dynamic_f_search_set,DS_COMB_ARM,fc,sampling_carrier_twist,NaN, xc);

        R_th1=chi2inv(1-(10.0^(-thresh1_n_nines)), 2*n_comb_xc*(2*DS_COMB_ARM+1));
        Z_th1=ex_gain.*R_th1*sp_incoherent/rx_cutoff/137/2/n_comb_xc/(2*DS_COMB_ARM+1);

        peaks=peak_search(xc_incoherent_collapsed_pow,xc_incoherent_collapsed_frq,Z_th1,dynamic_f_search_set,fc, sampling_carrier_twist,NaN);
    end
    
    for j = 1 : length(peaks)
        peaks(j).extra_info.num_peaks_raw = length(peaks)/2;
    end
    disp(['Hit  num peaks ' num2str(length(peaks)/2) ]);
    tdd_flags = kron(ones(1, length(peaks)/2), [0 1]); % even: tdd_flag 0; odd : tdd_flag 1
    
    detect_flag = zeros(1, length(peaks));
    for i=1:length(peaks)
        disp(['try peak ' num2str(floor((i+1)/2)) ' ' tdd_fdd_str{mod(i,2)+1} ' mode']);
        tdd_flag = tdd_flags(i);
        peak = sss_detect(peaks(i),capbuf_pbch,THRESH2_N_SIGMA,fc,sampling_carrier_twist,tdd_flag);
        if ~isnan( peak.n_id_1 )
            peak=pss_sss_foe(peak,capbuf_pbch,fc,sampling_carrier_twist,tdd_flag);
            [tfg, tfg_timestamp]=extract_tfg(peak,capbuf_pbch,fc,sampling_carrier_twist, 6);
            [tfg_comp, ~, peak]=tfoec(peak,tfg,tfg_timestamp,fc,sampling_carrier_twist, 6);
            peak=decode_mib(peak,tfg_comp);
            if isnan( peak.n_rb_dl)
                continue;
            end
            if tdd_flag == 1
                disp('  Detected a TDD cell!');
            else
                disp('  Detected a FDD cell!');
            end
            disp(['    cell ID: ' num2str(peak.n_id_cell)]);
            disp(['    PSS  ID: ' num2str(peak.n_id_2+1)]);
            disp(['    RX power level: ' num2str(10*log10(peak.pow))]);
            disp(['    residual frequency offset: ' num2str(peak.freq_superfine)]);
            peaks(i) = peak;
            detect_flag(i) = 1;
        end
    end
    if sum(detect_flag)==0
        disp('No LTE cells were found...');
    else
        break;
    end
end

cell_info_return = peaks;

cell_info = [];

% show all Cells information
disp(' ');
disp('-------------------------------Cells information summary-------------------------------');

if isempty(detect_flag)
    disp('No valid PSS is found at pre-proc phase! Please try again.');
else
    if sum(detect_flag)
        hit_idx = find(detect_flag);
        for i=1:length(hit_idx);
            peak = peaks(hit_idx(i));
            
            cell_info = [cell_info peak];
            
            if peak.duplex_mode == 1
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
