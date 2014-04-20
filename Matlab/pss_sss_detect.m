function [peak, s] = pss_sss_detect(s, sampling_rate, f_search_set)

pbch_sampling_rate = 1.92e6;

downsampling_ratio = sampling_rate/pbch_sampling_rate;

% LTE PSS/SSS/PBCH/ channle filter. 10MHz to 1.08MHz; 15.36Msps to 1.92Msps
coef = fir1(158, (0.18e6*6+150e3)/sampling_rate); %freqz(coef, 1, 1024);

num_radioframe = 8; % each radio frame length 10ms. MIB period is 4 radio frame

num_sample = num_radioframe*10e-3*sampling_rate;

fc = inf; % fake because sampling_carrier_twist = 0
sampling_carrier_twist = 0;
% f_search_set = -50e3:5e3:45e3;
[~, td_pss] = pss_gen;
pss_fo_set = pss_fo_set_gen(td_pss, f_search_set);

radioframe_idx = 1;
sp = (radioframe_idx-1)*num_sample + 1;
ep = sp + num_sample - 1;
while ep<length(s)
    disp(' ');
    disp(['pss_sss_detect: ------------------radioframe_idx = ' num2str(radioframe_idx)]);
    
    r = filter_wo_tail(s(sp:ep), coef, downsampling_ratio);
    
    [period_ppm, dynamic_f_search_set, xc] = sampling_ppm_f_search_set_by_pss(r, f_search_set, pss_fo_set, sampling_carrier_twist);
    if period_ppm ~= inf
        disp(['pss_sss_detect: pss detected! radioframe_idx = ' num2str(radioframe_idx)]);
        break;
    end
    radioframe_idx = radioframe_idx + 1;
    sp = (radioframe_idx-1)*num_sample + 1;
    ep = sp + num_sample - 1;
end

if period_ppm ~= inf
    s = s(sp:end).';
    k_factor_set=(1+period_ppm.*1e-6);
else
    peak = inf;
    return;
end

DS_COMB_ARM = 2;
FS_LTE = 30720000;
thresh1_n_nines=12;
rx_cutoff=(6*12*15e3/2+4*15e3)/(FS_LTE/16/2);
THRESH2_N_SIGMA = 3;

[xc_incoherent_collapsed_pow, xc_incoherent_collapsed_frq, n_comb_xc, n_comb_sp, xc_incoherent_single, xc_incoherent, sp_incoherent, sp]= ...
xcorr_pss(s,dynamic_f_search_set,DS_COMB_ARM,fc,sampling_carrier_twist,k_factor_set, xc(:,:,1));

R_th1=chi2inv(1-(10.0^(-thresh1_n_nines)), 2*n_comb_xc*(2*DS_COMB_ARM+1));
Z_th1=R_th1*sp_incoherent/rx_cutoff/137/2/n_comb_xc/(2*DS_COMB_ARM+1);

peaks=peak_search(xc_incoherent_collapsed_pow,xc_incoherent_collapsed_frq,Z_th1,dynamic_f_search_set,fc,sampling_carrier_twist,k_factor_set);

peaks_tmp = [];
for i=1:length(peaks)
    for tdd_flag=0:1
        peak = sss_detect(peaks(i),s,THRESH2_N_SIGMA,fc,sampling_carrier_twist,tdd_flag);
        if ~isnan( peak.n_id_1 )
            break;
        end
    end
    if isnan( peak.n_id_1 )
        continue;
    else
        peaks_tmp = [peaks_tmp, peak];
        disp(['pss_sss_detect: sss detected! i = ' num2str(i) '; tdd_flag = ' num2str(tdd_flag)]);
    end
end

if length(peaks_tmp) == 0
    disp('pss_sss_detect: No sss detected.');
    peak = inf;
    return;
end

peak = peaks_tmp;
