function peak = pss_sss_detect(s, sampling_rate)

pbch_sampling_rate = 1.92e6;

downsampling_ratio = sampling_rate/pbch_sampling_rate;

% LTE PSS/SSS/PBCH/ channle filter. 10MHz to 1.08MHz; 15.36Msps to 1.92Msps
coef = fir1(158, (0.18e6*6+150e3)/sampling_rate); %freqz(coef, 1, 1024);

num_radioframe = 8; % each radio frame length 10ms. MIB period is 4 radio frame

num_sample = num_radioframe*10e-3*sampling_rate;

sampling_carrier_twist = 0;
f_search_set = -50e3:5e3:45e3;
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
%         break;
    end
    radioframe_idx = radioframe_idx + 1;
    sp = (radioframe_idx-1)*num_sample + 1;
    ep = sp + num_sample - 1;
end

peak = inf;
