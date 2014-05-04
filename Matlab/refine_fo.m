function [freq_new, k_factor_idx] = refine_fo(capbuf, cp_type, n_id_2, freq, fs, frame_start, k_factor_vec)
len_pss = 137;

fo_set = freq + (-3:2:3).*1e3;

len = length(capbuf);

corr_val = zeros(1, 4);
[~, td_pss] = pss_gen;
td_pss = td_pss(:, n_id_2+1).';
for i=1:4
    k_factor_tmp = k_factor_vec(i);
    if strcmpi(cp_type,'extended')
        pss_from_frame_start = k_factor_tmp*(1920+2*(128+32));
    elseif strcmpi(cp_type, 'normal')
        pss_from_frame_start = k_factor_tmp*(1920 + 2*(128+9) + 1);
    end
    
    pss_sp = frame_start + pss_from_frame_start + 3;
    
    pss_fo = conj( td_pss.*exp( 1i.*(0:(len_pss-1)).*fo_set(i).*2.*pi./fs ) );
    
    corr_val(i) = 0;
    pss_count = 0;
    
    while (pss_sp + len_pss + 1) <= len
        
        pss_idx = round(pss_sp);
        chn_tmp = capbuf(pss_idx : (pss_idx + len_pss - 1) );
        corr_tmp = abs(sum(chn_tmp.*pss_fo)).^2;
%         plot(angle(chn_tmp.*pss_fo));
        corr_val(i) = corr_val(i) + corr_tmp;

        chn_tmp = capbuf(pss_idx+1 : (pss_idx+1 + len_pss - 1) );
        corr_tmp = abs(sum(chn_tmp.*pss_fo)).^2;
%         plot(angle(chn_tmp.*pss_fo));
        corr_val(i) = corr_val(i) + corr_tmp;
        
        chn_tmp = capbuf(pss_idx-1 : (pss_idx-1 + len_pss - 1) );
        corr_tmp = abs(sum(chn_tmp.*pss_fo)).^2;
%         plot(angle(chn_tmp.*pss_fo));
        corr_val(i) = corr_val(i) + corr_tmp;
        
        pss_count = pss_count + 1;
        pss_sp = pss_sp + k_factor_tmp*5*1920;
    end
    
    corr_val(i) = corr_val(i)/pss_count;
end

[~, k_factor_idx] = max(corr_val);
freq_new = fo_set(k_factor_idx);

disp(['refine_fo corr_val ' num2str(corr_val)]);
disp(['fo refined from ' num2str(freq) ' to ' num2str(freq_new)]);