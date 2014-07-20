function [corr_store, fo_search_set_new] = fft_corr(s, td_pss, fo_search_set)
% % s column vectors -- 153600
% % td_pss column vectors
% % fo_search_set row vectors

fo_search_set_new = fo_search_set; % it happens to be equal, because 1.92e6/153600 = 12.5Hz under 5kHz step size
sampling_rate = 1.92e6;

len = length(s);

freq_step = sampling_rate/len;

len_pss = size(td_pss, 1);
num_pss = size(td_pss, 2);
num_fo = length(fo_search_set);
num_fo_pss = num_fo*num_pss;
% len_short = len - (len_pss-1);

corr_store = zeros(len, num_fo_pss);

s = fft(s);
fd_pss = fft(conj(td_pss(end:-1:1,:)./len_pss), len, 1);

for i = 1 : num_fo_pss
    pss_idx = floor((i-1)/num_fo) + 1;
    fo_idx = i - (pss_idx-1)*num_fo;
    fo = fo_search_set(fo_idx);
    fd_fo_shift_len = fo/freq_step;
    corr_store(:, i) = s.*circshift(fd_pss(:,pss_idx), [fd_fo_shift_len,0]);
%     corr_store(:, i) = abs( ifft(s.*circshift(fd_pss(:,pss_idx), [fd_fo_shift_len,0]), [], 1) ).^2;
%     corr_store(:, i) = real( corr_store(:, i).*conj(corr_store(:, i)) );
end
corr_store = abs( ifft(corr_store, len, 1) ).^2;
corr_store = corr_store(len_pss:end,:);
