% function test_fft_corr
clear all;
close all;
[~, td_pss] = pss_gen;
td_pss = td_pss(:,1).';
len_pss = length(td_pss);

fo_search_set = 5e3;

sampling_rate = 1.92e6;
td_pss_fo = td_pss.*exp(1i.*2.*pi.*(1./sampling_rate).*(0:(len_pss-1)).*fo_search_set);

conj_td_pss_fo = conj(td_pss_fo);

len = 153600;
s = randn(1, len) + 1i.*randn(1, len);

freq_step = sampling_rate/len;

shift_len = fo_search_set/freq_step;

fft_conj_td_pss_fo = fft(conj_td_pss_fo, len);

fft_conj_td_pss = fft(conj(td_pss), len);

plot(angle(fft_conj_td_pss_fo)); hold on;
plot(angle(circshift(fft_conj_td_pss, [0, -shift_len])), 'r');

fft_s = fft(s);

% a = filter(conj_td_pss_fo(end:-1:1), 1, s);
% a = conv(s, conj_td_pss_fo);
% a(1: len_pss-1) = a(1: len_pss-1) + a(len+1 : end);
a = filter(conj_td_pss_fo, 1, s);

b = ifft( fft_s.*circshift(fft_conj_td_pss, [0, -shift_len]) );
% b = ifft( fft_s.*fft_conj_td_pss_fo );

% b = b(len_pss:end);

figure;
plot(angle(a)); hold on;
plot(angle(b),'r');
