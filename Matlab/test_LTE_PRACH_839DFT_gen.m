% function test_LTE_PRACH_839DFT_gen
clear all; close all;

n = 0:838;
u = 719;
zc = exp(-1i.*pi.*u.*n.*(n+1)./839);
fft_zc = fft(zc);
fft_zc = exp(1i.*angle(fft_zc));

n = 4096;
fft_zc1 = fft(zc, n);
fft_zc1 = interp1(0 : (n-1), fft_zc1, (0:838).*n./839, 'linear', 'extrap');
fft_zc1 = exp(1i.*angle(fft_zc1));

figure(1);
plot(angle(fft_zc)); hold on;
plot(angle(fft_zc1), 'r');

figure(3);
plot(abs(fft_zc-fft_zc1).^2);

10.*log10(sum(abs(fft_zc-fft_zc1).^2)/839)
