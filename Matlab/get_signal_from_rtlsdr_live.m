function s = get_signal_from_rtlsdr_live(freq, sampling_rate, num_second, gain)

bin_filename = 'rtlsdr_live_tmp.bin';

cmd_str = ['rtl_sdr -f ' num2str(freq) ' -s ' num2str(sampling_rate) ' -n ' num2str(num_second*sampling_rate) ' -g ' num2str(gain) ' '  bin_filename];

system(cmd_str);

s = get_signal_from_bin(bin_filename, inf, 'rtlsdr');
plot(abs(s));

s = s - mean(s);
% subplot(2,1,1); plot(real(s));
% subplot(2,1,2); plot(imag(s));


% figure;
% plot(angle(s));
