function show_signal_time_frequency(s, sampling_rate, rbw)
clf;
s = s(:).';

subplot(2,1,1); plot(1e6.*(0:(length(s)-1))./sampling_rate, real(s));hold on;
subplot(2,1,1); plot(1e6.*(0:(length(s)-1))./sampling_rate, imag(s), 'r'); grid on;

spec = abs(fft(s)).^2;
spec = [spec((end/2)+1 : end) spec(1:(end/2))];

len_time = length(s)/sampling_rate;
spec = filter_wo_tail(spec.', ones(1, floor(rbw*len_time)+1), 1);

subplot(2,1,2); plot( (1e-6.*(0:(length(s)-1))./len_time ) - (1e-6*sampling_rate/2), 10.*log10(spec)); grid on;


