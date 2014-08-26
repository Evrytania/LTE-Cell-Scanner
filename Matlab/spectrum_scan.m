function spectrum_scan(freq_start, freq_end, lna_gain, vga_gain)
close all;
[~, lna_gain_new, vga_gain_new] = hackrf_gain_regulation(0, lna_gain, vga_gain);

freq_start = freq_start*1e6;
freq_end = freq_end*1e6;

bw = 10e6;

num_bw = ceil( (freq_end - freq_start)/bw );
freq_end = freq_start + bw*num_bw;

disp(['Actual range ' num2str(freq_start/1e6) 'MHz to ' num2str(freq_end/1e6) 'MHz']);
freq_set = (freq_start+(bw/2)):bw:(freq_end-(bw/2));
disp(['freq set (MHz) ' num2str(freq_set./1e6)]);

filename_raw = 'hackrf_live_tmp.bin';

len_time = 500e-3;
len_time_pre = 5e-3;
len_time_total = len_time + len_time_pre;

sampling_rate = 15e6;
half_num_sample = (len_time*bw)/2;

rbw = 40e3;
spec = zeros(1, length(freq_set)*half_num_sample*2);

clf;
for i = 1 : length(freq_set)
    freq = freq_set(i);
    delete(filename_raw);
    cmd_str = ['hackrf_transfer -r ' filename_raw ' -f ' num2str(freq) ' -s ' num2str(sampling_rate) ' -b ' num2str(sampling_rate) ' -n ' num2str(len_time_total*sampling_rate) ' -l ' num2str(lna_gain_new) ' -g ' num2str(vga_gain_new) ];
    system(cmd_str);
    
    fid_raw = fopen(filename_raw, 'r');
    if fid_raw == -1
        disp('Open hackrf_live_tmp.bin failed!');
        return;
    end
    a = fread(fid_raw, inf, 'int8');
    fclose(fid_raw);
    
    b = a((2*(len_time_pre*sampling_rate)+1) : end).';
    c = b(1:2:end) + 1i.*b(1:2:end);
    c = c - mean(c);
    
    figure(1);
    sp = (i-1)*len_time*sampling_rate + 1;
    ep = sp + len_time*sampling_rate - 1;
    plot(sp:ep, real(c)); hold on; drawnow;
    
    d = [ c( (end-half_num_sample+1) : end), c(1 : half_num_sample) ];
    sp = (i-1)*half_num_sample*2 + 1;
    ep = sp + half_num_sample*2 - 1;
    spec(sp:ep) = abs(fft(d)).^2;
end

spec = filter_wo_tail(spec.', ones(1, floor(rbw*len_time)+1), 1);

figure(2);
plot( ( freq_start + (0:(length(freq_set)*half_num_sample*2-1))./len_time )./1e6, 10.*log10(spec)); hold on;

