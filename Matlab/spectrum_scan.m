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

len_time = 160e-3;
len_time_pre = 10e-3;
len_time_total = len_time + len_time_pre;

sampling_rate = 15e6;
half_num_sample = (((20e-3)*sampling_rate)*2/3)/2;

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
    
    b = a((2*(len_time_pre*sampling_rate)+1) : end);
    c = b(1:2:end) + 1i.*b(1:2:end);
    c = c - mean(c);
    d = vec2mat(c, (20e-3)*sampling_rate);
    
    sub_spec = zeros(1, (20e-3)*bw);
    for j = 1 : 8
        tmp_spec = 10.*log10(abs(fft(d(j,:))).^2);
        tmp_spec(1) = inf;
        tmp_spec = [tmp_spec((end-half_num_sample+1) : end) tmp_spec(1:half_num_sample)];
        sub_spec = sub_spec + tmp_spec;
    end
    sub_spec = sub_spec./8;
%     spec = pwelch(c);
    
    plot(linspace((freq-(bw/2))./1e6, (freq+(bw/2))./1e6, length(sub_spec)), sub_spec); hold on; drawnow;
end

