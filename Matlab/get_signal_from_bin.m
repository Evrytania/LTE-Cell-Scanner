function s = get_signal_from_bin(filename, num_sample_read)

% filename = '../test/f1860_s15.36_bw10_1s.bin'; % FDD 20MHZ
% filename = '../test/f1890_s15.36_bw10_1s.bin'; % TDD 20MHz
fid = fopen(filename);

if fid == -1
    disp('get_signal_from_hackrf_bin: Can not open file!');
    return;
end

[s, count] = fread(fid, num_sample_read*2, 'int8');
fclose(fid);

if num_sample_read~=inf && count ~= (num_sample_read*2)
    s = -1;
    clear s;
    disp('get_signal_from_hackrf_bin: No enough samples in the file!');
    return;
end

s = s(1:2:end) + 1i.*s(2:2:end);

% len_s = length(s);
% 
% s = s((len_s/2)+1:end);