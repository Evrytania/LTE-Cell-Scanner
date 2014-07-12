function s = get_signal_from_bin(filename, num_sample_read)

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

s = (s(1:2:end) + 1i.*s(2:2:end))./128;
