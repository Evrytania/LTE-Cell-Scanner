% Jiao Xianjun (putaoshu@msn.com; putaoshu@gmail.com)
% Extract part of rtl-sdr captured bin into a new file.
% A script of project: https://github.com/JiaoXianjun/rtl-sdr-LTE

function extract_part_from_bin_file(input_file_name, num_skip_bytes, num_extract_bytes, output_file_name)

dev = get_dev_name_from_filename(input_file_name);

fid = fopen(input_file_name);

if fid==-1
    disp('Can not open input file!');
    return;
end

if strcmpi(dev, 'rtlsdr')
    s = fread(fid, inf, 'uint8');
elseif strcmpi(dev, 'hackrf')
    s = fread(fid, inf, 'int8');
else
    disp('Not supported dev/file keyword/format!');
    return;
end

fclose(fid);

if (num_skip_bytes+num_extract_bytes)>length(s)
    disp('num_skip_bytes + num_extract_bytes > len_input');
    return;
end

s = s((num_skip_bytes+1) : (num_skip_bytes+num_extract_bytes));

fid = fopen(output_file_name, 'w');

if fid==-1
    disp('Can not open output file!');
    return;
end

if strcmpi(dev, 'rtlsdr')
    count = fwrite(fid, s, 'uint8');
elseif strcmpi(dev, 'hackrf')
    count = fwrite(fid, s, 'int8');
else
    disp('Not supported dev/file keyword/format!');
    return;
end

fclose(fid);

if count ~= length(s)
    disp('Number of written bytes is not as expectation!');
end
