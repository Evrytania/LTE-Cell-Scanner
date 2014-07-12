function source_info = regression_test_source(source_dir)
tmp_info = struct('filename', -1, 'fc', -1, 'fs', -1, 'dev', -1);

dir_info = dir([source_dir '/f*_s*.bin']);

num_file = length(dir_info);
if num_file == 0
    return;
end

source_info(1:num_file) = tmp_info;

for i = 1 : num_file
    filename = [source_dir, '/', dir_info(i).name];
%     disp(filename);
    
    fc = get_frequency_carrier(dir_info(i).name);
    fs = get_sampling_rate(dir_info(i).name);
    dev = get_dev_name(dir_info(i).name);
    
    source_info(i).filename = filename;
    source_info(i).fc = fc;
    source_info(i).fs = fs;
    source_info(i).dev = dev;
end

function dev = get_dev_name(filename)
dev = -1;
tmp1 = strfind(filename, 'rtlsdr');
tmp2 = strfind(filename, 'hackrf');
tmp3 = strfind(filename, 'usrp');
tmp4 = strfind(filename, 'bladerf');

tmp1 = ~isempty(tmp1);
tmp2 = ~isempty(tmp2);
tmp3 = ~isempty(tmp3);
tmp4 = ~isempty(tmp4);

tmp_sum = tmp1+tmp2+tmp3+tmp4;
if tmp_sum > 1 || tmp_sum == 0
    return;
end

if tmp1
    dev = 'rtlsdr';
elseif tmp2
    dev = 'hackrf';
elseif tmp3
    dev = 'usrp';
elseif tmp4
    dev = 'bladerf';
end

function fc = get_frequency_carrier(filename)
fc = -1;
sp = 2;
ep = strfind(filename, '_');
if isempty(ep)
    return;
end

ep = ep(1) - 1;
fc = str2double(filename(sp:ep))*1e6;

function sampling_rate = get_sampling_rate(filename)
sampling_rate = -1;

sp_set = strfind(filename, '_s');

for i = 1 : length(sp_set)
    number_sp = sp_set(i)+2;
    number_ep = strfind(filename(number_sp:end), '_');
    number_ep = number_sp + number_ep - 1;
    
    if ~isempty(number_ep)
        number_ep = number_ep(1) - 1;
    else
        disp('Abnormal when parse sampling_rate');
        return;
    end
    
    sampling_rate = str2double(filename(number_sp:number_ep))*1e6;
    if ~isnan(sampling_rate)
        break;
    end
end

if isnan(sampling_rate)
    sampling_rate = -1;
    disp('sampling_rate is NaN!');
end

