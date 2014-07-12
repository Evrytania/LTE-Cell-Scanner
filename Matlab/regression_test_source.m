function [filename, sampling_rate] = regression_test_source(source_dir_idx, varargin)

filename = -1;
sampling_rate = -1;

switch source_dir_idx
    case 1
        source_dir = '~/git/rtl-sdr-LTE/scan-capture/frequency-1850-1880MHz/';
    case 2
        source_dir = '~/git/rtl-sdr-LTE/scan-capture/frequency-1880-1900MHz/';
    case 3
        source_dir = '~/git/rtl-sdr-LTE/scan-capture/frequency-2555-2575MHz/';
    case 4
        source_dir = '~/git/rtl-sdr-LTE/scan-capture/frequency-2575-2595MHz/';
    case 5
        source_dir = '~/git/rtl-sdr-LTE/scan-capture/frequency-2595-2615MHz/';
    case 6
        source_dir = '~/git/rtl-sdr-LTE/scan-capture/frequency-2635-2655MHz/';
    case 7
        source_dir = '~/git/rtl-sdr-LTE/scan-capture/frequency-2635-2655MHz-know-PPM/';
    case 8
        source_dir = '../test/';
    otherwise
        disp('source_dir_idx must be with in [1, 8]!');
        return;
end

% dir_info1 = dir([source_dir '*1s*.bin']);
% dir_info2 = dir([source_dir '*0.8s*.bin']);
% dir_info = [dir_info1; dir_info2];
dir_info = dir([source_dir 'f*.bin']);

num_file = length(dir_info);
if num_file == 0
    filename = -1;
    return;
end

if nargin == 1
    for i = 1 : num_file
        disp([source_dir, dir_info(i).name]);
    end
    filename = num_file;
    return;
end

sub_idx = varargin{1};

if sub_idx > num_file
    disp(['sub_idx exceeds maximum number of files in ' source_dir]);
    return;
end

filename = [source_dir, dir_info(sub_idx).name];
sampling_rate = get_sampling_rate(dir_info(sub_idx).name);

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

