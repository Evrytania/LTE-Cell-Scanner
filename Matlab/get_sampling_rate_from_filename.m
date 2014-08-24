
function sampling_rate = get_sampling_rate_from_filename(filename)
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

