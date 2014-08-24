
function fc = get_frequency_carrier_from_filename(filename)
fc = -1;
sp = 2;
ep = strfind(filename, '_');
if isempty(ep)
    return;
end

ep = ep(1) - 1;
fc = str2double(filename(sp:ep))*1e6;
