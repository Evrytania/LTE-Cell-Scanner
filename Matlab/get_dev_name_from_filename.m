function dev = get_dev_name_from_filename(filename)
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
