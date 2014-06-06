function hex_out = convert_SIB_bin_to_HEX(a)
% a = load('f2360_s19.2_bw20_1s_hackrf_bda_SIB.txt');

num_hex = length(a)/4;

hex_out = zeros(1, num_hex);
for i = 1 : num_hex
    sp = (i-1)*4 + 1;
    ep = sp + 3;
    tmp_bin = bi2de(a(sp:ep), 'left-msb');
    hex_out(i) = dec2hex(tmp_bin);
end
% disp(char(hex_out))
