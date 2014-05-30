a = load('f2360_s19.2_bw20_1s_hackrf_bda_SIB.txt');

num_hex = length(a)/4;

hex_out = [];
for i = 1 : num_hex
    sp = (i-1)*4 + 1;
    ep = sp + 3;
    tmp_bin = bi2de(a(sp:ep), 'left-msb');
    hex_out = [hex_out dec2hex(tmp_bin)];
end
disp(hex_out)
