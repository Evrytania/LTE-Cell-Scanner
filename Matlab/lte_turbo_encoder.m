function bits_encoded = lte_turbo_encoder(bits)
code_length = length(bits);

intlv_table = intlvLTE(code_length);

bits1 = bits(intlv_table);

reg0 = zeros(1, 3);
reg1 = zeros(1, 3);

bits_encoded = zeros(1, 3*code_length+12);
j = 1;
for i = 1 : code_length
    bits_encoded(j) = bits(i);
    
    tmp_val = mod( bits(i) + reg0(2) + reg0(3), 2 );
    j = j + 1;
    bits_encoded(j) = mod( tmp_val + reg0(1) + reg0(3), 2 );
    reg0(2:3) = reg0(1:2);
    reg0(1) = tmp_val;
    
    tmp_val = mod( bits1(i) + reg1(2) + reg1(3), 2 );
    j = j + 1;
    bits_encoded(j) = mod( tmp_val + reg1(1) + reg1(3), 2 );
    reg1(2:3) = reg1(1:2);
    reg1(1) = tmp_val;
    
    j = j + 1;
end

bits_encoded(j) = xor(reg0(2), reg0(3));
j = j + 1;
bits_encoded(j) = xor(reg0(1), reg0(3));
reg0(2:3) = reg0(1:2);
reg0(1) = 0;

j = j + 1;
bits_encoded(j) = xor(reg0(2), reg0(3));
j = j + 1;
bits_encoded(j) = xor(reg0(1), reg0(3));
reg0(2:3) = reg0(1:2);
reg0(1) = 0;

j = j + 1;
bits_encoded(j) = xor(reg0(2), reg0(3));
j = j + 1;
bits_encoded(j) = xor(reg0(1), reg0(3));

j = j + 1;
bits_encoded(j) = xor(reg1(2), reg1(3));
j = j + 1;
bits_encoded(j) = xor(reg1(1), reg1(3));
reg1(2:3) = reg1(1:2);
reg1(1) = 0;

j = j + 1;
bits_encoded(j) = xor(reg1(2), reg1(3));
j = j + 1;
bits_encoded(j) = xor(reg1(1), reg1(3));
reg1(2:3) = reg1(1:2);
reg1(1) = 0;

j = j + 1;
bits_encoded(j) = xor(reg1(2), reg1(3));
j = j + 1;
bits_encoded(j) = xor(reg1(1), reg1(3));



