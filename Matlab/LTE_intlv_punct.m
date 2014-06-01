% the caller use: punct_pattern = LTE_intlv_punct(1:3.*L_total, currentCBS);
function e = LTE_intlv_punct(d, E)
d = vec2mat(d, 3);
d0 = d(:,1).';
d1 = d(:,2).';
d2 = d(:,3).';
[v0, R] = sub_block_intlv(d0, 0);
v1 = sub_block_intlv(d1, 1);
v2 = sub_block_intlv(d2, 2);

e = bit_C_S_T(v0, v1, v2, R, E);
