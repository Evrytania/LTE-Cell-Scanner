function [v, R] = sub_block_intlv(d, i)
D = length(d);
C = 32;
R = ceil(D./C);
Kpi = R.*C;
ND = Kpi - D;
y = [-63534.*ones(1, ND), d];
ymat = vec2mat(y, C);

P = [0, 16, 8, 24, 4, 20, 12, 28, 2, 18, 10, 26, 6, 22, 14, 30, 1, 17, 9, 25, 5, 21, 13, 29, 3, 19, 11, 27, 7, 23, 15, 31];

if i<=1
    ymat = ymat(:,P+1);
    v = ymat(:).';
else
    k = 0 : Kpi-1;
    temp1 = P(floor(k./R) + 1);
    temp2 = C.*mod(k, R);
    pi = mod( temp1 + temp2 +1, Kpi);
    v = y(pi + 1);
end
