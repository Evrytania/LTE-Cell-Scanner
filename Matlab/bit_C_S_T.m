function e = bit_C_S_T(v0, v1, v2, R, E)
Kpi = length(v0);
Kw = 3.*Kpi;
w = zeros(1, Kw);
w(1:Kpi) = v0;
w(Kpi+1:2:end) = v1;
w(Kpi+2:2:end) = v2;

Ncb = Kw;

k0 = 2.*R;

k = 0; j = 0; e = zeros(1, E);
while k<E
    temp = w( mod(k0+j, Ncb)+1 );
    if  temp~= -63534
        e(k+1) = temp;
        k = k + 1;
    end
    j = j + 1;
end
