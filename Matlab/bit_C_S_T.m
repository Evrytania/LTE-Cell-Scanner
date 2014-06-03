function e = bit_C_S_T(v0, v1, v2, R, E, max_num_HARQ, rv)
Kpi = length(v0);
Kw = 3.*Kpi;
w = zeros(1, Kw);
w(1:Kpi) = v0;
w(Kpi+1:2:end) = v1;
w(Kpi+2:2:end) = v2;

C = 1; % attention. need to change according to Code block segmentation and code block CRC attachment

N_soft = 250368; % cat 1
% N_soft = 1237248; % cat 2
% N_soft = 1237248; % cat 3
% N_soft = 1827072; % cat 4
% N_soft = 3667200; % cat 5
% N_soft = 3667200; % cat 5
% N_soft = 3654144; % cat 6
% N_soft = 3654144; % cat 7
% N_soft = 35982720; % cat 8

if N_soft == 35982720
    Kc = 5;
elseif  N_soft == 3654144
    Kc = 2;
else
    Kc = 1;
end

K_MIMO = 1;

M_limit = 8;

M_DL_HARQ = max_num_HARQ;

N_IR = floor( N_soft/( Kc*K_MIMO*min(M_DL_HARQ, M_limit) ) );

Ncb = min(floor(N_IR/C), Kw);

k0 = R*( 2 + 2*rv*ceil(Ncb/(8*R)) );

k = 0; j = 0; e = zeros(1, E);
while k<E
    temp = w( mod(k0+j, Ncb)+1 );
    if  temp~= -63534
        e(k+1) = temp;
        k = k + 1;
    end
    j = j + 1;
end
