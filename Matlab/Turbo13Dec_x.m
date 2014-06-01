function [d, L_code_out] = Turbo13Dec_x(rec_s, trellis, alpha, niter)
L_total = size(rec_s, 2)/2;

adj = 1e-100;

L_di_store = rec_s(1,1:2:end);
L_e = L_di_store;
L_e = L_e(alpha);
L_a = zeros(1,L_total);

% rec_s(1,1:2:end) = 0;
[pbit0, pbit1] = llr2pb(rec_s(1,:));
pbit0 = ex_odd_even(pbit0); pbit1 = ex_odd_even(pbit1);
pci1 = p_bit2sym(pbit0, pbit1, 2);
% rec_s(1,1:2:end) = 0;
% pcibit1(1,:) = exp(rec_s(1,:))./(1+exp(rec_s(1,:))); pcibit1(2,:) = 1-pcibit1(1,:);
% pci1_removeD=ProbBits2Symbols(pcibit1, 2);

rec_s(2,1:2:end) = 0;
[pbit0, pbit1] = llr2pb(rec_s(2,:));
pbit0 = ex_odd_even(pbit0); pbit1 = ex_odd_even(pbit1);
pci2 = p_bit2sym(pbit0, pbit1, 2);

for iter = 1:niter
    % Decoder one
    L_a(alpha) = L_e;  % a priori info. 

    [pbit0, pbit1] = llr2pb(L_a);
    pdi = [pbit0; pbit1];
    [pco1, pdo] = siso_c_jxjnew(trellis, pci1, pdi);
    
    L_all = pb2llr(pdo(2,:), pdo(1,:), adj);
    
    L_e = L_all - L_a;  % extrinsic info.

    % Decoder two         
    L_a = L_e(alpha);  % a priori info.

    [pbit0, pbit1] = llr2pb(L_a);
    pdi = [pbit0; pbit1];
    [pco2, pdo]=siso_c_jxjnew(trellis, pci2, pdi);

    L_all = pb2llr(pdo(2,:), pdo(1,:), adj);
    
    L_e = L_all - L_a;
    % Estimate the info. bits        
%     rbit_code(alpha) = (1-sign(L_all))/2;
end
d(alpha) = (1+sign(L_all))/2;
L_code_out = 1;
% adj_out = 1e-100;
% 
% L_code_outtemp = zeros(3, L_total);
% L_code_outtemp(1, alpha) = L_all;
% 
% [pbit0, pbit1] = p_sym2bit(pco1);
% Lbo = pb2llr(pbit1, pbit0, adj_out);
% Lbo = ex_odd_even(Lbo); 
% L_code_outtemp(2, :) = Lbo(2:2:end);
% 
% [pbit0, pbit1] = p_sym2bit(pco2);
% Lbo = pb2llr(pbit1, pbit0, adj_out);
% Lbo = ex_odd_even(Lbo); 
% L_code_outtemp(3, :) = Lbo(2:2:end);
% 
% % L_code_outtemp = zeros(3, L_total);
% % L_code_outtemp(1, alpha) = L_e;
% % 
% % pbo1 = ProbSymbols2Bits(pco1, 2);
% % Lbo = log((pbo1(1,:)+adj_out)./(pbo1(2,:)+adj_out));
% % pbi1 = ProbSymbols2Bits(pci1, 2);
% % Lbi = log((pbi1(1,:)+adj_out)./(pbi1(2,:)+adj_out));
% % L_code_outtemp(2, :) = Lbo(2:2:end) - Lbi(2:2:end);
% % 
% % pbo2 = ProbSymbols2Bits(pco2, 2);
% % Lbo = log((pbo2(1,:)+adj_out)./(pbo2(2,:)+adj_out));
% % pbi2 = ProbSymbols2Bits(pci2, 2);
% % Lbi = log((pbi2(1,:)+adj_out)./(pbi2(2,:)+adj_out));
% % L_code_outtemp(3, :) = Lbo(2:2:end) - Lbi(2:2:end);
% 
% L_code_out = L_code_outtemp(:).';
