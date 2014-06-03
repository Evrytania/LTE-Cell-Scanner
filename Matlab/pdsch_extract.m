function [syms, ce_extract] = pdsch_extract(peak, reg_info, dci_info, subframe_idx, tfg, ce)

n_rb_dl = peak.n_rb_dl;
n_ports = peak.n_ports;
cp_type = peak.cp_type;
cp_type_val = peak.cp_type_val;
n_id_cell = peak.n_id_cell;

if cp_type_val == 0
  n_symb_dl=7;
elseif cp_type_val == 1
  n_symb_dl=6;
end

RB_start = dci_info.RB_start;

L_CRBs = dci_info.L_CRBs;
% L_CRBs = dci_info.N_1A_PRB; % only care about 1A SI-RNTI PDSCH now!  according to 36.212 5.3.3.1.3

sc_sp = RB_start*12;
sc_ep = (RB_start + L_CRBs - 1)*12 + 11;

n_pdcch_symb = reg_info.n_pdcch_symb;

num_ofdm_sym_total = size(tfg, 1);

% remove CRS
slot_num = subframe_idx*2;
for i = n_pdcch_symb : (num_ofdm_sym_total-1)
    
    if i == n_symb_dl
        slot_num = slot_num + 1;
    end
    
    for port_idx = 0 : (n_ports-1)
        [~, shift]=rs_dl(slot_num, mod(i,n_symb_dl), port_idx, n_id_cell, n_rb_dl, cp_type);
        if ~isnan(shift)
            sc_idx = (1+shift) : 6 : (n_rb_dl*12);
            tfg(i+1, sc_idx) = NaN;
        end
    end
end

syms = tfg( (n_pdcch_symb+1) : end, (sc_sp+1):(sc_ep+1)).';
syms = syms(:).';
ef_logical = ~isnan(syms);
syms = syms(ef_logical);

ce_extract = zeros(n_ports, length(syms));

for i = 1 : n_ports
    ce_tmp = ce( (n_pdcch_symb+1) : end, (sc_sp+1):(sc_ep+1), i).';
    ce_tmp = ce_tmp(:).';
    ce_tmp = ce_tmp(ef_logical);
    ce_extract(i, :) = ce_tmp;
end
