function syms = decode_pdsch(peak, reg_info, dci_info, subframe_idx, tfg, ce, np_ce, mmse_flag)

n_ports = peak.n_ports;

[syms, ce] = pdsch_extract(peak, reg_info, dci_info, subframe_idx, tfg, ce);

% mmse_flag = 0;
[syms, ~] = LTE_SFBC_demod(syms, ce, np_ce, n_ports, mmse_flag);
