function [pdcch_sym, pdcch_ce, pdcch_sc_idx_store]=pdcch_extract(peak, subframe_idx, tfg, ce, n_phich_symb, n_pdcch_symb)
% only operate on the first slot

n_ports = peak.n_ports;
n_symb_dl = peak.n_symb_dl;
n_id_cell = peak.n_id_cell;
n_rb_dl = peak.n_rb_dl;
nSC = 12*n_rb_dl;
% uldl_cfg = peak.uldl_cfg;

pcfich_idx_set = get_pcfich_sc_idx(n_id_cell, n_rb_dl, n_symb_dl); % only valid for the first ofdm symbol

slot_num = 0;
sym_num = 0;
port_num = 0;
idx_rs_occupied1 = get_cell_specific_rs_occupied_idx(slot_num, sym_num, n_ports, port_num, n_id_cell, n_rb_dl, n_symb_dl);
idx_rs_occupied1 = kron(ones(n_ports,1), idx_rs_occupied1);
for i = 2 : n_ports
    port_num = i-1;
    idx_rs_occupied1(i,:) = get_cell_specific_rs_occupied_idx(slot_num, sym_num, n_ports, port_num, n_id_cell, n_rb_dl, n_symb_dl);
end
idx_rs_occupied1 = idx_rs_occupied1.';
idx_rs_occupied1 = idx_rs_occupied1(:).';

sym_num = n_symb_dl-3;
port_num = 0;
idx_rs_occupied2 = get_cell_specific_rs_occupied_idx(slot_num, sym_num, n_ports, port_num, n_id_cell, n_rb_dl, n_symb_dl);
idx_rs_occupied2 = kron(ones(n_ports,1), idx_rs_occupied2);
for i = 2 : n_ports
    port_num = i-1;
    idx_rs_occupied2(i,:) = get_cell_specific_rs_occupied_idx(slot_num, sym_num, n_ports, port_num, n_id_cell, n_rb_dl, n_symb_dl);
end
idx_rs_occupied2 = idx_rs_occupied2.';
idx_rs_occupied2 = idx_rs_occupied2(:).';

pdcch_sc_idx_store = cell(1, n_pdcch_symb);
pdcch_sym_count = 0;
for i = 0 : (n_pdcch_symb-1)

    if i == 0
        idx_rs_occupied = idx_rs_occupied1;
    elseif i == (n_symb_dl-3)
        idx_rs_occupied = idx_rs_occupied2;
    else
        idx_rs_occupied = [];
    end
    
    if i==0
        pcfich_sc_idx = pcfich_idx_set;
    else
        pcfich_sc_idx = [];
    end

    if i < n_phich_symb
        phich_sc_idx = get_phich_sc_idx(peak, subframe_idx, i);
    else
        phich_sc_idx = [];
    end

    pdcch_sc_idx = ones(1, nSC);
    pdcch_sc_idx([idx_rs_occupied, pcfich_sc_idx, phich_sc_idx]) = 0;
    pdcch_sc_idx = find(pdcch_sc_idx);
    
    pdcch_sc_idx_store{i+1} = pdcch_sc_idx;
    pdcch_sym_count = pdcch_sym_count + length(pdcch_sc_idx);
end

pdcch_sym = zeros(1, pdcch_sym_count);
pdcch_ce = zeros(n_ports, pdcch_sym_count);

pdcch_sym_count = 0;
for i=1:n_pdcch_symb
    pdcch_sc_idx = pdcch_sc_idx_store{i};
    tmp_count = length(pdcch_sc_idx);
    pdcch_sym( (pdcch_sym_count+1):(pdcch_sym_count+tmp_count)  ) = tfg(i, pdcch_sc_idx);
    
    for j = 1 : n_ports
        pdcch_ce(j, (pdcch_sym_count+1):(pdcch_sym_count+tmp_count) ) = ce(i, pdcch_sc_idx, j);
    end
    
    pdcch_sym_count = pdcch_sym_count + tmp_count;
end
