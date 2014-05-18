function [pdcch_sym, pdcch_ce, reg_info]=pdcch_extract(peak, subframe_idx, tfg, ce, n_phich_symb, n_pdcch_symb)
% only operate on the first slot

n_ports = peak.n_ports;

pdcch_sc_idx_store = cell(1, n_pdcch_symb);
pdcch_sym_count = 0;

phich_abs_reg_idx_pre_calc = get_phich_abs_reg_idx(peak, subframe_idx, n_phich_symb);

reg_info.reg_map = cell(1, n_pdcch_symb);
reg_info.num_reg_all = zeros(1, n_pdcch_symb);
reg_info.pcfich_abs_reg_idx = cell(1, n_pdcch_symb);
reg_info.phich_abs_reg_idx = cell(1, n_pdcch_symb);
reg_info.pdcch_abs_reg_idx = cell(1, n_pdcch_symb);
for i = 0 : (n_pdcch_symb-1)
    num_reg_all = get_num_ctrl_reg_all(peak, i);
    
    pcfich_abs_reg_idx = get_pcfich_abs_reg_idx(peak, i);
    if i<n_phich_symb
        phich_abs_reg_idx = phich_abs_reg_idx_pre_calc{i+1};
    else
        phich_abs_reg_idx = [];
    end
    
    pdcch_abs_reg_idx = 0 : (num_reg_all-1);
    pdcch_abs_reg_idx([pcfich_abs_reg_idx, phich_abs_reg_idx]+1) = [];
    pdcch_sc_idx = conv_abs_reg_idx_to_sc_idx( peak, i, pdcch_abs_reg_idx );
    
    pdcch_sc_idx_store{i+1} = pdcch_sc_idx;
    pdcch_sym_count = pdcch_sym_count + length(pdcch_sc_idx);

    reg_map = zeros(1, num_reg_all);
    reg_map(pcfich_abs_reg_idx+1) = 1;
    reg_map(phich_abs_reg_idx+1) = 2;
    reg_info.reg_map{i+1} = reg_map;
    reg_info.num_reg_all(i+1) = num_reg_all;
    reg_info.pcfich_abs_reg_idx{i+1} = pcfich_abs_reg_idx;
    reg_info.phich_abs_reg_idx{i+1} = phich_abs_reg_idx;
    reg_info.pdcch_abs_reg_idx{i+1} = pdcch_abs_reg_idx;
end

pdcch_sym = zeros(1, pdcch_sym_count);
pdcch_ce = zeros(n_ports, pdcch_sym_count);

pdcch_sym_count = 0;
for i=1:n_pdcch_symb
    pdcch_sc_idx = pdcch_sc_idx_store{i};
    tmp_count = length(pdcch_sc_idx);
    pdcch_sym( (pdcch_sym_count+1):(pdcch_sym_count+tmp_count)  ) = tfg(i, pdcch_sc_idx+1);
    
    for j = 1 : n_ports
        pdcch_ce(j, (pdcch_sym_count+1):(pdcch_sym_count+tmp_count) ) = ce(i, pdcch_sc_idx+1, j);
    end
    
    pdcch_sym_count = pdcch_sym_count + tmp_count;
end
