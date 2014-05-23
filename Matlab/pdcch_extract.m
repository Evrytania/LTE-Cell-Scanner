function [pdcch_sym, pdcch_ce]=pdcch_extract(peak, reg_info, tfg, ce, subframe_idx)
% only operate on the first slot

n_ports = peak.n_ports;
n_rb_dl = peak.n_rb_dl;
nSC = 12*n_rb_dl;
n_pdcch_symb = reg_info.n_pdcch_symb;

if n_pdcch_symb == 0
    pdcch_sym = [];
    pdcch_ce = [];
    return;
end

pdcch_sym_count = 0;

max_num_reg = -1;
for i = 0 : (n_pdcch_symb-1)
    num_pdcch_reg = length( reg_info.pdcch_abs_reg_idx{i+1} );
    pdcch_sym_count = pdcch_sym_count + num_pdcch_reg*4;
    if num_pdcch_reg > max_num_reg
        max_num_reg = num_pdcch_reg;
    end
end

pdcch_sym = zeros(1, pdcch_sym_count);
pdcch_ce = zeros(n_ports, pdcch_sym_count);

abs_reg_idx_mat = -ones(n_pdcch_symb, max_num_reg+1);
reg_count = zeros(1, n_pdcch_symb);
for i = 1 : n_pdcch_symb
    pdcch_abs_reg_idx = reg_info.pdcch_abs_reg_idx{i};
    tmp_len = length(pdcch_abs_reg_idx);
    abs_reg_idx_mat(i, 1:tmp_len) = pdcch_abs_reg_idx;
end

pdcch_sym_count = 0;
for k = 0 : (nSC - 1)
    for i = 0 : (n_pdcch_symb-1)

        abs_reg_idx = abs_reg_idx_mat(i+1, reg_count(i+1)+1 );
        [sc_idx, sc_idx_all] = conv_abs_reg_idx_to_sc_idx( peak, i, abs_reg_idx );
        if sum(k==sc_idx_all) > 0
            pdcch_sym( (pdcch_sym_count+1):(pdcch_sym_count+4)  ) = tfg(i+1, sc_idx+1);

            for j = 1 : n_ports
                pdcch_ce(j, (pdcch_sym_count+1):(pdcch_sym_count+4) ) = ce(i+1, sc_idx+1, j);
            end
            
            pdcch_sym_count = pdcch_sym_count + 4;
            
            reg_count(i+1) = reg_count(i+1) + 1;
        end

    end
end
