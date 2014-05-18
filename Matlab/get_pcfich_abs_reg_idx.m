function pcfich_abs_reg_idx = get_pcfich_abs_reg_idx(peak, sym_idx)
% reg idx starts from 0
% only when sym_idx == 0 there is pcfich

if sym_idx == 0
    n_rb_dl = peak.n_rb_dl;
    n_id_cell = peak.n_id_cell;

    nSC = n_rb_dl*12;

    k_bar = 6*mod(n_id_cell, 2*n_rb_dl);

    pcfich_abs_reg_idx =zeros(1, 4);

    for i=1:4
        reg_idx = mod( k_bar + floor((n_rb_dl/2)*(i-1))*6, nSC );
        pcfich_abs_reg_idx(i) = reg_idx/6;
    end
else
    pcfich_abs_reg_idx = [];
end