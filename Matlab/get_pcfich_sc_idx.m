function idx = get_pcfich_sc_idx(n_id_cell, n_rb_dl, n_symb_dl)

n_ports = 1; % dummy
port_num = 0; % dummy
sym_num = 0; % dummy
slot_num = 0; % dummy
idx_rs_occupied = get_cell_specific_rs_occupied_idx(slot_num, sym_num, n_ports, port_num, n_id_cell, n_rb_dl, n_symb_dl);

nSC = n_rb_dl*12;

k_bar = 6*mod(n_id_cell, 2*n_rb_dl);

a = zeros(1, nSC);
for i=1:4
    reg_idx = mod( k_bar + floor((n_rb_dl/2)*(i-1))*6, nSC );
    a( (reg_idx+1) : (reg_idx+6)) = 1;
end

a(idx_rs_occupied) = 0;

idx = find(a);
