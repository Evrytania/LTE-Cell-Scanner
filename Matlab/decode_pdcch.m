function pdcch_info = decode_pdcch(peak, reg_info, subframe_idx, tfg, ce_tfg, np_ce)
pdcch_info = 0;

% Derive some values
n_id_cell = peak.n_id_cell;
n_ports = peak.n_ports;

[pdcch_sym, pdcch_ce] = pdcch_extract(peak, reg_info, tfg, ce_tfg);

if isempty(pdcch_sym)
    pdcch_info = [];
    return;
end

% group into quadruplet -----------------
sym = vec2mat( pdcch_sym, 4 );
M_quad = size(sym,1);
ce = zeros(M_quad, size(sym,2), n_ports);
for i = 1 : n_ports
    ce(:,:,i) = vec2mat( pdcch_ce(i, :), 4 );
end

% de cyclic shifting -------------------
sym = pdcch_de_cyclic_shift(sym, n_id_cell); % quadruplet idx is in the 1st dim!
ce = pdcch_de_cyclic_shift(ce, n_id_cell); % quadruplet idx is in the 1st dim!

% de subblock interleave ---------------------
v = lte_subblock_interleave(1:M_quad);
v = v(~isnan(v));
deinterleave_order = zeros(1, M_quad);
deinterleave_order(v) = 1:M_quad;

sym = sym(deinterleave_order, :);
ce_deinterleave = ce;
for i = 1 : n_ports
    ce_deinterleave(:,:,i) = ce(deinterleave_order,:,i);
end

% de layer mapping and precoding blind search etc --------------
pdcch_info = pdcch_bit_level_proc(peak, subframe_idx*2, sym, ce, np_ce);
