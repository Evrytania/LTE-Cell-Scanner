function pdcch_info = decode_pdcch(peak, reg_info, subframe_idx, tfg, ce_tfg, np_ce)

% Derive some values
n_id_cell = peak.n_id_cell;
n_ports = peak.n_ports;

[pdcch_sym, pdcch_ce] = pdcch_extract(peak, reg_info, tfg, ce_tfg);

if isempty(pdcch_sym)
    pdcch_info = [];
    return;
end

% if isempty(pdcch_sym)
N_REG = length(pdcch_sym)/4;
N_CCE = floor(N_REG/9);

num_CCE = 16;  % common search space
L_set = [4 8]; % common search space
bits_set = [288 576];
Y = 0;         % common search space
pdcch_info = [];
for l = 1 : length(L_set)
    L = L_set(l);
    M = num_CCE/L;
    for m = 0 : (M-1)
        sc_idx = zeros(L, 36);
        for i = 0 : (L-1)
            CCE_idx = L*( mod( Y+m, floor(N_CCE/L) ) ) + i;
            sc_idx(i+1, :) = CCE_idx*36 : (CCE_idx*36 + 35);
        end
        sc_idx = sc_idx.';
        sc_idx = sc_idx(:).'; 
        
        sym = vec2mat( pdcch_sym(sc_idx+1), 4 ); % group into quadruplet
        M_quad = size(sym,1);
        ce = zeros(M_quad, size(sym,2), n_ports);
        for i = 1 : n_ports
            ce(:,:,i) = vec2mat( pdcch_ce(i, sc_idx+1), 4 );
        end
        
        sym = pdcch_de_cyclic_shift(sym, n_id_cell); % quadruplet idx is in the 1st dim!
        ce = pdcch_de_cyclic_shift(ce, n_id_cell); % quadruplet idx is in the 1st dim!
%         plot(abs(sym));
        tmp_info = pdcch_bit_level_proc(peak, subframe_idx*2, sym, ce, np_ce, bits_set(l));
        pdcch_info = [pdcch_info tmp_info];
    end
end
