function [pcfich_sym, pcfich_ce]=pcfich_extract(peak, tfg, ce)

% Extract only the PCFICH RE's from the TFG

% Local shortcuts
% n_id_cell=peak.n_id_cell;
% n_rb_dl = peak.n_rb_dl;
n_ports = peak.n_ports;

% Derive some values
n_ofdm=size(tfg,1);
n_symb_dl = peak.n_symb_dl;

symbol_idx = 1 : (2*n_symb_dl) : n_ofdm;

n_pcfich_symbol = length(symbol_idx);

pcfich_sym=NaN(n_pcfich_symbol, 16);
pcfich_ce=NaN(n_pcfich_symbol, 16, n_ports);

% pcfich_sc_idx = get_pcfich_sc_idx(n_id_cell, n_rb_dl, n_symb_dl);
pcfich_abs_reg_idx = get_pcfich_abs_reg_idx(peak, 0);
pcfich_sc_idx = conv_abs_reg_idx_to_sc_idx( peak, 0, pcfich_abs_reg_idx );

for i = 1 : n_pcfich_symbol
    symbol_idx_tmp = symbol_idx(i);

    pcfich_sym(i,:) = tfg(symbol_idx_tmp, pcfich_sc_idx+1);
    
    for j = 1 : n_ports
        pcfich_ce(i,:,j) = ce(symbol_idx_tmp, pcfich_sc_idx+1, j);
    end
end

