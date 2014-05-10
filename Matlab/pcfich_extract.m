function [pcfich_sym, pcfich_ce]=pcfich_extract(peak, tfg, ce)

% Extract only the PCFICH RE's from the TFG

% Local shortcuts
cp_type=peak.cp_type;
n_id_cell=peak.n_id_cell;
n_rb_dl = peak.n_rb_dl;
n_ports = peak.n_ports;

% Derive some values
n_ofdm=size(tfg,1);
if (strcmpi(cp_type,'normal'))
  n_symb_dl=7;
elseif (strcmpi(cp_type,'extended'))
  n_symb_dl=6;
else
 error('Check code...');
end


symbol_idx = 1 : (2*n_symb_dl) : n_ofdm;

n_pcfich_symbol = length(symbol_idx);

pcfich_sym=NaN(1, n_pcfich_symbol*16);
pcfich_ce=NaN(n_ports, n_pcfich_symbol*16);

pcfich_sc_idx = get_pcfich_sc_idx(n_id_cell, n_rb_dl, cp_type);

for i = 1 : n_pcfich_symbol
    symbol_idx_tmp = symbol_idx(i);
    
    sp = (i-1)*16 + 1;
    ep = sp + 16 - 1;
    
    pcfich_sym(sp:ep) = tfg(symbol_idx_tmp, pcfich_sc_idx);
    
    for j = 1 : n_ports
        pcfich_ce(j, sp:ep) = ce(symbol_idx_tmp, pcfich_sc_idx, j);
    end
end

