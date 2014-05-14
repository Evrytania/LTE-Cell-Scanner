function pdcch_info = decode_pdcch(peak, pcfich_info, subframe_idx, tfg, ce_tfg, np_ce)
pdcch_info = 0;

% subframe_idx = mod(subframe_idx, 10);
start_slot_idx = subframe_idx*2;
end_slot_idx = start_slot_idx + 1;

% Derive some values
n_ofdm = size(tfg,1);
n_rb_dl = peak.n_rb_dl;
nSC = n_rb_dl*12;
n_id_cell = peak.n_id_cell;
n_symb_dl = peak.n_symb_dl;
n_ports = peak.n_ports;
uldl_cfg = peak.uldl_cfg;

% decide number of ofdm symbols of phich and pdcch
if peak.phich_dur_value == 0 % normal
    n_phich_symb = 1;
elseif peak.phich_dur_value == 1 % extended
    if (peak.duplex_mode == 1) && ( subframe_idx==1 || subframe_idx==6 ) % TDD
        n_phich_symb = 2;
    else
        n_phich_symb = 3;
    end
else
    disp('Invalid peak.phich_dur_value!');
    return;
end

if n_rb_dl > 10
    if peak.phich_dur_value == 1 % extended
        n_pdcch_symb = n_phich_symb;
    else
        n_pdcch_symb = pcfich_info;
    end
else
    n_pdcch_symb = pcfich_info + 1;
end

[pdcch_sym, pdcch_ce] = pdcch_extract(peak, subframe_idx, tfg, ce_tfg, n_phich_symb, n_pdcch_symb);
