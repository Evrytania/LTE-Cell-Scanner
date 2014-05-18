function [pdcch_info, reg_info] = decode_pdcch(peak, pcfich_info, subframe_idx, tfg, ce_tfg, np_ce)
pdcch_info = 0;
if pcfich_info == 0 % no PDCCH in this subframe
    reg_info.reg_map = [];
    reg_info.num_reg_all = -1;
    reg_info.pcfich_abs_reg_idx = [];
    reg_info.phich_abs_reg_idx = [];
    reg_info.pdcch_abs_reg_idx = [];

    pdcch_info = NaN;
    return;
end

% subframe_idx = mod(subframe_idx, 10);
% start_slot_idx = subframe_idx*2;
% end_slot_idx = start_slot_idx + 1;

% Derive some values
% n_ofdm = size(tfg,1);
n_rb_dl = peak.n_rb_dl;
% nSC = n_rb_dl*12;
% n_id_cell = peak.n_id_cell;
% n_symb_dl = peak.n_symb_dl;
% n_ports = peak.n_ports;
% uldl_cfg = peak.uldl_cfg;

% decide number of ofdm symbols of phich and pdcch
if peak.phich_dur_val == 0 % normal
    n_phich_symb = 1;
elseif peak.phich_dur_val == 1 % extended
    if (peak.duplex_mode == 1) && ( subframe_idx==1 || subframe_idx==6 ) % TDD
        n_phich_symb = 2;
    else
        n_phich_symb = 3;
    end
else
    disp('Invalid peak.phich_dur_val! It must be 0 or 1!');
    return;
end

if n_rb_dl > 10
    if peak.phich_dur_val == 1 % extended
        n_pdcch_symb = n_phich_symb;
    else
        n_pdcch_symb = pcfich_info;
    end
else
    n_pdcch_symb = pcfich_info + 1;
end

[pdcch_sym, pdcch_ce, reg_info] = pdcch_extract(peak, subframe_idx, tfg, ce_tfg, n_phich_symb, n_pdcch_symb);

% scatterplot(pdcch_sym./pdcch_ce);
