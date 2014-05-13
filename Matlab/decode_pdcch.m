function [pdcch_info, pcfich_info, pcfich_corr] = decode_pdcch(peak, subframe_idx, tfg)
pdcch_info = 0;

% subframe_idx = mod(subframe_idx, 10);
% start_slot_idx = subframe_idx*2;
% end_slot_idx = start_slot_idx + 1;

[pcfich_info, pcfich_corr] = decode_pcfich(peak, subframe_idx, tfg);

n_rb_dl = peak.n_rb_dl;

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

% phich_sc_idx = get_phich_sc_idx(n_id_cell, n_rb_dl, n_symb_dl);
