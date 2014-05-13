function [pdcch_info, pcfich_info, pcfich_corr] = decode_pdcch(peak, subframe_idx, tfg)
pdcch_info = 0;

% subframe_idx = mod(subframe_idx, 10);
% start_slot_idx = subframe_idx*2;
% end_slot_idx = start_slot_idx + 1;

[pcfich_info, pcfich_corr] = decode_pcfich(peak, subframe_idx, tfg);
