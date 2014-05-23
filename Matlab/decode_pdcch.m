function pdcch_info = decode_pdcch(peak, reg_info, subframe_idx, tfg, ce_tfg, np_ce)
pdcch_info = 0;

% Derive some values
n_id_cell = peak.n_id_cell;
n_ports = peak.n_ports;

[pdcch_sym, pdcch_ce] = pdcch_extract(peak, reg_info, tfg, ce_tfg, subframe_idx);

if isempty(pdcch_sym)
    pdcch_info = [];
    return;
end

% % -----------------------------------
% enb.CellRefP = n_ports;
% enb.NCellID = n_id_cell;
% enb.NSubframe = subframe_idx;
% enb.NDLRB = peak.n_rb_dl;
% if peak.cp_type_val == 0
%     enb.CyclicPrefix = 'Normal';
% else
%     enb.CyclicPrefix = 'Extended';
% end
% enb.CFI = 3;
% if peak.phich_res == 1
%     enb.Ng = 'One';
% elseif peak.phich_res == 2
%     enb.Ng = 'Two';
% elseif peak.phich_res == 0.5
%     enb.Ng = 'Half';
% elseif peak.phich_res == 1/6
%     enb.Ng = 'Sixth';
% end
% if peak.phich_dur_val == 0
%     enb.PHICHDuration = 'Normal';
% else
%     enb.PHICHDuration = 'Extended';
% end
% enb.NFrame = peak.sfn;
% enb.NSubframe = subframe_idx;
% if peak.duplex_mode == 0
%     enb.DuplexMode = 'FDD';
% elseif peak.duplex_mode == 1
%     enb.DuplexMode = 'TDD';
% end
% chs.RNTI = 0;
% [softbits,symbols] = ltePDCCHDecode(enb,pdcch_sym(:),pdcch_ce(:),mean(np_ce));
% [dcistr,dcibits] = ltePDCCHSearch(enb,chs,softbits);
% pdcch_info = dcistr;
% return;

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

% out = ltePDCCHDeinterleave(enb, pdcch_sym(:));
% sym_tmp = sym.';
% sym_tmp = sym_tmp(:);
% figure; plot(abs(out - sym_tmp));

% de layer mapping and precoding blind search etc --------------
pdcch_info = pdcch_bit_level_proc(peak, subframe_idx*2, sym, ce, np_ce);
