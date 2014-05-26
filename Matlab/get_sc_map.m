function [sc_map, reg_info] = get_sc_map(peak, pcfich_info, subframe_idx)
% sc_map is 1200*n_pdcch_symb matrix
% values definition:
% 1 -- rs
% 2 -- pcfich
% 3 -- phich
% 4 -- pdcch

% if pcfich_info<=0
%     pcfich_info = 3; % fake one to see
% end

n_rb_dl = peak.n_rb_dl;
n_ports = peak.n_ports;
cp_type = peak.cp_type;
n_id_cell = peak.n_id_cell;

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
% ---------------------------------------

phich_abs_reg_idx_pre_calc = get_phich_abs_reg_idx(peak, subframe_idx, n_phich_symb);

reg_info.n_phich_symb = n_phich_symb;
reg_info.n_pdcch_symb = n_pdcch_symb;
reg_info.reg_map = cell(1, n_pdcch_symb);
reg_info.num_reg_all = zeros(1, n_pdcch_symb);
reg_info.pcfich_abs_reg_idx = cell(1, n_pdcch_symb);
reg_info.phich_abs_reg_idx = cell(1, n_pdcch_symb);
reg_info.pdcch_abs_reg_idx = cell(1, n_pdcch_symb);

sc_map = zeros(n_rb_dl*12, n_pdcch_symb);
slot_num = subframe_idx*2; % control domain only exist in the 1st slot of each subframe
for i = 0 : (n_pdcch_symb-1)
    
    % rs
    for port_idx = 0 : (n_ports-1)
        [~, shift]=rs_dl(slot_num, i, port_idx, n_id_cell, n_rb_dl, cp_type);
        if ~isnan(shift)
            sc_idx = (1+shift) : 6 : (n_rb_dl*12);
            sc_map(sc_idx, i+1) = 1;
        end
    end
    
    % pcfich only in the 1st ofdm symbol
    pcfich_abs_reg_idx = get_pcfich_abs_reg_idx(peak, i);
    sc_idx = conv_abs_reg_idx_to_sc_idx( peak, i, pcfich_abs_reg_idx );
    sc_map(sc_idx+1, i+1) = 2;
    
    % phich
    if i < n_phich_symb
        phich_abs_reg_idx = phich_abs_reg_idx_pre_calc{i+1};
        sc_idx = conv_abs_reg_idx_to_sc_idx( peak, i, phich_abs_reg_idx );
        sc_map(sc_idx+1, i+1) = 3;
    else
        phich_abs_reg_idx = [];
    end
    
    % pdcch
    num_reg_all = get_num_ctrl_reg_all(peak, i);
    pdcch_abs_reg_idx = 0 : (num_reg_all-1);
    pdcch_abs_reg_idx([pcfich_abs_reg_idx, phich_abs_reg_idx]+1) = [];
    sc_idx = conv_abs_reg_idx_to_sc_idx( peak, i, pdcch_abs_reg_idx );
    sc_map(sc_idx+1, i+1) = 4;
    
    % reg_info
    reg_map = zeros(1, num_reg_all);
    reg_map(pcfich_abs_reg_idx+1) = 1;
    reg_map(phich_abs_reg_idx+1) = 2;
    reg_info.reg_map{i+1} = reg_map;
    reg_info.num_reg_all(i+1) = num_reg_all;
    reg_info.pcfich_abs_reg_idx{i+1} = pcfich_abs_reg_idx;
    reg_info.phich_abs_reg_idx{i+1} = phich_abs_reg_idx;
    reg_info.pdcch_abs_reg_idx{i+1} = pdcch_abs_reg_idx;
end
