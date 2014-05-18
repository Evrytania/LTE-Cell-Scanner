function phich_sc_idx = get_phich_sc_idx(peak, subframe_idx, sym_idx)

% shortcuts
duplex_mode = peak.duplex_mode;
phich_res = peak.phich_res;
cp_type_val = peak.cp_type_val;
n_rb_dl = peak.n_rb_dl;
nSC = n_rb_dl*12;
uldl_cfg = peak.uldl_cfg;
n_id_cell = peak.n_id_cell;
phich_dur_val = peak.phich_dur_val;

Ng = phich_res;

tmp_val = ceil(Ng*(n_rb_dl/8));
if cp_type_val == 0
    N_group = tmp_val;
elseif cp_type_val == 1
    N_group = 2*tmp_val;
else
    disp('Invalid cp_type_val! Must be 0 or 1!');
    return;
end

if duplex_mode == 1 % TDD
    mi_table = [ ...
        2  1 -1 -1 -1  2  1 -1 -1 -1; ...
        0  1 -1 -1  1  0  1 -1 -1  1; ...
        0  0 -1  1  0  0  0 -1  1  0; ...
        1  0 -1 -1 -1  0  0  0  1  1; ...
        0  0 -1 -1  0  0  0  0  1  1; ...
        0  0 -1  0  0  0  0  0  1  0; ...
        1  1 -1 -1 -1  1  1 -1 -1  1
        ];
    mi = mi_table(uldl_cfg+1, subframe_idx+1);
end

if cp_type_val == 0
    if duplex_mode == 0
        max_mp = N_group - 1;
    elseif duplex_mode == 1
        max_mp = (mi*N_group) - 1;
    else
        disp('Invalid duplex_mode! It must be 0(FDD) or 1(TDD)!');
        return;
    end
else
    if duplex_mode == 0
        max_mp = (N_group/2) - 1;
    elseif duplex_mode == 1
        max_mp = (mi*N_group/2) - 1;
    else
        disp('Invalid duplex_mode! It must be 0(FDD) or 1(TDD)!');
        return;
    end
end

all_sc_idx = zeros(1, nSC);
num_reg_all = get_num_ctrl_reg_all(peak, sym_idx);
pcfich_abs_reg_idx = get_pcfich_abs_reg_idx(peak, sym_idx);

n_reg = [num_reg_all-4, num_reg_all, num_reg_all];

% table used to convert relative phich reg idx to absolute reg idx
phich_abs_reg_table = 0:(num_reg_all-1);
phich_abs_reg_table(pcfich_abs_reg_idx+1) = [];

for mp = 0 : max_mp
    for i = 0 : 2
        if phich_dur_val == 0
            lip = 0;
        elseif (phich_dur_val == 1) && (duplex_mode == 1) && (subframe_idx==1 || subframe_idx==6)
            lip = mod( floor(mp/2) + i + 1 , 2);
        else
            lip = i;
        end
        if sym_idx == lip % if it falls into current ofdm symbol
            n_reg_lip = n_reg(lip+1);
            if (phich_dur_val == 1) && (duplex_mode == 1) && (subframe_idx==1 || subframe_idx==6)
                n_reg_0 = n_reg(0+1);
                nibar = mod( floor(n_id_cell*n_reg_lip/n_reg_0) + mp + floor(i*n_reg_lip/3) , n_reg_lip);
            else
                n_reg_1 = n_reg(1+1);
                nibar = mod( floor(n_id_cell*n_reg_lip/n_reg_1) + mp + floor(i*n_reg_lip/3) , n_reg_lip);
            end
            sc_idx = conv_abs_reg_idx_to_sc_idx( peak, sym_idx, phich_abs_reg_table(nibar) );
            all_sc_idx(sc_idx+1) = 1;
        end
    end
end

phich_sc_idx = find(all_sc_idx);
