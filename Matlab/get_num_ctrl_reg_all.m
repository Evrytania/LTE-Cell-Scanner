function num_reg_all = get_num_ctrl_reg_all(peak, sym_idx)
% this function only valid for the 1st slot of subframe
% because control domain only exists in the 1st slot of subframe

n_rb_dl = peak.n_rb_dl;
n_ports = peak.n_ports;
cp_type_val = peak.cp_type_val;

if sym_idx == 0
    reg_x = 2;
elseif sym_idx == 1
    if n_ports == 1 || n_ports == 2
        reg_x = 3;
    elseif n_ports == 4
        reg_x = 2;
    else
        disp('Invalid n_ports! It must be 1 or 2 or 4!');
        return;
    end
elseif sym_idx == 2
    reg_x = 3;
elseif sym_idx == 3
    if cp_type_val == 0
        reg_x = 3;
    elseif cp_type_val == 1
        reg_x = 2;
    else
        disp('Invalid cp_type_val! It must be 0 or 1!');
        return;
    end
else
    disp('Invalid sym_idx! It must be 0~3!');
    return;
end

num_reg_all = reg_x*n_rb_dl;
