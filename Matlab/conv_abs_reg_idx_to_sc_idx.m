function sc_idx = conv_abs_reg_idx_to_sc_idx( peak, sym_idx, abs_reg_idx )
% this function only valid for the 1st slot of subframe
% because control domain only exists in the 1st slot of subframe

if isempty(abs_reg_idx)
    sc_idx = [];
    return;
end

if abs_reg_idx<0
    sc_idx = [];
    return;
end

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

num_abs_reg_idx = length(abs_reg_idx);
if reg_x == 3
    sp = abs_reg_idx(:).*4;
    sc_idx = kron(ones(1,4), sp) + kron( ones(num_abs_reg_idx,1), 0:3 );
%     for i = 1 : num_abs_reg_idx
%         sc_idx(i,:) = sp(i) : (sp(i) + 3);
%     end
elseif reg_x == 2
    sp = abs_reg_idx(:).*6;
    sc_idx = kron(ones(1,4), sp) + kron( ones(num_abs_reg_idx,1), [1 2 4 5] );
%     for i = 1 : num_abs_reg_idx
%         sc_idx(i,:) = [sp(i)+1, sp(i)+2, sp(i)+4, sp(i)+5];
%     end
end

sc_idx = sc_idx.';
sc_idx = sc_idx(:).';
