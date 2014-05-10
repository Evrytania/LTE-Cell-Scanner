function idx = get_cell_specific_rs_occupied_idx(slot_num, sym_num, n_ports, port_num, n_id_cell, n_rb_dl, cp_type)
% include both real RS and dummy RS

nSC = n_rb_dl*12;

if (strcmpi(cp_type,'normal'))
  n_symb_dl=7;
elseif (strcmpi(cp_type,'extended'))
  n_symb_dl=6;
else
  error('Unrecognized cp_type specified');
end

if ((port_num==0)||(port_num==1))
  if ((sym_num~=0)&&(sym_num~=n_symb_dl-3))
      disp('Ant port 0 or 1 only have cell specific RS in ofdm symbol 0 and n_symb_dl-3!');
      return;
  end
end

if ((port_num==2)||(port_num==3))
  if (sym_num~=1)
      disp('Ant port 2 or 3 only have cell specific RS in ofdm symbol 1!');
      return;
  end
end

slot_num_mod2 = mod(slot_num,2);

if ((port_num==0)&&(sym_num==0))
  v=0;
elseif ((port_num==0)&&(sym_num~=0))
  v=3;
elseif ((port_num==1)&&(sym_num==0))
  v=3;
elseif ((port_num==1)&&(sym_num~=0))
  v=0;
elseif (port_num==2)
  v=3*slot_num_mod2;
elseif (port_num==3)
  v=3+3*slot_num_mod2;
end

v_shift=mod(n_id_cell,6);

shift=mod(v+v_shift,6);

idx1 = (1+shift):6:nSC;

if (n_ports == 1) && ( slot_num_mod2==1 || ( slot_num_mod2==0 && sym_num~=0 ) )
    idx2 = [];
else
    idx2 = mod( (idx1-1) + 3, nSC ) + 1;
end

idx = sort([idx1, idx2]);
