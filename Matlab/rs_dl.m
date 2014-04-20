function [r shift]=rs_dl(slot_num,sym_num,port_num,n_id_cell,n_rb_dl,cp_type)

% [r shift]=rs_dl(slot_num,sym_num,port_num,n_id_cell,n_rb_dl,cp_type)
%
% This function supplies the reference symbols to be used on OFDM symbol
% 'sym_num' of slot 'slot_num' for antenna port 'port_num'.
%
% The reference symbols should be located on subcarriers (1+shift):6:n_rb_dl*12
%
% If a particular OFDM symbol of a particular port does not contain any
% reference symbols, this function will return the empty matrix in 'r'
% and NaN in 'shift'.
%
% 'cp_type' can be either 'normal' or 'extended'.

% Copyright 2012 Evrytania LLC (http://www.evrytania.com)
%
% Written by James Peroulas <james@evrytania.com>
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Constants
n_rb_maxdl=110;

error(nargchk(6,6,nargin));
error(chk_param(slot_num,'slot_num','scalar','real','integer','>=',0,'<=',19));
error(chk_param(sym_num,'sym_num','scalar','real','integer','>=',0,'<=',6));
error(chk_param(port_num,'port_num','scalar','real','integer','>=',0,'<=',3));
error(chk_param(n_id_cell,'n_id_cell','scalar','real','integer','>=',0,'<=',503));
error(chk_param(n_rb_dl,'n_rb_dl','scalar','real','integer','>=',0,'<=',110));
error(chk_param(cp_type,'cp_type','string'));

% Second level processing
if (strcmpi(cp_type,'normal'))
  n_symb_dl=7;
  n_cp=1;
elseif (strcmpi(cp_type,'extended'))
  n_symb_dl=6;
  n_cp=0;
else
  error('Unrecognized cp_type specified');
end

% Quick return in some cases
r=[];
shift=NaN;
if ((port_num==0)||(port_num==1))
  if ((sym_num~=0)&&(sym_num~=n_symb_dl-3))
    return
  end
end
if ((port_num==2)||(port_num==3))
  if (sym_num~=1)
    return
  end
end

% Create the source sequence.
c_init=2^10*(7*(slot_num+1)+sym_num+1)*(2*n_id_cell+1)+2*n_id_cell+n_cp;
c=lte_pn(c_init,4*n_rb_maxdl);

% Create symbols to be used if n_rb_maxdl RB's were used on the DL.
r_l_ns=(1/sqrt(2))*complex(1-2*c(1:2:end),1-2*c(2:2:end));

r=r_l_ns(1+n_rb_maxdl-n_rb_dl:2*n_rb_dl+n_rb_maxdl-n_rb_dl);

if ((port_num==0)&&(sym_num==0))
  v=0;
elseif ((port_num==0)&&(sym_num~=0))
  v=3;
elseif ((port_num==1)&&(sym_num==0))
  v=3;
elseif ((port_num==1)&&(sym_num~=0))
  v=0;
elseif (port_num==2)
  v=3*mod(slot_num,2);
elseif (port_num==3)
  v=3+3*mod(slot_num,2);
end

v_shift=mod(n_id_cell,6);

shift=mod(v+v_shift,6);

