function [pbch_sym pbch_ce]=pbch_extract(peak,tfg,ce);

% Extract only the MIB RE's from the TFG

% Copyright 2012 Evrytania LLC (http://www.evrytania.com)
%
% Written by James Peroulas <james@evrytania.com>
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Affero General Public License for more details.
%
% You should have received a copy of the GNU Affero General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Local shortcuts
n_id_1=peak.n_id_1;
n_id_2=peak.n_id_2;
cp_type=peak.cp_type;

% Derive some values
n_ofdm=size(tfg,1);
if (strcmpi(cp_type,'normal'))
  n_symb_dl=7;
  m_bit=1920;
elseif (strcmpi(cp_type,'extended'))
  n_symb_dl=6;
  m_bit=1728;
else
 error('Check code...');
end
n_id_cell=n_id_2+3*n_id_1;
v_shift=mod(n_id_cell,6);

pbch_sym=NaN(1,m_bit/2);
pbch_ce=NaN(4,m_bit/2);
idx=1;
v_shift_m3=mod(v_shift,3);
for fr=0:3
  for sym=0:3
    for sc=0:71
      % Skip (possible) RS locations
      if ((mod(sc,3)==v_shift_m3)&&((sym==0)||(sym==1)||((sym==3)&&(n_symb_dl==6))))
        continue
      end
      sym_num=fr*10*2*n_symb_dl+n_symb_dl+sym+1;
      pbch_sym(idx)=tfg(sym_num,sc+1);
      pbch_ce(1,idx)=ce(sym_num,sc+1,1);
      pbch_ce(2,idx)=ce(sym_num,sc+1,2);
      pbch_ce(3,idx)=ce(sym_num,sc+1,3);
      pbch_ce(4,idx)=ce(sym_num,sc+1,4);
      idx=idx+1;
    end
  end
end
assert(idx-1==m_bit/2);

