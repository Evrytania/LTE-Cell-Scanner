function sig_tx=create_dl_sig(cp_type,n_subframes,slot_start,n_id_1,n_id_2,load_factor)

% Create an eNodeB DL signal containing PSS & SSS

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

error(nargchk(6,6,nargin));
error(chk_param(cp_type,'cp_type','string'));
error(chk_param(n_subframes,'n_subframes','scalar','real','integer','>',0));
error(chk_param(slot_start,'slot_start','scalar','real','integer','>=',0,'<=',19));
error(chk_param(n_id_1,'n_id_1','scalar','real','integer','>=',0,'<=',167));
error(chk_param(n_id_2,'n_id_2','scalar','real','integer','>=',0,'<=',2));
error(chk_param(load_factor,'load_factor','scalar','real','>=',0,'<=',1));

% These are constants that cannot change
n_dft=128;
n_rb=6;

% Derive some values
n_sc=n_rb*12;
if strcmpi(cp_type,'normal')
  n_ofdm=7;
elseif strcmpi(cp_type,'extended')
  n_ofdm=6;
else
  error('Check code...');
end
n_id_cell=n_id_2+3*n_id_1;

% Create the wideband transmitted signal
sig_tx=NaN(1,n_subframes/1000*(fs_lte/16));
offset=1;
for t=1:2*n_subframes
  slot_num=mod(slot_start+t-1,20);

  for k=0:n_ofdm-1
    % Select RS
    [p0 s0]=rs_dl(slot_num,k,0,n_id_cell,6,cp_type);
    [p1 s1]=rs_dl(slot_num,k,1,n_id_cell,6,cp_type);
    if (~isnan(s0))
      rs_ind=[(s0+1:6:n_sc) (s1+1:6:n_sc)];
    else
      rs_ind=[];
    end
    rs=[p0 p1];
    nrs_cand=setdiff(1:n_sc,rs_ind);
    nrs_ind=randperm(length(nrs_cand));
    nrs_ind=nrs_cand(nrs_ind(1:round(length(nrs_ind)*load_factor)));

    % Create subcarrier symbols for this OFDM symbol
    syms=zeros(1,n_sc);
    syms(rs_ind)=rs;
    temp=qam(bitstream(2*length(nrs_ind)),'qam');
    syms(nrs_ind)=temp;

    % Prepare the input to the idft operation
    idft_in=zeros(1,n_dft);
    % DC subcarrier remains at zero power
    idft_in(2:2+n_sc/2-1)=syms(n_sc/2+1:end);
    idft_in(end-n_sc/2+1:end)=syms(1:n_sc/2);

    % Replace some subcarriers with PSS or SSS as necessary
    if ((mod(slot_num,10)==0)&&(k>=n_ofdm-2))
      if (k==n_ofdm-1)
        overwrite=pss(n_id_2);
      else
        overwrite=sss(n_id_1,n_id_2,slot_num);
      end
      idft_in(2:37)=[overwrite(32:end) zeros(1,5)];
      idft_in(end-36+1:end)=[zeros(1,5) overwrite(1:31)];
    end

    % Create time domain signal and add CP
    td=idft(idft_in);
    if (strcmpi(cp_type,'extended'))
      cp_length=32;
    elseif (strcmpi(cp_type,'normal'))
      if (k==0)
        cp_length=10;
      else
        cp_length=9;
      end
    else
      error('Check code');
    end
    td=[td(end-cp_length+1:end) td];
    %if((mod(slot_num,10)==0)&&(k==6))
    %  return
    %end

    % Add onto the transmitted signal vector.
    sig_tx(offset:offset+length(td)-1)=td;
    offset=offset+length(td);
  end

end

n_samp_tx=length(sig_tx);
assert(n_samp_tx==n_subframes*.001*(fs_lte/16));

