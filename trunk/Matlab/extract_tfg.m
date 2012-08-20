function [tfg tfg_timestamp]=extract_tfg(peak,capbuf,fc);

% Convert from time domain to frequency domain and create the time/frequency
% grid.

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

frame_start=peak.frame_start;
cp_type=peak.cp_type;
freq_fine=peak.freq_fine;

% fc*k_factor is the receiver's actual RX center frequency.
k_factor=(fc-peak.freq_fine)/fc;

if (strcmpi(cp_type,'normal'))
  n_symb_dl=7;
  dft_location=frame_start+10;
elseif (strcmpi(cp_type,'extended'))
  n_symb_dl=6;
  dft_location=frame_start+16;
else
 error('Check code...');
end

% See if we can advance the frame start
if (dft_location-k_factor*fs_lte/16*.01>=0.5)
  dft_location=dft_location-k_factor*fs_lte/16*.01;
end

% Perform FOC
capbuf=fshift(capbuf,-freq_fine,fs_lte/16);

% Extract 6 frames + 2 slots worth of data
n_ofdm_sym=6*10*2*n_symb_dl+2*n_symb_dl;
tfg=NaN(n_ofdm_sym,128);
tfg_timestame=NaN(1,n_ofdm_sym);
sym_num=0;
for t=1:n_ofdm_sym
  indices=round(dft_location):round(dft_location)+127;
  % Perform 2 sample coarse TOC
  %indices=[indices(3:end) indices(1:2)];
  %indices=[indices(end-2:end) indices(1:end-3)];
  tfg(t,:)=dft(capbuf(indices));
  tfg_timestamp(t)=dft_location;
  if (n_symb_dl==6)
    dft_location=dft_location+k_factor*(128+16);
  else
    if (sym_num==6)
      dft_location=dft_location+k_factor*(128+10);
    else
      dft_location=dft_location+k_factor*(128+9);
    end
    sym_num=mod(sym_num+1,7);
  end
end

% Extract the columns of interest
tfg=[tfg(:,end-35:end) tfg(:,2:37)];

% Compensate for the residual time offset.
cn=[-36:-1 1:36];
for t=1:n_ofdm_sym
  ideal_offset=tfg_timestamp(t);
  actual_offset=round(ideal_offset);
  % How late were we in locating the DFT
  late=actual_offset-ideal_offset;
  % Compensate for the improper location of the DFT
  tfg(t,:)=tfg(t,:).*exp(-j*2*pi*cn*late/128);
end
% At this point, tfg(t,:) contains the results of a DFT that was performed
% at time offset tfg_timestamp(t). Note that tfg_timestamp(t) is not an
% integer!

