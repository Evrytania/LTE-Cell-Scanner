function [tfg_comp, tfg_comp_timestamp, peak_out]=tfoec(peak,tfg,tfg_timestamp,fc,sampling_carrier_twist, nRB)

% add 100RB support. Jiao Xianjun(putaohsu@gmail.com)
% Compensates for frequency offset, time offset, and also rotates the
% RS so that they properly reflect the channel response.

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

% nRB = 6;
if nRB == 6
    decimation_ratio = 16;
elseif nRB == 100
    decimation_ratio = 1;
else
    disp('nRB must be 6 or 100!');
    return;
end

fs = fs_lte/decimation_ratio;
fft_size = 2048/decimation_ratio;
nSC = nRB*12;
nRS = nRB*2;

% fc*k_factor is the receiver's actual RX center frequency.
if sampling_carrier_twist==1
    k_factor=(fc-peak.freq_fine)/fc;
else
    k_factor = peak.k_factor;
end

% Second order derivations
n_ofdm=size(tfg,1);
if (strcmpi(cp_type,'normal'))
  n_symb_dl=7;
elseif (strcmpi(cp_type,'extended'))
  n_symb_dl=6;
else
 error('Check code...');
end
n_id_cell=n_id_2+3*n_id_1;

peak_out=peak;

% Pre-compute the RS
rs0_start=NaN(20,nRS);
rs0_mid=NaN(20,nRS);
for slot_num=0:19
  [r, shift_start]=rs_dl(slot_num,0,0,n_id_cell,nRB,cp_type);
  rs0_start(slot_num+1,:)=r;
  [r, shift_mid]=rs_dl(slot_num,n_symb_dl-3,0,n_id_cell,nRB,cp_type);
  rs0_mid(slot_num+1,:)=r;
end
%shift_start
%shift_mid

% Rotate the RS based on the known transmitted sequence.
rs_set=sort([1:n_symb_dl:n_ofdm n_symb_dl-3+1:n_symb_dl:n_ofdm]);
slot_num=0;
for t=1:length(rs_set)
  idx=rs_set(t);
  if (mod(t,2)==1)
    tfg(idx,shift_start+1:6:end)=tfg(idx,shift_start+1:6:end).*conj(rs0_start(slot_num+1,:));
  else
    tfg(idx,shift_mid+1:6:end)=tfg(idx,shift_mid+1:6:end).*conj(rs0_mid(slot_num+1,:));
  end

  if (mod(t,2)==0)
    slot_num=mod(slot_num+1,20);
  end
end

% FOE using starting OFDM symbol
foe=0;
for t=1:nRS
  idx=(t-1)*6+shift_start+1;
  % Extract all the RS for this subcarrier
  rs_extracted=transpose(tfg(1:n_symb_dl:end,idx));
  % Compensate for the known RS
  %rs_known=transpose(rs0_start(:,t));
  %rs_known=repmat(rs_known,1,ceil(length(rs_extracted)/20));
  %rs_known=rs_known(1:length(rs_extracted));
  %rs_comp=rs_extracted.*conj(rs_known);
  rs_comp=rs_extracted;
  %plot(unwrap(angle(rs_comp)),'x-');
  %drawnow
  %pause
  %keyboard
  % Calculate the frequency offset
  foe=foe+sum(conj(rs_comp(1:end-1)).*rs_comp(2:end));
end
%residual_f=angle(foe)/(2*pi)/(k_factor*.0005);
% FOE using mid OFDM symbol
for t=1:nRS
  idx=(t-1)*6+shift_mid+1;
  % Extract all the RS for this subcarrier
  rs_extracted=transpose(tfg(n_symb_dl-3+1:n_symb_dl:end,idx));
  % Compensate for the known RS
  %rs_known=transpose(rs0_mid(:,t));
  %rs_known=repmat(rs_known,1,ceil(length(rs_extracted)/20));
  %rs_known=rs_known(1:length(rs_extracted));
  %rs_comp=rs_extracted.*conj(rs_known);
  rs_comp=rs_extracted;
  %plot(unwrap(angle(rs_comp)),'x-');
  %drawnow
  %pause
  % Calculate the frequency offset
  foe=foe+sum(conj(rs_comp(1:end-1)).*rs_comp(2:end));
end

residual_f=angle(foe)/(2*pi)/(k_factor*.0005);
peak_out.freq_superfine=peak_out.freq_fine+residual_f;
if sampling_carrier_twist
    peak_out.k_factor=(fc-peak_out.freq_superfine)/fc;
end

% Perform FOC. This does not compensate for the ICI, only the bulk
% frequency offset and time shift between OFDM symbols.
if sampling_carrier_twist==1
    k_factor_residual=(fc-residual_f)/fc;
else
    k_factor_residual=1;
end

tfg_comp=NaN(size(tfg));
tfg_comp_timestamp=1+k_factor_residual*(tfg_timestamp-1);
cn=[-(nSC/2):-1 1:(nSC/2)];
for t=1:n_ofdm
  tfg_comp(t,:)=tfg(t,:)*exp(1i*2*pi*-residual_f*(tfg_comp_timestamp(t)-1)/fs);
  % How late has the DFT been placed?
  late=tfg_timestamp(t)-tfg_comp_timestamp(t);
  % Compensate for the location of the DFT
  tfg_comp(t,:)=tfg_comp(t,:).*exp(-1i*2*pi*cn*late/fft_size);
end

% Perform TOE
% We are assuming that the time offset is not changing for the entire
% duration of the capture buffer.
toe=0;
offsets=[shift_start+1 shift_mid+1];
rs_set=sort([1:n_symb_dl:n_ofdm n_symb_dl-3+1:n_symb_dl:n_ofdm]);
for t=1:length(rs_set)-1
  % We are crossing over the DC subcarrier so this is not 100% correct...
  if (offsets(2)>offsets(1))
    r1=rs_set(t);
    r2=rs_set(t+1);
    o1=offsets(1);
    o2=offsets(2);
  else
    r1=rs_set(t+1);
    r2=rs_set(t);
    o1=offsets(2);
    o2=offsets(1);
  end
  toe_current1=sum(conj(tfg_comp(r1,o1:6:end)).*tfg_comp(r2,o2:6:end));
  toe_current2=sum(conj(tfg_comp(r2,o2:6:end-6)).*tfg_comp(r1,o1+6:6:end));
  %1
  %-angle(toe_current1)/3/(2*pi/128)
  %-angle(toe_current2)/3/(2*pi/128)
  toe_current=toe_current1+toe_current2;
  %2
  %-angle(toe_current)/3/(2*pi/128)
  %3
  %keyboard
  toe=toe+toe_current;

  %figure(1);
  %plot(abs(tfg_comp(rs_set(t),offsets(1):6:end)));
  %ylim([0 1]);
  %figure(2);
  %plot(unwrap(angle(tfg_comp(rs_set(t),offsets(1):6:end))));
  %ylim([2*-pi .5*pi]);
  %drawnow;
  %pause
  %if (t==length(rs_set)-1)
  %  keyboard
  %end

  offsets=fliplr(offsets);
  %keyboard
end
delay=-angle(toe)/3/(2*pi/fft_size);
% disp(['residual_f ' num2str(residual_f)]);
% disp(['delay ' num2str(delay)]);

% Finally, perform TOC.
tfg_comp=tfg_comp.*repmat(exp((1i*2*pi/fft_size*delay)*cn),n_ofdm,1);

%offsets=[shift_start+1 shift_mid+1];
%for t=1:length(rs_set)-1
%  figure(2);
%  plot(unwrap(angle(tfg_comp(rs_set(t),offsets(1):6:end))));
%  %ylim([0 pi]);
%  drawnow;
%  pause
%  offsets=fliplr(offsets);
%end

% Rotate the RS back to where they were
rs_set=sort([1:n_symb_dl:n_ofdm n_symb_dl-3+1:n_symb_dl:n_ofdm]);
slot_num=0;
for t=1:length(rs_set)
  idx=rs_set(t);
  if (mod(t,2)==1)
    tfg_comp(idx,shift_start+1:6:end)=tfg_comp(idx,shift_start+1:6:end).*rs0_start(slot_num+1,:);
  else
    tfg_comp(idx,shift_mid+1:6:end)=tfg_comp(idx,shift_mid+1:6:end).*rs0_mid(slot_num+1,:);
  end

  if (mod(t,2)==0)
    slot_num=mod(slot_num+1,20);
  end
end

