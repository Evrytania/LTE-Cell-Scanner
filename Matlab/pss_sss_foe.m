function peak_out=pss_sss_foe(peak,capbuf,fc,sampling_carrier_twist,tdd_flag)

% Use (only) the PSS and SSS to calculate the frequency offset.

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

% fc*k_factor is the receiver's actual RX center frequency.
if sampling_carrier_twist==1
    k_factor=(fc-peak.freq)/fc;
% else
%     k_factor=1;
else
    k_factor = peak.k_factor;
end

%%%%%%
% Calculate the frequency offset
%%%%%%
% Find the location of the first SSS where we can take a DFT.
if (strcmpi(peak.cp_type,'normal'))
    if tdd_flag==1
        pss_sss_dist=round((3*(128+9)+1)*k_factor); % TDD
        first_sss_dft_location=peak.frame_start+(1920-128)*k_factor; % TDD
    else
        pss_sss_dist=round((128+9)*k_factor); % FDD
        first_sss_dft_location=peak.frame_start+(960-128-9-128)*k_factor; % FDD
    end
elseif (strcmpi(peak.cp_type,'extended'))
    if tdd_flag==1
        pss_sss_dist=round(3*(128+32)*k_factor); %TDD
        first_sss_dft_location=peak.frame_start+(1920-128)*k_factor; %TDD
    else
        pss_sss_dist=round((128+32)*k_factor); %FDD
        first_sss_dft_location=peak.frame_start+(960-128-32-128)*k_factor; %FDD
    end
else
  error('Check code...');
end
first_sss_dft_location=wrap(first_sss_dft_location,0.5,9600*2+0.5);
if (first_sss_dft_location-9600*k_factor>0.5)
%if (first_sss_dft_location>=9601)
  first_sss_dft_location=first_sss_dft_location-9600*k_factor;
  sn=10;
else
  sn=0;
end
sss_dft_loc_set=first_sss_dft_location:9600*k_factor:length(capbuf)-127-pss_sss_dist-100;
n_sss=length(sss_dft_loc_set);

h_raw_fo_pss=NaN(n_sss,62);
sss_raw_fo=NaN(n_sss,62);
%M1=0;
%M2=0;
sn=(1-(sn/10))*10;
M=0;
for k=1:n_sss
  sn=(1-(sn/10))*10;
  sss_dft_location=round(sss_dft_loc_set(k));

  % Find where to take the pss dft
  pss_dft_location=sss_dft_location+pss_sss_dist;
  %if (pss_dft_location+127>length(capbuf))
  %  break;
  %end

  % Find the PSS and use it to estimate the channel.
  dft_in_pss=fshift(capbuf(pss_dft_location:pss_dft_location+127),-peak.freq,fs_lte/16);
  % TOC
  %disp('toc diabled');
  dft_in_pss=[dft_in_pss(3:end) dft_in_pss(1:2)];
  dft_out_pss=dft(dft_in_pss);
  h_raw_fo_pss(k,:)=[dft_out_pss(end-30:end) dft_out_pss(2:32)];
  h_raw_fo_pss(k,:)=h_raw_fo_pss(k,:).*conj(pss(peak.n_id_2));

  % Smoothening... Basic...
  for t=1:62
    %arm_length=min([6 t-1 62-t]);
    lt=max([1 t-6]);
    rt=min([62 t+6]);
    % Growing matrix...
    %h_sm(k,t)=mean(h_raw_fo_pss(k,t-arm_length:t+arm_length));
    h_sm(k,t)=mean(h_raw_fo_pss(k,lt:rt));
  end
  %h_sm(k,:)=1;

  % Estimate noise power.
  % Growing matrix...
  pss_np(k)=sigpower(h_sm(k,:)-h_raw_fo_pss(k,:));

  %M2=M2+sum(conj(h_sm(k,:)).*h_raw_fo_pss(k,:))/pss_np(k);

  % Calculate the SSS in the frequency domain
  dft_in_sss=fshift(capbuf(sss_dft_location:sss_dft_location+127),-peak.freq,fs_lte/16)*exp(j*pi*-peak.freq/(fs_lte/16/2)*-pss_sss_dist);
  % TOC
  %disp('toc diabled');
  dft_in_sss=[dft_in_sss(3:end) dft_in_sss(1:2)];
  dft_out_sss=dft(dft_in_sss);
  sss_raw_fo(k,:)=[dft_out_sss(end-30:end) dft_out_sss(2:32)];
  sss_raw_fo(k,:)=sss_raw_fo(k,:).*conj(sss(peak.n_id_1,peak.n_id_2,sn));

  %M1=M1+sum(conj(h_sm(k,:)).*sss_raw_fo(k,:))/pss_np(k);
  M=M+sum(conj(sss_raw_fo(k,:)).*h_raw_fo_pss(k,:).*absx2(h_sm(k,:))./(2*absx2(h_sm(k,:))*pss_np(k)+pss_np(k)^2));
end
%f_off_est=peak.freq+angle(conj(M1)*M2)/(2*pi)/(1/(fs_lte/16)*pss_sss_dist);
f_off_est=peak.freq+angle(M)/(2*pi)/(1/(fs_lte/16)*pss_sss_dist);

peak_out=peak;
peak_out.freq_fine=f_off_est;

