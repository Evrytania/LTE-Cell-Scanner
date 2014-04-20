function peak_out = sss_detect(peak,capbuf,thresh2_n_sigma,fc,sampling_carrier_twist,tdd_flag)

% Perform maximum likelihood estimation of the SSS.

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

peak_loc=peak.ind;
peak_freq=peak.freq;
n_id_2_est=peak.n_id_2;

% % fc*k_factor is the receiver's actual RX center frequency.
if sampling_carrier_twist==1
    k_factor=(fc-peak.freq)/fc;
% else
%     k_factor=1;
else
    k_factor = peak.k_factor;
end

if tdd_flag == 1
    min_idx = 3*(128+32)+32;
    sss_ext_offset = 3*(128+32);
    sss_nrm_offset = 412;
else
    min_idx = 163-9;
    sss_ext_offset = 128+32;
    sss_nrm_offset = 128+9;
end
if (peak_loc<min_idx) % TDD
  peak_loc=peak_loc+9600*k_factor;
end
pss_loc_set=peak_loc:9600*k_factor:length(capbuf)-125-9;
% pss_loc_set=peak_loc + (0:9600:7*9600);
% pss_loc_set=peak_loc + (8*9600:9600:15*9600);
n_pss=length(pss_loc_set);
pss_np=NaN(1,n_pss);
h_raw=NaN(n_pss,62);
h_sm=NaN(n_pss,62);
sss_nrm_raw=NaN(n_pss,62);
sss_ext_raw=NaN(n_pss,62);
for k=1:n_pss
  pss_loc=round(pss_loc_set(k));

  % Find where to take the pss dft
  pss_dft_location=pss_loc+9-2;
  %if (pss_dft_location+127>length(capbuf))
  %  break;
  %end

  % Calculate the channel response
  dft_in=fshift(capbuf(pss_dft_location:pss_dft_location+127),-peak_freq,fs_lte/16);
  % TOC
  dft_in=[dft_in(3:end) dft_in(1:2)];
  dft_out=dft(dft_in);
  h_raw(k,:)=[dft_out(end-30:end) dft_out(2:32)];
  h_raw(k,:)=h_raw(k,:).*conj(pss(n_id_2_est));
  %plot(angle(h_raw(k,:)));
  %ylim([-pi pi]);
  %drawnow;
  %pause

  % Smoothening... Basic...
  for t=1:62
    %arm_length=min([6 t-1 62-t]);
    lt=max([1 t-6]);
    rt=min([62 t+6]);
    % Growing matrix...
    %h_sm(k,t)=mean(h_raw(k,t-arm_length:t+arm_length));
    h_sm(k,t)=mean(h_raw(k,lt:rt));
  end

  % Estimate noise power.
  pss_np(k)=sigpower(h_sm(k,:)-h_raw(k,:));

  % Calculate the SSS in the frequency domain (ext)
  sss_ext_dft_location=pss_dft_location-sss_ext_offset;
  dft_in=fshift(capbuf(sss_ext_dft_location:sss_ext_dft_location+127),-peak_freq,fs_lte/16);
  % TOC
  dft_in=[dft_in(3:end) dft_in(1:2)];
  dft_out=dft(dft_in);
  sss_ext_raw(k,1:62)=[dft_out(end-30:end) dft_out(2:32)];

  % Calculate the SSS in the frequency domain (nrm)
  sss_nrm_dft_location=pss_dft_location-sss_nrm_offset;
  dft_in=fshift(capbuf(sss_nrm_dft_location:sss_nrm_dft_location+127),-peak_freq,fs_lte/16);
  % TOC
  dft_in=[dft_in(3:end) dft_in(1:2)];
  dft_out=dft(dft_in);
  sss_nrm_raw(k,1:62)=[dft_out(end-30:end) dft_out(2:32)];
end

% % interpolation along time to get accurate response at sss.
% h_sm_ext_interp = zeros(n_pss, 62);
% h_sm_nrm_interp = zeros(n_pss, 62);
% for t=1:62
%     h_sm_ext_interp(:,t) = interp1(pss_loc_set, h_sm(:,t), pss_loc_set-sss_ext_offset, 'linear','extrap');
%     h_sm_nrm_interp(:,t) = interp1(pss_loc_set, h_sm(:,t), pss_loc_set-sss_nrm_offset, 'linear','extrap');
% end
% 
% h_raw_ext_interp = zeros(n_pss, 62);
% h_raw_nrm_interp = zeros(n_pss, 62);
% for t=1:62
%     h_raw_ext_interp(:,t) = interp1(pss_loc_set, h_raw(:,t), pss_loc_set-sss_ext_offset, 'linear','extrap');
%     h_raw_nrm_interp(:,t) = interp1(pss_loc_set, h_raw(:,t), pss_loc_set-sss_nrm_offset, 'linear','extrap');
% end
% 
% pss_np_ext=zeros(1,n_pss);
% pss_np_nrm=zeros(1,n_pss);
% for k=1:n_pss
%     pss_np_ext(k)=sigpower(h_sm_ext_interp(k,:)-h_raw_ext_interp(k,:));
%     pss_np_nrm(k)=sigpower(h_sm_nrm_interp(k,:)-h_raw_nrm_interp(k,:));
% end

% ----recover original algorithm by using following 4 lines
h_sm_ext_interp = h_sm;
h_sm_nrm_interp = h_sm;
pss_np_ext = pss_np;
pss_np_nrm = pss_np;

% Combine results from different slots
sss_h1_nrm_np_est=NaN(1,62);
sss_h2_nrm_np_est=NaN(1,62);
sss_h1_ext_np_est=NaN(1,62);
sss_h2_ext_np_est=NaN(1,62);

sss_h1_nrm_est=NaN(1,62);
sss_h2_nrm_est=NaN(1,62);
sss_h1_ext_est=NaN(1,62);
sss_h2_ext_est=NaN(1,62);
for t=1:62
  sss_h1_nrm_np_est(t)=real((1+ctranspose(h_sm_nrm_interp(1:2:end,t))*diag(1./pss_np_nrm(1:2:end))*h_sm_nrm_interp(1:2:end,t))^-1);
  sss_h2_nrm_np_est(t)=real((1+ctranspose(h_sm_nrm_interp(2:2:end,t))*diag(1./pss_np_nrm(2:2:end))*h_sm_nrm_interp(2:2:end,t))^-1);

  sss_h1_ext_np_est(t)=real((1+ctranspose(h_sm_ext_interp(1:2:end,t))*diag(1./pss_np_ext(1:2:end))*h_sm_ext_interp(1:2:end,t))^-1);
  sss_h2_ext_np_est(t)=real((1+ctranspose(h_sm_ext_interp(2:2:end,t))*diag(1./pss_np_ext(2:2:end))*h_sm_ext_interp(2:2:end,t))^-1);

  sss_h1_nrm_est(t)=sss_h1_nrm_np_est(t)*ctranspose(h_sm_nrm_interp(1:2:end,t))*diag(1./pss_np_nrm(1:2:end))*sss_nrm_raw(1:2:end,t);
  sss_h2_nrm_est(t)=sss_h2_nrm_np_est(t)*ctranspose(h_sm_nrm_interp(2:2:end,t))*diag(1./pss_np_nrm(2:2:end))*sss_nrm_raw(2:2:end,t);

  sss_h1_ext_est(t)=sss_h1_ext_np_est(t)*ctranspose(h_sm_ext_interp(1:2:end,t))*diag(1./pss_np_ext(1:2:end))*sss_ext_raw(1:2:end,t);
  sss_h2_ext_est(t)=sss_h2_ext_np_est(t)*ctranspose(h_sm_ext_interp(2:2:end,t))*diag(1./pss_np_ext(2:2:end))*sss_ext_raw(2:2:end,t);
end

% Maximum likelihood detection of SSS
log_lik_nrm=NaN(168,2);
log_lik_ext=NaN(168,2);
for t=0:167
  sss_h1_try=sss(t,n_id_2_est,0);
  sss_h2_try=sss(t,n_id_2_est,10);

  % Rotate the candiate sequence to match the received sequence.
  ang=angle(sum(conj([sss_h1_nrm_est sss_h2_nrm_est]).*[sss_h1_try sss_h2_try]));
  sss_h1_try=sss_h1_try*exp(j*-ang);
  sss_h2_try=sss_h2_try*exp(j*-ang);
  df=[sss_h1_try sss_h2_try]-[sss_h1_nrm_est sss_h2_nrm_est];
  log_lik_nrm(t+1,1)=sum(-[real(df) imag(df)].^2./repmat([sss_h1_nrm_np_est sss_h2_nrm_np_est],1,2));

  % Exchange h1 and h2 and re-do
  temp=sss_h1_try;
  sss_h1_try=sss_h2_try;
  sss_h2_try=temp;
  ang=angle(sum(conj([sss_h1_nrm_est sss_h2_nrm_est]).*[sss_h1_try sss_h2_try]));
  sss_h1_try=sss_h1_try*exp(j*-ang);
  sss_h2_try=sss_h2_try*exp(j*-ang);
  df=[sss_h1_try sss_h2_try]-[sss_h1_nrm_est sss_h2_nrm_est];
  log_lik_nrm(t+1,2)=sum(-[real(df) imag(df)].^2./repmat([sss_h1_nrm_np_est sss_h2_nrm_np_est],1,2));

  % Re-do for extended prefix
  % Rotate the candiate sequence to match the received sequence.
  temp=sss_h1_try;
  sss_h1_try=sss_h2_try;
  sss_h2_try=temp;
  ang=angle(sum(conj([sss_h1_ext_est sss_h2_ext_est]).*[sss_h1_try sss_h2_try]));
  sss_h1_try=sss_h1_try*exp(j*-ang);
  sss_h2_try=sss_h2_try*exp(j*-ang);
  df=[sss_h1_try sss_h2_try]-[sss_h1_ext_est sss_h2_ext_est];
  log_lik_ext(t+1,1)=sum(-[real(df) imag(df)].^2./repmat([sss_h1_ext_np_est sss_h2_ext_np_est],1,2));

  % Exchange h1 and h2 and re-do
  temp=sss_h1_try;
  sss_h1_try=sss_h2_try;
  sss_h2_try=temp;
  ang=angle(sum(conj([sss_h1_ext_est sss_h2_ext_est]).*[sss_h1_try sss_h2_try]));
  sss_h1_try=sss_h1_try*exp(j*-ang);
  sss_h2_try=sss_h2_try*exp(j*-ang);
  df=[sss_h1_try sss_h2_try]-[sss_h1_ext_est sss_h2_ext_est];
  log_lik_ext(t+1,2)=sum(-[real(df) imag(df)].^2./repmat([sss_h1_ext_np_est sss_h2_ext_np_est],1,2));
end

%warning('Check code here!!!!');
cp_type_flag = 0;
if (max(log_lik_nrm(:))>max(log_lik_ext(:)))
  cp_type_est='normal';
  cp_type_flag = 0;
  log_lik=log_lik_nrm;
else
  cp_type_est='extended';
  cp_type_flag = 1;
  log_lik=log_lik_ext;
end
% frame_start is the location of the 'start' of the cp of the frame.
% The first DFT for the frame should be located at frame_start+cp_length
if tdd_flag==1
    if cp_type_flag == 0
        frame_start=peak_loc+(-(2*(128+9)+1)-1920-2)*k_factor; % TDD NORMAL CP
    else
        frame_start=peak_loc+(-(2*(128+32))-1920-2)*k_factor; % TDD EXTENDED CP
    end
else
    frame_start=peak_loc+(128+9-960-2)*k_factor; % FDD
end
if (max(log_lik(:,1))>max(log_lik(:,2)))
  ll=log_lik(:,1);
else
  frame_start=frame_start+9600*k_factor;
  ll=log_lik(:,2);
end
frame_start=wrap(frame_start,0.5,2*9600+.5);
[lik_final n_id_1_est]=max(ll);
n_id_1_est=n_id_1_est-1;

% Second plot and second threshold check
L=[log_lik_nrm log_lik_ext];
L_mean=mean(L(:));
L_var=var(L(:));
figure(2);
plot(0:167,[log_lik_nrm log_lik_ext],[0 167],repmat(L_mean,1,2),[0 167],repmat(L_mean+sqrt(L_var)*thresh2_n_sigma,1,2));
zgo;
drawnow;
peak_out=peak;
if (lik_final<L_mean+sqrt(L_var)*thresh2_n_sigma)
  %disp('Thresh2 fail');
  peak_out.n_id_1=NaN;
  peak_out.cp_type='';
  peak_out.frame_start=NaN;
  peak_out.duplex_mode=NaN;
else
  peak_out.n_id_1=n_id_1_est;
  peak_out.cp_type=cp_type_est;
  peak_out.frame_start=frame_start;
  peak_out.duplex_mode=tdd_flag;
end

