function [xc_incoherent_collapsed_pow xc_incoherent_collapsed_frq n_comb_xc n_comb_sp xc_incoherent_single xc_incoherent sp_incoherent sp] = ...
    xcorr_pss(capbuf,f_search_set,ds_comb_arm,fc,sampling_carrier_twist,k_factor,xc)

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

% Perform the main correlations necessary to detect the PSS

error(nargchk(7,7,nargin));
error(chk_param(capbuf,'capbuf','vector','horizontal'));
error(chk_param(f_search_set,'f_search_set','vector','real'));
error(chk_param(ds_comb_arm,'ds_comb_arm','scalar','real','integer','>=',0));
error(chk_param(fc,'fc','scalar','real','>',0));

% % Create the time domain pss signals
% persistent pss_td
% persistent pss_td_populated
% if (isempty(pss_td_populated))
%   pss_td=NaN(3,128+9);
%   for t=0:2
%     temp=pss(t);
%     %temp=temp(6:67);
%     temp=[0 temp(32:end) zeros(1,65) temp(1:31)];
%     temp_td=idft(temp)*sqrt(128/62);
%     pss_td(t+1,:)=[temp_td(end-8:end) temp_td];
%   end
%   pss_td_populated=1;
% end

n_cap=length(capbuf);
n_f=length(f_search_set);

% % Correlate against the PSS
% xc=NaN(3,n_cap-136,n_f);
% for foi=1:n_f
% %   f_off=f_search_set(foi);
%   for t=1:3
% %     temp=conj(fshift(pss_td(t,:),f_off,fs_lte/16)/137);
% %     for k=1:n_cap-136
% %       xc(t,k,foi)=sum(temp.*capbuf(k:k+136));
% %     end
%     col_idx = (t-1)*n_f + foi;
%     xc(t,:,foi)=corr_store(1:(n_cap-136),col_idx);
%   end
% end

% Calculate the received power in the vicinity of the possible PSS
sp=NaN(1,n_cap-136-137);
capbufx2=absx2(capbuf);
sp(1)=sum(capbufx2(1:137*2));
for k=2:n_cap-136-137
  sp(k)=sp(k-1)-capbufx2(k-1)+capbufx2(k+137*2-1);
end
sp=sp/(137*2);
assert(any(isnan(sp))==0);

% Perform incoherent combining
xc_incoherent=NaN(3,9600,n_f);
xc_incoherent_single=NaN(3,9600,n_f);
n_comb_xc=floor((size(xc,2)-100)/9600);
for foi=1:n_f
  
  if sampling_carrier_twist == 1
      f_off=f_search_set(foi);
      % fc*k_factor is the receiver's actual RX center frequency.
      k_factor=(fc-f_off)/fc;
%   else
%       k_factor = 1; % because it is already corrected outside
  end
  
  for t=1:3
    %xc_incoherent_single(t,:,foi)=sum(transpose(reshape(absx2(xc(t,1:n_comb_xc*9600,foi)),9600,n_comb_xc)),1)/n_comb_xc;
    %xc_incoherent_single(t,:,foi)=absx2(xc(t,1:9600,foi));
    xc_incoherent_single(t,:,foi)=zeros(1,9600);
    for m=1:n_comb_xc
      actual_time_offset=(m-1)*.005*k_factor;
      actual_start_index=round(actual_time_offset*fs_lte/16)+1;
%       xc_incoherent_single(t,:,foi)=xc_incoherent_single(t,:,foi)+absx2(xc(t,actual_start_index:actual_start_index+9599,foi));
      xc_incoherent_single(t,:,foi)=xc_incoherent_single(t,:,foi)+xc(t,actual_start_index:actual_start_index+9599,foi);
    end
  end
  xc_incoherent_single(:,:,foi)=xc_incoherent_single(:,:,foi)/n_comb_xc;
  % Combine adjacent samples that might be different taps from the same channel.
  xc_incoherent(:,:,foi)=xc_incoherent_single(:,:,foi);
  for t=1:ds_comb_arm
    for k=1:3
      xc_incoherent(k,:,foi)=xc_incoherent(k,:,foi)+tshift(xc_incoherent_single(k,:,foi),-t)+tshift(xc_incoherent_single(k,:,foi),t);
    end
  end
  xc_incoherent(:,:,foi)=xc_incoherent(:,:,foi)/(2*ds_comb_arm+1);
end
sp_incoherent=NaN(1,9600);
n_comb_sp=floor(length(sp)/9600);
sp_incoherent=sum(transpose(reshape(sp(1:n_comb_sp*9600),9600,n_comb_sp)),1)/n_comb_sp;

% Align the correlations and the signal power measurements.
% xc(:,1) represents the correlation from capture buffer offset 1 to 137.
% sp(1) represents the signal power for samples 1 through 274.
% Suppose capbuf(1) contains the first sample of an SSS and thus the PSS is
% located from samples 138 to 274. The pss correlation peak will be located
% at time sample 138. This PSS correlation peak must be compared against
% the signal power for samples 1 through 274. Thus, the signal power must
% be cyclically shifted right by 137 samples.
sp_incoherent=tshift(sp_incoherent,137);

% For each time offset and PSS index, find the peak correlation
% among all the frequencies.
xc_incoherent_collapsed_pow=NaN(3,9600);
xc_incoherent_collapsed_frq=NaN(3,9600);
for t=1:3
  [pow frq]=max(transpose(shiftdim(xc_incoherent(t,:,:),1)),[],1);
  xc_incoherent_collapsed_pow(t,:)=pow;
  xc_incoherent_collapsed_frq(t,:)=frq;
end

