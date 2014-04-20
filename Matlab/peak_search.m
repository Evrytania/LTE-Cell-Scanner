function peaks=peak_search(xc_incoherent_collapsed_pow,xc_incoherent_collapsed_frq,Z_th1,f_search_set,fc,sampling_carrier_twist,k_factor)

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

% The cross-correlation properties of the PSS signals are not very good.
% If PSS0 was transmitted at offset I and was detected with a power
% level of 0dB, this same received signal will correlate with PSS2 with
% a power level of -8.2dB.
%pss_xc_matrix=udb10([NaN -14 -8.2;-14 NaN -11.7;-8.2 -11.7 NaN]+4);
pss_xc_matrix=udb10([NaN -8 -8; -8 NaN -8; -8 -8 NaN]);

xc_incoherent_working=xc_incoherent_collapsed_pow;
peaks=[];
while (1)
  % Search for the largest peak (not the largest peak relative to Z_th1!)
  [peak_pow peak_ind]=max(transpose(xc_incoherent_working));
  [peak_pow peak_n_id_2]=max(peak_pow);
  peak_n_id_2=peak_n_id_2-1;
  peak_ind=peak_ind(peak_n_id_2+1);
  if (peak_pow<Z_th1(peak_ind))
    break;
  end
  % Record this peak for further SSS processing
  rec.pow=peak_pow;
  rec.ind=peak_ind;
  rec.freq=f_search_set(xc_incoherent_collapsed_frq(peak_n_id_2+1,peak_ind));
  rec.n_id_2=peak_n_id_2;
  % The following are filled in after the sss search
  rec.n_id_1=NaN;
  rec.cp_type='';
  rec.frame_start=NaN;
  rec.freq_fine=NaN;
  % The following are filled in after time / frequency estimation and
  % correction
  rec.freq_superfine=NaN;
  % The following are filled in after successful MIB decoding
  rec.n_rb_dl=NaN;
  rec.phich_dur='';
  rec.phich_res=NaN;
  rec.sfn=NaN;
  rec.n_id_cell=NaN;
  rec.n_ports=NaN;
  rec.duplex_mode=NaN;
  
  if sampling_carrier_twist
    rec.k_factor = (fc-rec.freq)/fc;
  else
    rec.k_factor = k_factor;
  end
  
  peaks=[peaks rec];
  % Cancel out certain peaks around this one.
  % It's assumed that if there is a peak for PSS N at offset I, there cannot
  % be any more PSS N sequences within 2*137 samples of I.
  cancel=wrap(peak_ind-137*2:peak_ind+137*2,1,9601);
  xc_incoherent_working(peak_n_id_2+1,cancel)=0;
  % If PSS N has a peak at offset I, the other PSS sequences near I are only
  % considered real sequences if their power is above a certain threshold
  % relative to the power of PSS sequence N.
  for t=setdiff(0:2,peak_n_id_2)
    xc_incoherent_working(t+1,cancel(xc_incoherent_working(t+1,cancel)<peak_pow*pss_xc_matrix(peak_n_id_2+1,t+1)))=0;
  end
  % Because of the repetitive nature of the CRS, a PSS at offset I with power
  % P will produce correlation peaks during all CRS OFDM symbols with power
  % P-14dB. Cancel them out.
  xc_incoherent_working(xc_incoherent_working<peak_pow*udb10(-12))=0;
end

