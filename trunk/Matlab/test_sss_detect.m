% This script generates test data for sss_detect and pss_sss_foe

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

load test_sss_detect.mat

clear peaks_out;
sss_h1_np_est=NaN(1,62);
sss_h2_np_est=NaN(1,62);
sss_h1_nrm_est=NaN(1,62);
sss_h2_nrm_est=NaN(1,62);
sss_h1_ext_est=NaN(1,62);
sss_h2_ext_est=NaN(1,62);
for t=1:length(peaks)
  [peaks_out_temp sss_h1_np_est_temp sss_h2_np_est_temp sss_h1_nrm_est_temp sss_h2_nrm_est_temp sss_h1_ext_est_temp sss_h2_ext_est_temp]=sss_detect(peaks(t),capbuf,thresh2_n_sigma,fc);
  if (isfinite(peaks_out_temp.n_id_1))
    peaks_out_temp=pss_sss_foe(peaks_out_temp,capbuf,fc);
  end
  peaks_out(t)=peaks_out_temp;
  sss_h1_np_est(t,:)=sss_h1_np_est_temp;
  sss_h2_np_est(t,:)=sss_h2_np_est_temp;
  sss_h1_nrm_est(t,:)=sss_h1_nrm_est_temp;
  sss_h2_nrm_est(t,:)=sss_h2_nrm_est_temp;
  sss_h1_ext_est(t,:)=sss_h1_ext_est_temp;
  sss_h2_ext_est(t,:)=sss_h2_ext_est_temp;
end
peaks_pow=[peaks(:).pow];
peaks_ind=[peaks(:).ind];
peaks_freq=[peaks(:).freq];
peaks_n_id_2=[peaks(:).n_id_2];

peaks_out_n_id_1=[peaks_out(:).n_id_1];
peaks_out_frame_start=[peaks_out(:).frame_start];
clear peaks_out_cp_type;
for t=1:length(peaks_out);
  peaks_out_cp_type(t)=strcmpi('extended',peaks_out(t).cp_type)+0;
end

peaks_out_freq_fine=[peaks_out(:).freq_fine];

itsave('test_sss_detect.it',peaks_pow,peaks_ind,peaks_freq,peaks_n_id_2,capbuf,thresh2_n_sigma,fc,peaks_out_n_id_1,peaks_out_frame_start,peaks_out_cp_type,peaks_out_freq_fine,sss_h1_np_est,sss_h2_np_est,sss_h1_nrm_est,sss_h2_nrm_est,sss_h1_ext_est,sss_h2_ext_est);

