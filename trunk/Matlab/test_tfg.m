% Create test vectors to test all the code starting with generation of
% the tfg.

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

load test_tfg.mat
peaks_in_n_id_1=peaks.n_id_1;
peaks_in_n_id_2=peaks.n_id_2;
peaks_in_cp_type=strcmpi(peaks.cp_type,'extended')+0;
peaks_in_frame_start=peaks.frame_start;
peaks_in_freq_fine=peaks.freq_fine;

% Create the time/frequency grid
[tfg tfg_timestamp]=extract_tfg(peaks,capbuf,fc);

% Perform FOE-C and TOE-C
[tfg_comp tfg_comp_timestamp peaks_out_tfoec]=tfoec(peaks,tfg,tfg_timestamp,fc);

% Decode the MIB!!! Finally!!!
peaks_out_mib=decode_mib(peaks_out_tfoec,tfg_comp);

peaks_out_freq_superfine=peaks_out_mib.freq_superfine;
peaks_out_n_rb_dl=peaks_out_mib.n_rb_dl;
peaks_out_phich_dur=strcmpi(peaks_out_mib.phich_dur,'extended')+0;
peaks_out_phich_res=peaks_out_mib.phich_res;
peaks_out_sfn=peaks_out_mib.sfn;

itsave('test_tfg.it',peaks_in_n_id_1,peaks_in_n_id_2,peaks_in_cp_type,peaks_in_frame_start,peaks_in_freq_fine,capbuf,fc,tfg,tfg_timestamp,tfg_comp,tfg_comp_timestamp,peaks_out_freq_superfine,peaks_out_n_rb_dl,peaks_out_phich_dur,peaks_out_phich_res,peaks_out_sfn);

