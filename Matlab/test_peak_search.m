% Create test vectors to test the C++ implementation of peak_search.

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

% Load input data
load test_peak_search.mat

% Call function
peaks=peak_search(xc_incoherent_collapsed_pow,xc_incoherent_collapsed_frq,Z_th1,f_search_set);

peaks_pow=[peaks(:).pow];
peaks_ind=[peaks(:).ind];
peaks_freq=[peaks(:).freq];
peaks_n_id_2=[peaks(:).n_id_2];

itsave('test_peak_search.it',xc_incoherent_collapsed_pow,xc_incoherent_collapsed_frq,Z_th1,f_search_set,peaks_pow,peaks_ind,peaks_freq,peaks_n_id_2);

