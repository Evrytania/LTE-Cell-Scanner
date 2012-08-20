% Create test vectors to test the C++ implementation of xcorr_pss.

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

% Load capture buffer
load test_xcorr_pss.mat

ds_comb_arm=2;
fc=739e6;
f_search_set=35e3:5e3:45e3;

% Call xcorr_pss
[xc_incoherent_collapsed_pow xc_incoherent_collapsed_frq n_comb_xc n_comb_sp xc_incoherent_single xc_incoherent sp_incoherent xc sp]=xcorr_pss(capbuf,f_search_set,ds_comb_arm,fc);

capbuf=capbuf(:);
xc_incoherent_collapsed_pow=xc_incoherent_collapsed_pow(:);
xc_incoherent_collapsed_frq=xc_incoherent_collapsed_frq(:);
xc_incoherent_single=xc_incoherent_single(:);
xc_incoherent=xc_incoherent(:);
sp_incoherent=sp_incoherent(:);
xc=xc(:);
sp=sp(:);

itsave('test_xcorr_pss.it',capbuf,f_search_set,ds_comb_arm,fc,xc_incoherent_collapsed_pow,xc_incoherent_collapsed_frq,n_comb_xc,n_comb_sp,xc_incoherent_single,xc_incoherent,sp_incoherent,xc,sp);

