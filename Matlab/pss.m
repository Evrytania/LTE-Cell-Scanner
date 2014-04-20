function y=pss(n_id_2)

% y=pss(n_id_2)
%
% Return the PSS signal indexed by n_id_2.
%
% y is a length 62 vector containing the symbols that should be transmitted
% on the subcarriers occupied by the PSS. The function calling this function
% must insert the DC subcarrier and also the 5 subcarriers required for padding
% on the left and right side of the PSS.

% Copyright 2012 Evrytania LLC (http://www.evrytania.com)
%
% Written by James Peroulas <james@evrytania.com>
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

error(nargchk(1,1,nargin));
error(chk_param(n_id_2,'n_id_2','scalar','real','integer','>=',0,'<=',2));

lut=[25 29 34];
y=zadoff_chu(63,lut(n_id_2+1));
y(32)=[];

