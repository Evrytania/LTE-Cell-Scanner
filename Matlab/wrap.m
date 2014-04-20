function y=wrap(x,sm,lg)

% y=wrap(x,sm,lg)
%
% Wrap the values of x so that they are >=sm and <lg

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

error(nargchk(3,3,nargin));
error(chk_param(x,'x','numeric','real'));
error(chk_param(sm,'sm','scalar','real'));
error(chk_param(sm,'lg','scalar','real'));
if (lg<=sm)
  error('lg must be > sm');
end

y=mod(x-sm,lg-sm)+sm;

