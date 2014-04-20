function y=fshift(x,f,fs)

% y=fshift(x,f,fs)
%
% Shift time domain signal x, sampled at a sampling rate of fs by
% f. fs defaults to 2.

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

error(nargchk(2,3,nargin));
error(chk_param(x,'x','vector','horizontal'));
error(chk_param(f,'f','scalar','real'));
if (nargin<3)
  fs=2;
end
error(chk_param(fs,'fs','scalar','real','>',0));

y=x.*exp((j*pi*f/(fs/2))*(0:length(x)-1));

