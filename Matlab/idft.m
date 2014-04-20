function y=idft(p1,p2,p3,p4)

% y=idft(x,...);
%
% Exactly the same as the ifft function except that the output of this function
% is scaled so that sigpower(x) is the same as sigpower(dft(x));

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

error(nargchk(1,4,nargin));

scale=sqrt(length(p1));
switch nargin
  case 1
    y=ifft(p1)*scale;
  case 2
    y=ifft(p1,p2)*scale;
  case 3
    y=ifft(p1,p2,p3)*scale;
  case 4
    y=ifft(p1,p2,p3,p4)*scale;
  otherwise
    error('Check code...');
end

