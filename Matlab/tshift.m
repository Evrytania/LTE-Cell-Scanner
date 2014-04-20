function y=tshift(x,t);

% y=tshift(x,t);
%
% Cyclically shift sequence x to the rigth by t samples.
% t is allowed to be a non-integer.

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

error(nargchk(2,2,nargin));
error(chk_param(x,'x','vector','horizontal'));
error(chk_param(t,'t','scalar','real'));

% Shortcut if possible
if (t==floor(t))
  t=mod(t,length(x));
  y(1+t:length(x))=x(1:length(x)-t);
  y(1:t)=x(end-t+1:end);
  return
end

n_samp=length(x);

x_f=fft(x);
x_f=fftshift(x_f);

f=linspace(0,2,n_samp+1);
f=f(1:end-1);
f(f>=1)=f(f>=1)-2;
f=fftshift(f);

phase_shift=exp((-j*t*pi)*f);
% Fix the phase compensation if x is even length
if (bitand(n_samp,1)==0)
  phase_shift(1)=real(phase_shift(1));
end

x_f=x_f.*phase_shift;
y=ifft(ifftshift(x_f));

if (isreal(x))
  y=real(y);
end

