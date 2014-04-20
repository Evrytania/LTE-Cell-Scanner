function e=lte_conv_ratematch(d,n_e);

% e=lte_conv_ratematch(d,n_e);
%
% Perform LTE ratematching for convolutionally encoded bits d to fit
% an e vector of length n_e.

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
error(chk_param(d,'d','integer','rows==',3));
error(chk_param(n_e,'n_e','scalar','integer','real','>=',1));

% Create the matrix before permutation.
n_c=32;
n_r=ceil(size(d,2)/n_c);

% Subblock interleaving
perm_pattern=1+[1 17 9 25 5 21 13 29 3 19 11 27 7 23 15 31 0 16 8 24 4 20 12 28 2 18 10 26 6 22 14 30];
v=NaN(3,n_r*n_c);
for t=1:3
  y=[NaN(1,n_r*n_c-size(d,2)) d(t,:)];
  y=transpose(reshape(y,n_c,n_r));

  % Permute the columns
  y_perm=y(:,perm_pattern);
  v(t,:)=y_perm(:);
end

% Bit collection
w=reshape(transpose(v),n_r*n_c*3,1);

% Selection
e=NaN(1,n_e);
k=1;
j=1;
while (k<=n_e)
  if (isfinite(w(j)))
    e(k)=w(j);
    k=k+1;
  end
  j=mod(j,3*n_r*n_c)+1;
end

