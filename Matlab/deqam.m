function [bits const_prob]=deqam(syms,np,modulation,map)

% bits=deqam(syms,np,modulation,map)
%
% Convert between received symbols and bit probabilities.
%
% 'syms' contains the received symbols and it is assumed that channel
% compensation has already been performed. ie, in the absence of noise,
% the 'syms' vector is only composed of points on the constellation of
% the modulation map.
%
% 'np(k)' is the power of the gaussian noise assumed to be present in syms(k).
% If 'np' is a scalar, it is assumed that all symbols have the same amount
% of noise power.
%
% Descriptions of 'modulation' and 'map' can be found in the documentation
% for qam().

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

error(nargchk(3,4,nargin));
error(chk_param(syms,'syms','vector'));
error(chk_param(np,'np','vector','real','>=',0));
error(chk_param(modulation,'modulation','string'));
if (nargin<4)
  map='LTE';
end
if (~ischar(map))
  error(chk_param(map,'map','real','integer','>=',0));
end

n_sym=length(syms);

% Second order checking and derivations
if (length(np)==1)
  np=repmat(np,1,n_sym);
end
if (length(np)~=length(syms))
  error('np must either be of length 1 or equal to the length of syms');
end
switch (lower(modulation))
  case 'qam'
    bps=2;
  case 'qam16'
    bps=4;
  case 'qam64'
    bps=6;
  case 'qam256'
    bps=8;
  otherwise
    error('Unrecognized modulation specified.');
end
map_dim=2^(bps/2);
if (ischar(map))
  if (strcmpi(map,'LTE'))
    switch (bps)
      case 2
        map=[2 0 ; 3 1];
      case 4
        map=[11 9 1 3; 10 8 0 2; 14 12 4 6; 15 13 5 7];
      case 6
        map=[47 45 37 39 7 5 13 15; 46 44 36 38 6 4 12 14; 42 40 32 34 2 0 8 10; 43 41 33 35 3 1 9 11; 59 57 49 51 19 17 25 27; 58 56 48 50 18 16 24 26; 62 60 52 54 22 20 28 30; 63 61 53 55 23 21 29 31];
      otherwise
        error(['no lte mapping is known for ' modulation]);
    end
  else
    error('unrecognized mapping specified.');
  end
else
  error(chk_param(map,'map','integer','real','>=',0,'<=',2^bps-1));
end
if ((length(map(:))~=2^bps)||(size(map,1)~=map_dim)||size(map,2)~=map_dim)
  error(sprintf('map matrix dimensions must be %ix%i',map_dim,map_dim));
end
if (length(unique(map(:)))~=2^bps)
  error('map matrix contains duplicate entries');
end

% Create the mapping between integers and complex modulation symbols.
const=complex(repmat(1:map_dim,map_dim,1),repmat(transpose(map_dim:-1:1),1,map_dim));
const=const-mean(const(:));
const=const/sqrt(sigpower(const(:)));
map_flat=NaN(1,2^bps);
for t=1:map_dim
  for m=1:map_dim
    map_flat(map(t,m)+1)=const(t,m);
  end
end

% For each received symbol, calculate the probablity that the received symbol
% would be received if the actual transmitted symbol had been one of the
% constellation points.
const_prob=exp(absx2(repmat(map_flat,n_sym,1)-repmat(transpose(syms),1,2^bps))./repmat(transpose(-np),1,2^bps));

% to avoid NaN
for i=1: (2^bps)
    a = const_prob(:,i);
    a(a<1e-20) = 1e-20;
    const_prob(:,i) = a;
end

% Normalize so that each row sums to 1.
const_prob=const_prob./repmat(sum(const_prob,2),1,2^bps);

% Now, translate into bits. Find the probability that any individual bit
% is 1.
bit_prob=NaN(n_sym,bps);
for t=1:bps
  % Find the constellation points where this bit is 1.
  cols=find(bitget(0:2^bps-1,bps-t+1));
  bit_prob(:,t)=sum(const_prob(:,cols),2);
end

bits=reshape(transpose(bit_prob),1,n_sym*bps);

% Correct machine truncation errors
bits(bits<0)=0;
bits(bits>1)=1;

