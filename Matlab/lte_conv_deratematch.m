function d_est=lte_conv_deratematch(e_est,n_c);

% d_est=lte_conv_deratematch(e_est,n_c);
%
% Undoes the rate-matching operation which mapped the 'd_est' bits to the
% 'e_est' bits.
%
% n_c is the number of bits that were originally input to the convolutional
% encoder.
%
% e_est(k) is the probability that a '1' was transmitted in bit position k.

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
error(chk_param(e_est,'e_est','real','vector','>=',0,'<=',1));

% Need to know the mapping of which d_est bits ended up where in the e_est
% sequence.
% (There are other ways to do this!)
d_probe=complex(repmat([1;2;3],1,n_c),repmat(1:n_c,3,1));
e_probe=lte_conv_ratematch(d_probe,length(e_est));

% Convert from probabilities to values on the x axis assuming the symbols are
% +1 and -1 and the received noise power is 1.
e_est(e_est<=eps)=eps;
e_est(e_est>=1-eps(1))=1-eps;
e_x=(log(e_est)-log(1-e_est))/2;

% Collect all the observations
d_x=zeros(3,n_c);
d_x_count=zeros(3,n_c);
for t=1:length(e_est)
  r=real(e_probe(t));
  c=imag(e_probe(t));
  d_x(r,c)=d_x(r,c)+e_x(t);
  d_x_count(r,c)=d_x_count(r,c)+1;
end
d_x=d_x(:);
d_x_count=d_x_count(:);
d_x(d_x_count~=0)=d_x(d_x_count~=0)./d_x_count(d_x_count~=0);
d_x=reshape(d_x,3,n_c);
d_x_count=reshape(d_x_count,3,n_c);
%d_x

% Convert back to probabilities
d_est=exp(-(d_x-1).^2/2)./(exp(-(d_x+1).^2/2)+exp(-(d_x-1).^2/2));

