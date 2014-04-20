function s=sss(n_id_1,n_id_2,slot_num);

% s=sss(n_id_1,n_id_2,slot_num);
%
% Return the sss for slot slot_num for the specified n_id_1 and n_id_2.
%
% s is of length 62 and only includes the non-zero subcarriers. The calling
% function must insert the DC subcarrier and the 5 zero subcarriers that
% surround the sss.

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
error(chk_param(n_id_1,'n_id_1','scalar','real','integer','>=',0,'<=',167));
error(chk_param(n_id_2,'n_id_2','scalar','real','integer','>=',0,'<=',2));
error(chk_param(slot_num,'slot_num','scalar','real','integer','>=',0,'<=',10));
if ((slot_num~=0)&&(slot_num~=10))
  error('slot_num must be either 0 or 10');
end

% Calculate m0 and m1 from n_id_1
qp=floor(n_id_1/30);
q=floor((n_id_1+qp*(qp+1)/2)/30);
mp=n_id_1+q*(q+1)/2;
m0=mod(mp,31);
m1=mod(m0+floor(mp/31)+1,31);

%s_td=[0 0 0 0 1];
%for t=1:26
%  s_td=[s_td mod(s_td(end-2)+s_td(end-4),2)];
%end
s_td=[0 0 0 0 1 0 0 1 0 1 1 0 0 1 1 1 1 1 0 0 0 1 1 0 1 1 1 0 1 0 1];
s_td=1-2*s_td;

%c_td=[0 0 0 0 1];
%for t=1:26
%  c_td=[c_td mod(c_td(end-1)+c_td(end-4),2)];
%end
c_td=[0 0 0 0 1 0 1 0 1 1 1 0 1 1 0 0 0 1 1 1 1 1 0 0 1 1 0 1 0 0 1];
c_td=1-2*c_td;

%z_td=[0 0 0 0 1];
%for t=1:26
%  z_td=[z_td mod(z_td(end)+z_td(end-2)+z_td(end-3)+z_td(end-4),2)];
%end
z_td=[0 0 0 0 1 1 1 0 0 1 1 0 1 1 1 1 1 0 1 0 0 0 1 0 0 1 0 1 0 1 1];
z_td=1-2*z_td;

s0_m0=s_td(mod(m0:30+m0,31)+1);
s1_m1=s_td(mod(m1:30+m1,31)+1);

c0=c_td(mod(n_id_2:30+n_id_2,31)+1);
c1=c_td(mod(n_id_2+3:30+n_id_2+3,31)+1);

z1_m0=z_td(mod((0:30)+mod(m0,8),31)+1);
z1_m1=z_td(mod((0:30)+mod(m1,8),31)+1);

if (slot_num==0)
  s(2:2:62)=s1_m1.*c1.*z1_m0;
  s(1:2:62)=s0_m0.*c0;
elseif (slot_num==10)
  s(2:2:62)=s0_m0.*c1.*z1_m1;
  s(1:2:62)=s1_m1.*c0;
else
  error('Check code...');
end

