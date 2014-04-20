function c_est=lte_conv_decode(d_est);

% c_est=lte_conv_decode(d_est);
%
% Apply the Viterbi algorithm to soft bits d_est to produce the
% c_est hard decisions.
%
% d_est must contain three rows.

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
error(chk_param(d_est,'d_est','real','rows==',3,'>=',0,'<=',1));

n_c=size(d_est,2);

% Make sure d_est is neither +1 nor -1
d_est(d_est<16*eps)=16*eps;
d_est(d_est>1-16*eps)=1-16*eps;

% Initialize the branch structure
rec.c_est=NaN(1,n_c);
rec.metric=0;
for t=63:-1:0
  rec.state=t;
  rec.init_state=rec.state;
  branches(t+1)=rec;
end
% Pre-declare new branches variable
branches_next(128)=branches(1);

% Calculate the state transition table
for t=63:-1:0
  for try_bit=1:-1:0
    [str.new_state str.output]=state_trans(t,try_bit);
    state_table(t+1,try_bit+1)=str;
  end
end
%state_table(1,1)
%state_table(1,2)
%state_table(2,1)
%state_table(2,2)
%state_table(23,1)
%state_table(23,2)

% Calculate all the possible new states
for t=1:n_c
  % Create all possible transitions
  for bn=1:64
    for try_bit=0:1
      rec=branches(bn);

      state_trans_lookup=state_table(rec.state+1,try_bit+1);
      new_state=state_trans_lookup.new_state;
      output=state_trans_lookup.output;

      rec.c_est(t)=try_bit;
      for k=1:3
        rec.metric=rec.metric+log((1-output(k))+(output(k)*2-1)*d_est(k,t));
      end
      rec.state=new_state;
      branches_next((bn-1)*2+try_bit+1)=rec;
    end
  end

  % Sort branches based on their new state
  [dummy ind]=sort([branches_next(:).state]);
  branches_next=branches_next(ind);
  % Verify that nothing 'strange' happened. Can be commented out.
  temp=[branches_next(:).state];
  assert(all(temp(1:2:end)==0:63));
  assert(all(temp(2:2:end)==0:63));

  % Keep only the branch with the highest metric
  for bn=1:2:128
    if (branches_next(bn).metric>branches_next(bn+1).metric)
      branches((bn-1)/2+1)=branches_next(bn);
    else
      branches((bn-1)/2+1)=branches_next(bn+1);
    end
  end
  %[branches(:).state]

end

% Find the path with the highest metric that ended up in the same state
% it started in.
best_index=NaN;
best_metric=-Inf;
best_return_index=NaN;
best_return_metric=-Inf;
for t=1:64
  if (branches(t).metric>best_metric)
    best_index=t;
    best_metric=branches(t).metric;
  end
  if ((branches(t).init_state==branches(t).state)&&(branches(t).metric>best_return_metric))
    best_return_index=t;
    best_return_metric=branches(t).metric;
  end
end

if (isfinite(best_return_index))
  c_est=branches(best_return_index).c_est;
else
  c_est=branches(best_index).c_est;
end

end

function [s_new d]=state_trans(s,inbit);

d=NaN(1,3);
%d(1)=xor(inbit,xor(s(2),xor(s(3),xor(s(5),s(6)))));
%d(2)=xor(inbit,xor(s(1),xor(s(2),xor(s(3),s(6)))));
%d(3)=xor(inbit,xor(s(1),xor(s(2),xor(s(4),s(6)))));
d(1)=mod(inbit+sum(bitget(s,[5 4 2 1])),2);
d(2)=mod(inbit+sum(bitget(s,[6 5 4 1])),2);
d(3)=mod(inbit+sum(bitget(s,[6 5 3 1])),2);

%s_new=[inbit s(1:5)];
s_new=floor(s/2)+inbit*2^5;

end

