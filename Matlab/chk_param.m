function y=chk_param(val,name,varargin)

% y=chk_param(val,name,varargin)
%
% Check that 'val' satisfies all the requirements given in 'varargin'.
% If any of the requirements are not met, the parameter's 'name' is printed
% along with information on which requirement failed.
%
% Example usage:
%   error(chk_param(speed,'speed','>',0,'<=',100));
% The above will throw an error if speed is <=0 or >100. Otherwise,
% nothing will be printed.
%
% Supported parameter checks are:
%   'string'
%   TBD

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

global CHK_PARAM_BYPASS
if (CHK_PARAM_BYPASS==1)
  y=[];
  return
end

failed=0;
t=1;
while (t<=length(varargin))
  switch lower(varargin{t})
    case 'scalar'
      if (~isscalar(val)) failed=1; break; end
    case 'vector'
      if (~isvector(val)) failed=1; break; end
    case 'horizontal'
      if (size(val,1)~=1) failed=1; break; end
    case 'vertical'
      if (size(val,2)~=1) failed=1; break; end
    case 'string'
      if (~ischar(val)) failed=1; break; end
    case 'numeric'
      if (~isnumeric(val)) failed=1; break; end
    case 'struct'
      if (~isstruct(val)) failed=1; break; end
    case 'real'
      if (~isreal(val)) failed=1; break; end
    case '>'
      t=t+1;
      if (t>length(varargin)) error('> requires an argument'); end
      if (~all(val(:)>varargin{t})) failed=1; break; end
    case '<'
      t=t+1;
      if (t>length(varargin)) error('< requires an argument'); end
      if (~all(val(:)<varargin{t})) failed=1; break; end
    case '<='
      t=t+1;
      if (t>length(varargin)) error('<= requires an argument'); end
      if (~all(val(:)<=varargin{t})) failed=1; break; end
    case '>='
      t=t+1;
      if (t>length(varargin)) error('>= requires an argument'); end
      if (~all(val(:)>=varargin{t})) failed=1; break; end
    case 'rows>'
      t=t+1;
      if (t>length(varargin)) error('rows> requires an argument'); end
      if (size(val,1)<=varargin{t}) failed=1; break; end
    case 'rows<'
      t=t+1;
      if (t>length(varargin)) error('rows< requires an argument'); end
      if (size(val,1)>=varargin{t}) failed=1; break; end
    case 'rows<='
      t=t+1;
      if (t>length(varargin)) error('rows<= requires an argument'); end
      if (size(val,1)>varargin{t}) failed=1; break; end
    case 'rows>='
      t=t+1;
      if (t>length(varargin)) error('rows>= requires an argument'); end
      if (size(val,1)<varargin{t}) failed=1; break; end
    case 'rows=='
      t=t+1;
      if (t>length(varargin)) error('rows>= requires an argument'); end
      if (size(val,1)~=varargin{t}) failed=1; break; end
    case 'cols>'
      t=t+1;
      if (t>length(varargin)) error('cols> requires an argument'); end
      if (size(val,2)<=varargin{t}) failed=1; break; end
    case 'cols<'
      t=t+1;
      if (t>length(varargin)) error('cols< requires an argument'); end
      if (size(val,2)>=varargin{t}) failed=1; break; end
    case 'cols<='
      t=t+1;
      if (t>length(varargin)) error('cols<= requires an argument'); end
      if (size(val,2)>varargin{t}) failed=1; break; end
    case 'cols>='
      t=t+1;
      if (t>length(varargin)) error('cols>= requires an argument'); end
      if (size(val,2)<varargin{t}) failed=1; break; end
    case 'cols=='
      t=t+1;
      if (t>length(varargin)) error('cols== requires an argument'); end
      if (size(val,2)~=varargin{t}) failed=1; break; end
    case 'length>'
      t=t+1;
      if (t>length(varargin)) error('length> requires an argument'); end
      if (length(val)<=varargin{t}) failed=1; break; end
    case 'length<'
      t=t+1;
      if (t>length(varargin)) error('length< requires an argument'); end
      if (length(val)>=varargin{t}) failed=1; break; end
    case 'length<='
      t=t+1;
      if (t>length(varargin)) error('length<= requires an argument'); end
      if (length(val)>varargin{t}) failed=1; break; end
    case 'length>='
      t=t+1;
      if (t>length(varargin)) error('length>= requires an argument'); end
      if (length(val)<varargin{t}) failed=1; break; end
    case 'length=='
      t=t+1;
      if (t>length(varargin)) error('length== requires an argument'); end
      if (length(val)~=varargin{t}) failed=1; break; end
    case 'integer'
      if (any(val(:)~=floor(val(:)))) failed=1; break; end
    otherwise
      if (ischar(varargin{t}))
        error(sprintf('Unrecognized requirement passed to chk_param: %s',varargin{t}));
      else
        error(sprintf('Unrecognized requirement passed to chk_param index %i',t));
      end
  end
  t=t+1;
end

if (failed)
  y=['parameter ''' name ''' must satisfy the following requirements:'];
  for t=1:length(varargin)
    if (t~=1)
      y=[y ', '];
    else
      y=[y ' '];
    end
    if (ischar(varargin{t}))
      y=[y varargin{t}];
    else
      y=[y num2str(varargin{t})];
    end
  end
else
  y=[];
end

