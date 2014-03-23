% File:   itload.m
% Brief:  Load an IT++ itfile content to Matlab/Octave workspace
% Author: Tony Ottosson and Adam Piatyszek
%
% Usage: itload("fname.it")
%
% This functions loads all variables from an IT++ file format to the
% Matlab/Octave workspace.
%
% -------------------------------------------------------------------------
%
% Copyright (C) 1995-2010  (see AUTHORS file for a list of contributors)
%
% This file is part of IT++ - a C++ library of mathematical, signal
% processing, speech processing, and communications classes and functions.
%
% IT++ is free software: you can redistribute it and/or modify it under the
% terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any
% later version.
%
% IT++ is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along
% with IT++.  If not, see <http://www.gnu.org/licenses/>.
%
% -------------------------------------------------------------------------

function itload(fname)

[fid, err_msg] = fopen(fname, 'rb', 'ieee-le');
if (fid == -1)
  fname = [fname '.it'];
  [fid, err_msg2] = fopen(fname, 'rb', 'ieee-le');
  if (fid == -1)
    error(err_msg);
  end
end

% Check file size
fseek(fid, 0, 'eof');
file_size = ftell(fid);
fseek(fid, 0, 'bof');

% Read "IT++" magic string
[d, n] = fread(fid, 5, 'char');
if (n ~= 5 | d(1:4) ~= [73 84 43 43]')
  error('Not an IT++ file!');
end

% Check the IT++ file version
if (d(5) ~= 3)
  error('Only IT++ file version 3 is supported by this function!');
end

while (1)
  p = ftell(fid); % Save current file position
  [d1, n] = fread(fid, 3, 'uint64'); % Read header, data, and total block sizes
  name = fgetstr(fid); % Read current variable name
  type = fgetstr(fid); % Read current variable type
  fseek(fid, p+d1(1), 'bof'); % Skip header bytes

  if (length(type) == 0) % A deleted entry -> skip it

  % --- bin ---
  elseif (strcmp(type, 'bin'))
    [d, n] = fread(fid, 1, 'char');
    assignin('caller', name, d);
  % --- int8 (char) ---
  elseif (strcmp(type, 'int8'))
    [d, n] = fread(fid, 1, 'int8');
    assignin('caller', name, d);
  % --- int16 (short) ---
  elseif (strcmp(type, 'int16'))
    [d, n] = fread(fid, 1, 'int16');
    assignin('caller', name, d);
  % --- int32 (int) ---
  elseif (strcmp(type, 'int32'))
    [d, n] = fread(fid, 1, 'int32');
    assignin('caller', name, d);
  % --- float32 (float) ---
  elseif (strcmp(type, 'float32'))
    [d, n] = fread(fid, 1, 'float32');
    assignin('caller', name, d);
  % --- float64 (double) ---
  elseif (strcmp(type, 'float64'))
    [d, n] = fread(fid, 1, 'float64');
    assignin('caller', name, d);
  % --- cfloat32 (complex<float>) ---
  elseif (strcmp(type, 'cfloat32'))
    [d, n] = fread(fid, 2, 'float32');
    d = complex(d(1), d(2));
    assignin('caller', name, d);
  % --- cfloat64 (complex<double>) ---
  elseif (strcmp(type, 'cfloat64'))
    [d, n] = fread(fid, 2, 'float64');
    d = complex(d(1), d(2));
    assignin('caller', name, d);

  % --- bvec ---
  elseif (strcmp(type, 'bvec'))
    [size, n] = fread(fid, 1, 'uint64');
    [d, n] = fread(fid, size, 'char');
    assignin('caller', name, d);
  % --- svec ---
  elseif (strcmp(type, 'svec'))
    [size, n] = fread(fid, 1, 'uint64');
    [d, n] = fread(fid, size, 'int16');
    assignin('caller', name, d);
  % --- ivec ---
  elseif (strcmp(type, 'ivec'))
    [size, n] = fread(fid, 1, 'uint64');
    [d, n] = fread(fid, size, 'int32');
    assignin('caller', name, d);
  % --- fvec ---
  elseif (strcmp(type, 'fvec'))
    [size, n] = fread(fid, 1, 'uint64');
    [d, n] = fread(fid, size, 'float32');
    assignin('caller', name, d);
  % --- dvec ---
  elseif (strcmp(type, 'dvec'))
    [size, n] = fread(fid, 1, 'uint64');
    [d, n] = fread(fid, size, 'float64');
    assignin('caller', name, d);
  % --- fcvec ---
  elseif (strcmp(type, 'fcvec'))
    [size, n] = fread(fid, 1, 'uint64');
    [d, n] = fread(fid, size*2, 'float32');
    d = complex(d(1:2:end), d(2:2:end));
    assignin('caller', name, d);
  % --- dcvec ---
  elseif (strcmp(type, 'dcvec'))
    [size, n] = fread(fid, 1, 'uint64');
    [d, n] = fread(fid, size*2, 'float64');
    d = complex(d(1:2:end), d(2:2:end));
    assignin('caller', name, d);
  % --- string ---
  elseif (strcmp(type, 'string'))
    [size, n] = fread(fid, 1, 'uint64');
    [d, n] = fread(fid, size, 'char');
    d = char(d);
    assignin('caller', name, d);

  % --- bmat ---
  elseif (strcmp(type, 'bmat'))
    [size, n] = fread(fid, 2, 'uint64');
    [d, n] = fread(fid, size', 'char');
    assignin('caller', name, d);
  % --- smat ---
  elseif (strcmp(type, 'smat'))
    [size, n] = fread(fid, 2, 'uint64');
    [d, n] = fread(fid, size', 'int16');
    assignin('caller', name, d);
  % --- imat ---
  elseif (strcmp(type, 'imat'))
    [size, n] = fread(fid, 2, 'uint64');
    [d, n] = fread(fid, size', 'int32');
    assignin('caller', name, d);
  % --- fmat ---
  elseif (strcmp(type, 'fmat'))
    [size, n] = fread(fid, 2, 'uint64');
    [d, n] = fread(fid, size', 'float32');
    assignin('caller', name, d);
  % --- dmat ---
  elseif (strcmp(type, 'dmat'))
    [size, n] = fread(fid, 2, 'uint64');
    [d, n] = fread(fid, size', 'float64');
    assignin('caller', name, d);
  % --- fcmat ---
  elseif (strcmp(type, 'fcmat'))
    [size, n] = fread(fid, 2, 'uint64');
    [d, n] = fread(fid, size(1)*size(2)*2, 'float32');
    d = reshape(complex(d(1:2:end), d(2:2:end)), size(1), size(2));
    assignin('caller', name, d);
  % --- dcmat ---
  elseif (strcmp(type, 'dcmat'))
    [size, n] = fread(fid, 2, 'uint64');
    [d, n] = fread(fid, size(1)*size(2)*2, 'float64');
    d = reshape(complex(d(1:2:end), d(2:2:end)), size(1), size(2));
    assignin('caller', name, d);

  % --- bArray ---
  elseif (strcmp(type, 'bArray'))
    [size, n] = fread(fid, 1, 'uint64');
    [d, n] = fread(fid, size, 'char');
    assignin('caller', name, d);
  % --- sArray ---
  elseif (strcmp(type, 'sArray'))
    [size, n] = fread(fid, 1, 'uint64');
    [d, n] = fread(fid, size, 'int16');
    assignin('caller', name, d);
  % --- iArray ---
  elseif (strcmp(type, 'iArray'))
    [size, n] = fread(fid, 1, 'uint64');
    [d, n] = fread(fid, size, 'int32');
    assignin('caller', name, d);
  % --- fArray ---
  elseif (strcmp(type, 'fArray'))
    [size, n] = fread(fid, 1, 'uint64');
    [d, n] = fread(fid, size, 'float32');
    assignin('caller', name, d);
  % --- dArray ---
  elseif (strcmp(type, 'dArray'))
    [size, n] = fread(fid, 1, 'uint64');
    [d, n] = fread(fid, size, 'float64');
    assignin('caller', name, d);
  % --- fcArray ---
  elseif (strcmp(type, 'fcArray'))
    [size, n] = fread(fid, 1, 'uint64');
    [d, n] = fread(fid, size*2, 'float32');
    d = complex(d(1:2:end), d(2:2:end));
    assignin('caller', name, d);
  % --- dcArray ---
  elseif (strcmp(type, 'dcArray'))
    [size, n] = fread(fid, 1, 'uint64');
    [d, n] = fread(fid, size*2, 'float64');
    d = complex(d(1:2:end), d(2:2:end));
    assignin('caller', name, d);

  % --- bvecArray ---
  elseif (strcmp(type, 'bvecArray'))
    [size, n] = fread(fid, 1, 'uint64');
    clear d2;
    for i=1:size;
      [size2, n] = fread(fid, 1, 'uint64');
      [d, n] = fread(fid, size2, 'char');
      d2{i} = d;
    end
    assignin('caller', name, d2);
  % --- svecArray ---
  elseif (strcmp(type, 'svecArray'))
    [size, n] = fread(fid, 1, 'uint64');
    clear d2;
    for i=1:size;
      [size2, n] = fread(fid, 1, 'uint64');
      [d, n] = fread(fid, size2, 'int16');
      d2{i} = d;
    end
    assignin('caller', name, d2);
  % --- ivecArray ---
  elseif (strcmp(type, 'ivecArray'))
    [size, n] = fread(fid, 1, 'uint64');
    clear d2;
    for i=1:size;
      [size2, n] = fread(fid, 1, 'uint64');
      [d, n] = fread(fid, size2, 'int32');
      d2{i} = d;
    end
    assignin('caller', name, d2);
  % --- vecArray ---
  elseif (strcmp(type, 'vecArray'))
    [size, n] = fread(fid, 1, 'uint64');
    clear d2;
    for i=1:size;
      [size2, n] = fread(fid, 1, 'uint64');
      [d, n] = fread(fid, size2, 'float64');
      d2{i} = d;
    end
    assignin('caller', name, d2);
  % --- cvecArray ---
  elseif (strcmp(type, 'cvecArray'))
    [size, n] = fread(fid, 1, 'uint64');
    clear d2;
    for i=1:size;
      [size2, n] = fread(fid, 1, 'uint64');
      [d, n] = fread(fid, size2*2, 'float64');
      d2{i} = complex(d(1:2:end), d(2:2:end));
    end
    assignin('caller', name, d2);
  % --- stringArray ---
  elseif (strcmp(type, 'stringArray'))
    [size, n] = fread(fid, 1, 'uint64');
    clear d2;
    for i=1:size;
      [size2, n] = fread(fid, 1, 'uint64');
      [d, n] = fread(fid, size2, 'char');
      d2{i} = char(d);
    end
    assignin('caller', name, d2);

  % --- bmatArray ---
  elseif (strcmp(type, 'bmatArray'))
    [size, n] = fread(fid, 1, 'uint64');
    clear d2;
    for i=1:size;
      [size2, n] = fread(fid, 2, 'uint64');
      [d, n] = fread(fid, size2', 'char');
      d2{i} = d;
    end
    assignin('caller', name, d2);
  % --- smatArray ---
  elseif (strcmp(type, 'smatArray'))
    [size, n] = fread(fid, 1, 'uint64');
    clear d2;
    for i=1:size;
      [size2, n] = fread(fid, 2, 'uint64');
      [d, n] = fread(fid, size2', 'int16');
      d2{i} = d;
    end
    assignin('caller', name, d2);
  % --- imatArray ---
  elseif (strcmp(type, 'imatArray'))
    [size, n] = fread(fid, 1, 'uint64');
    clear d2;
    for i=1:size;
      [size2, n] = fread(fid, 2, 'uint64');
      [d, n] = fread(fid, size2', 'int32');
      d2{i} = d;
    end
    assignin('caller', name, d2);
  % --- matArray ---
  elseif (strcmp(type, 'matArray'))
    [size, n] = fread(fid, 1, 'uint64');
    clear d2;
    for i=1:size;
      [size2, n] = fread(fid, 2, 'uint64');
      [d, n] = fread(fid, size2', 'float64');
      d2{i} = d;
    end
    assignin('caller', name, d2);
  % --- cmatArray ---
  elseif (strcmp(type, 'cmatArray'))
    [size, n] = fread(fid, 1, 'uint64');
    clear d2;
    for i=1:size;
      [size2, n] = fread(fid, 2, 'uint64');
      [d, n] = fread(fid, size2(1)*size2(2)*2, 'float64');
      d2{i} = reshape(complex(d(1:2:end), d(2:2:end)), size2(1), size2(2));
    end
    assignin('caller', name, d2);

  % --- else ---
  else
    warning(['Not supported type: ' type]);
  end

  if (p + d1(3) >= file_size)
    break;
  else
    fseek(fid, p+d1(3), 'bof');
  end

end

fclose(fid);



function str = fgetstr(fid)
str = '';
while (1)
  [d, n] = fread(fid, 1, 'char');
  if (d == 0)
    break;
  end
  str = [str char(d)];
end
