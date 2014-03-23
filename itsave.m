% File:   itsave.m
% Brief:  Saves Matlab/Octave workspace variables to an IT++ itfile
% Author: Tony Ottosson and Adam Piatyszek
%
% Usage: itsave("fname.it", var1, [var2], ...)
%
% This function saves a set of Matlab/Octave workspace variables to an IT++
% file format. Currently, only vectors and 2-D matrices can be saved. The
% type of data saved is detected automatically and can be one of the
% following types: bvec, bmat, ivec, imat, vec, mat, cvec, cmat.
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


function itsave(fname, varargin)

if (nargin > 1)
  vars = varargin;
else
  error('Usage: itsave("fname.it", var1, [var2], ...)');
end

nargs = size(vars,2);

% The current file-version of it_file
file_version = 3;

[fid, err_msg] = fopen(fname, 'wb', 'ieee-le');
if (fid == -1)
  error(err_msg);
end

% Write a file header consisting of "IT++" and a char containing the file
% version number
fprintf(fid, 'IT++%c', file_version);

for ai=1:nargs
  if (exist('OCTAVE_VERSION')) % check for octave
      vname = deblank(argn(ai+1,:)); % octave way of getting parameter name
      is_octave=1; %used by function itsizeof to identify octave
  else
      vname = inputname(ai+1); % matlab way of getting parameter name
      is_octave=0; %used by function itsizeof to identify matlab
  end
  v = vars{ai};

  is_scalar = all(size(v) == 1); % true if scalar (for future use)
  is_vector = (sum(size(v) > 1) <= 1); % true if a vector (or a scalar)
  is_intbin = min(min(v == round(v))); % true if integer or binary

  if ( isreal(v) && is_intbin ) % binary or integer type
    if (max(max(v)) == 1 && min(min(v)) == 0) % binary type
      % Calculate sizes
      if (is_vector)
 	hdr_bytes = 3 * itsizeof(uint64(0),is_octave) + size(vname,2)+1 ...
	    + itsizeof('bvec',is_octave)+1 +1;
	data_bytes = itsizeof(uint64(0),is_octave) + itsizeof(char(0),is_octave) * prod(size(v));
      else % a matrix
	hdr_bytes = 3 * itsizeof(uint64(0),is_octave) + size(vname,2)+1 ...
	    + itsizeof('bmat',is_octave)+1 +1;
	data_bytes = 2 * itsizeof(uint64(0),is_octave) + itsizeof(char(0),is_octave) * prod(size(v));
      end
      block_bytes = hdr_bytes + data_bytes;

      % Write header sizes
      fwrite(fid, [hdr_bytes data_bytes block_bytes], 'uint64');
      % Write variable name as string
      fprintf(fid, '%s%c', vname); fwrite(fid, 0, 'char');

      % Write data type string, empty description string and data size
      if (is_vector)
	fprintf(fid, 'bvec'); fwrite(fid, 0, 'char'); fwrite(fid, 0, 'char');
	fwrite(fid, size(v(:),1), 'uint64');
      else % a matrix
	fprintf(fid, 'bmat'); fwrite(fid, 0, 'char'); fwrite(fid, 0, 'char');
	fwrite(fid, size(v), 'uint64');
      end
      % Write data
      fwrite(fid, v, 'char');

    else % integer type
      % Calculate sizes
      if (is_vector)
	hdr_bytes = 3 * itsizeof(uint64(0),is_octave) + size(vname,2)+1 ...
	    + itsizeof('ivec',is_octave)+1 +1;
	data_bytes = itsizeof(uint64(0),is_octave) + itsizeof(int32(0),is_octave) * prod(size(v));
      else % a matrix
	hdr_bytes = 3 * itsizeof(uint64(0),is_octave) + size(vname,2)+1 ...
	    + itsizeof('imat',is_octave)+1 +1;
	data_bytes = 2 * itsizeof(uint64(0),is_octave) + itsizeof(int32(0),is_octave) * prod(size(v));
      end
      block_bytes = hdr_bytes + data_bytes;

      % Write header sizes
      fwrite(fid, [hdr_bytes data_bytes block_bytes], 'uint64');
      % Write variable name as string
      fprintf(fid, '%s%c', vname); fwrite(fid, 0, 'char');

      % Write data type string, empty description string and data size
      if (is_vector)
	fprintf(fid, 'ivec'); fwrite(fid, 0, 'char'); fwrite(fid, 0, 'char');
	fwrite(fid, size(v(:),1), 'uint64');
      else % a matrix
	fprintf(fid, 'imat'); fwrite(fid, 0, 'char'); fwrite(fid, 0, 'char');
	fwrite(fid, size(v), 'uint64');
      end
      % Write data
      fwrite(fid, v, 'int32');
    end % binary or integer

  elseif (isa(v, 'double')) % double precision floating point type
    if (isreal(v)) % Check if real values
      % Calculate sizes
      if (is_vector)
	hdr_bytes = 3 * itsizeof(uint64(0),is_octave) + size(vname,2)+1 ...
	    + itsizeof('dvec',is_octave)+1 + 1;
	data_bytes = itsizeof(uint64(0),is_octave) + itsizeof(double(0),is_octave) * prod(size(v));
      else % a matrix
	hdr_bytes = 3 * itsizeof(uint64(0),is_octave) + size(vname,2)+1 ...
	    + itsizeof('dmat',is_octave)+1 + 1;
	data_bytes = 2 * itsizeof(uint64(0),is_octave) + itsizeof(double(0),is_octave) ...
	    * prod(size(v));
      end
      block_bytes = hdr_bytes + data_bytes;

      % Write a header
      fwrite(fid, [hdr_bytes data_bytes block_bytes], 'uint64');
      % Writes variable name as string
      fprintf(fid, '%s%c', vname); fwrite(fid, 0, 'char');

      if (is_vector)
	fprintf(fid, 'dvec'); fwrite(fid, 0, 'char'); fwrite(fid, 0, 'char');
	fwrite(fid, size(v(:),1), 'uint64');
      else % a matrix
	fprintf(fid, 'dmat'); fwrite(fid, 0, 'char'); fwrite(fid, 0, 'char');
	fwrite(fid, size(v), 'uint64');
      end
      fwrite(fid, v, 'float64');

    else % complex values
      % Calculate sizes
      if (is_vector)
	hdr_bytes = 3 * itsizeof(uint64(0),is_octave) + size(vname,2)+1 ...
	    + itsizeof('dcvec',is_octave)+1 + 1;
	data_bytes = itsizeof(uint64(0),is_octave) + 2 * itsizeof(double(0),is_octave) ...
	    * prod(size(v));
      else % a matrix
	hdr_bytes = 3 * itsizeof(uint64(0),is_octave) + size(vname,2)+1 ...
	    + itsizeof('dcmat',is_octave)+1 + 1;
	data_bytes = 2 * itsizeof(uint64(0),is_octave) + 2 * itsizeof(double(0),is_octave) ...
	    * prod(size(v));
      end
      block_bytes = hdr_bytes + data_bytes;

      % Writes header sizes
      fwrite(fid, [hdr_bytes data_bytes block_bytes], 'uint64');
      % Write variable name as string
      fprintf(fid, '%s', vname);  fwrite(fid, 0, 'char');

      if (is_vector)
	% Write data type string, empty description string and data size
	fprintf(fid, 'dcvec'); fwrite(fid, 0, 'char'); fwrite(fid, 0, 'char');
	fwrite(fid, size(v(:),1), 'uint64');
	% Write data
	for i=1:size(v(:),1)
	  fwrite(fid, real(v(i)), 'float64');
	  fwrite(fid, imag(v(i)), 'float64');
	end
      else % a matrix
        % Write data type string, empty description string and data size
	fprintf(fid, 'dcmat'); fwrite(fid, 0, 'char'); fwrite(fid, 0, 'char');
	fwrite(fid, size(v), 'uint64');
	% Write data
	for j=1:size(v,2)
	  for i=1:size(v,1)
	    fwrite(fid, real(v(i,j)), 'float64');
	    fwrite(fid, imag(v(i,j)), 'float64');
	  end
	end
      end
    end % real or complex
  else
      warning(['Variable ''' vname ''' is neither a vector nor matrix. Not saved.']);
  end

end

fclose(fid);


% sizeof function (like C-style sizeof)
% returns no. bytes used by a variable
function nbytes=itsizeof(in,is_octave)
  if (~is_octave)
    tmp=whos('in');
    nbytes=tmp.bytes;

    % matlab uses 2 bytes by default for char
    % overwrite using 1 byte for file format
    if (ischar(in))
      nbytes=size(in,2);
    end;
  else
    nbytes=sizeof(in);
  end;
