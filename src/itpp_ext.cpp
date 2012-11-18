// Copyright 2012 Evrytania LLC (http://www.evrytania.com)
//
// Written by James Peroulas <james@evrytania.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <itpp/itbase.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <curses.h>
#include "rtl-sdr.h"
#include "common.h"
#include "macros.h"

// 'Extensions' to itpp that probably should have been included in itpp.
namespace itpp_ext {

using namespace itpp;
using namespace std;

// Flatten 3D vector m. Similar to Matlab's m(:) functionality.
//
// This should only be used for debugging since it is not very fast.
cvec flatten(const vcf3d & m) {
  uint32 d1_len=m.size();
  uint32 d2_len=m[0].size();
  uint32 d3_len=m[0][0].size();

  // Allocate space without initializing values.
  cvec r;
  r.set_size(d1_len*d2_len*d3_len);

  // Extract the elements. First rows, then cols, then the rest of
  // the dimensions.
  uint32 idx=0;
  for (uint32 d3=0;d3<d3_len;d3++) {
    //ASSERT(m[d3].size()==d3_len);
    for (uint32 d2=0;d2<d2_len;d2++) {
      //ASSERT(m[d2].size()==d2_len);
      for (uint32 d1=0;d1<d1_len;d1++) {
        //if (d3==0) {
        //  ASSERT(d3_len==m[d1][d2].size());
        //}
        r[idx++]=m[d1][d2][d3];
      }
    }
  }

  return r;
}

// This should only be used for debugging since it is not very fast.
// Code has been duplicated from flatten(vcf3d) ... :( :(
vec flatten(const vf3d & m) {
  uint32 d1_len=m.size();
  uint32 d2_len=m[0].size();
  uint32 d3_len=m[0][0].size();
  //cout << d1_len << " " << d2_len << " " << d3_len << endl;

  // Allocate space without initializing values.
  vec r;
  r.set_size(d1_len*d2_len*d3_len);

  // Extract the elements. First rows, then cols, then the rest of
  // the dimensions.
  uint32 idx=0;
  for (uint32 d3=0;d3<d3_len;d3++) {
    //ASSERT(m[d3].size()==d3_len);
    for (uint32 d2=0;d2<d2_len;d2++) {
      //ASSERT(m[d2].size()==d2_len);
      for (uint32 d1=0;d1<d1_len;d1++) {
        //if (d3==0) {
        //  ASSERT(d3_len==m[d1][d2].size());
        //}
        r[idx++]=m[d1][d2][d3];
      }
    }
  }

  return r;
}

// Implement the matlab range notation a:b:c. Ex, 1:2:6 is [1 3 5].
itpp::vec matlab_range(const double first, const double incr, const double last) {
  itpp::vec r;
  //it_assert_debug(incr!=0,"increment is zero...");
  ASSERT(incr!=0,"increment is zero...");
  if (sign(last-first)*sign(incr)>=0) {
    r.set_length(floor_i((last-first)/incr)+1,false);
    for (int t=0;t<length(r);t++) {
      r(t)=first+t*incr;
    }
  }
  return r;
}
// Implement the matlab range notation a:b. Ex, 1:6 is [1 2 3 4 5 6].
itpp::vec matlab_range(const double first, const double last) {
  return matlab_range(first,1,last);
}

// Implement the matlab range notation a:b:c. Ex, 1:2:6 is [1 3 5].
itpp::ivec matlab_range(const int32 first, const int32 incr, const int32 last) {
  itpp::ivec r;
  //it_assert_debug(incr!=0,"increment is zero...");
  ASSERT(incr!=0,"increment is zero...");
  if (sign(last-first)*sign(incr)>=0) {
    r.set_length(floor_i((last-first)/((double)incr))+1,false);
    for (int t=0;t<length(r);t++) {
      r(t)=first+t*incr;
    }
  }
  //std::cout << first << incr << last << std::endl;
  //std::cout << r << std::endl;
  return r;
}
// Implement the matlab range notation a:b. Ex, 1:6 is [1 2 3 4 5 6].
itpp::ivec matlab_range(const int32 first, const int32 last) {
  return matlab_range(first,1,last);
}

// Implement the matlab range notation a:b:c. Ex, 1:2:6 is [1 3 5].
itpp::ivec matlab_range(const uint32 first, const uint32 incr, const uint32 last) {
  itpp::ivec r;
  //it_assert_debug(incr!=0,"increment is zero...");
  ASSERT(incr!=0,"increment is zero...");
  ASSERT(last>=first);
  r.set_length(floor_i((last-first)/((double)incr))+1,false);
  for (int t=0;t<length(r);t++) {
    r(t)=first+t*incr;
  }
  return r;
}
// Implement the matlab range notation a:b. Ex, 1:6 is [1 2 3 4 5 6].
itpp::ivec matlab_range(const uint32 first, const uint32 last) {
  return matlab_range(first,(unsigned)1,last);
}

// Logical AND of all the bits.
// Is there another way to do this???
bool and_reduce(
  const itpp::bvec & v
) {
  bool retval=true;
  int32 t=0;
  while (t<length(v)) {
    if ((retval=retval&&v(t++))==0) break;
  }
  return retval;
}
bool and_reduce(
  const itpp::ivec & v
) {
  bool retval=true;
  int32 t=0;
  while (t<length(v)) {
    if ((retval=retval&&v(t++))==0) break;
  }
  return retval;
}


// Read data captured by the rtl_sdr program into a cvec.
void rtl_sdr_to_cvec(
  const string & filename,
  cvec & v
) {
  // Get filesize
  struct stat filestatus;
  stat(filename.c_str(),&filestatus);
  uint32 file_length=filestatus.st_size;
  //cout << "file length: " << file_length << " bytes\n";
  if (floor(file_length/2.0)!=file_length/2.0) {
    cout << "Warning: file contains an odd number of samples" << endl;
  }
  uint32 n_samp=floor(file_length/2.0);

  // Open file
  FILE *file;
  file=fopen(filename.c_str(),"rb");
  if (!file) {
    cerr << "Error: could not open input file" << endl;
    ABORT(-1);
  }

  // Read entire file, all at once!
  uint8 * buffer=(uint8 *)malloc(n_samp*2*sizeof(uint8));
  uint32 n_read=fread(buffer,1,n_samp*2,file);
  if (n_read!=2*n_samp) {
    cerr << "Error: error while reading file" << endl;
    ABORT(-1);
  }

  // Convert to cvec
  v.set_size(n_samp);
  for (uint32 t=0;t<n_read-1;t+=2) {
    complex <double> sample=complex <double>((buffer[t]-127.0)/128.0,(buffer[t+1]-127.0)/128.0);
    // Append to vector.
    v(t>>1)=sample;
  }

  // Cleanup
  free(buffer);
  fclose(file);
}

}

