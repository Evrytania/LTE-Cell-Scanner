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

#ifndef HAVE_ITPP_EXT_H
#define HAVE_ITPP_EXT_H

namespace itpp_ext {

// Implement the matlab range notation a:b:c. Ex, 1:2:6 is [1 3 5].
itpp::vec matlab_range(const double first, const double incr, const double last);
// Implement the matlab range notation a:b. Ex, 1:6 is [1 2 3 4 5 6].
itpp::vec matlab_range(const double first, const double last);
// Implement the matlab range notation a:b:c. Ex, 1:2:6 is [1 3 5].
itpp::ivec matlab_range(const int32 first, const int32 incr, const int32 last);
// Implement the matlab range notation a:b. Ex, 1:6 is [1 2 3 4 5 6].
itpp::ivec matlab_range(const int32 first, const int32 last);
// Implement the matlab range notation a:b:c. Ex, 1:2:6 is [1 3 5].
itpp::ivec matlab_range(const uint32 first, const uint32 incr, const uint32 last);
// Implement the matlab range notation a:b. Ex, 1:6 is [1 2 3 4 5 6].
itpp::ivec matlab_range(const uint32 first, const uint32 last);

// itpp has a mod function but it only operates like the Matlab mod
// function on integers. The matlab_mod function defined here should
// operate exactly like the matalab mod function on both integers and
// floats.
inline double matlab_mod(const double k,const double n) {
  return (n==0)?k:(k-n*itpp::floor_i(k/n));
}
inline double matlab_mod(const double k,const unsigned int n) {
  return (n==0)?k:(k-n*itpp::floor_i(k/n));
}
inline int matlab_mod(const int k,const int n) {
  return (n==0)?k:(k-n*itpp::floor_i((double)k/n));
}
inline unsigned int matlab_mod(const unsigned int k,const unsigned int n) {
  return (n==0)?k:(k-n*itpp::floor_i((double)k/n));
}
inline int matlab_mod(const unsigned int k,const int n) {
  return (n==0)?k:(k-n*itpp::floor_i((double)k/n));
}
inline unsigned int matlab_mod(const int k,const unsigned int n) {
  return (n==0)?k:(k-n*itpp::floor_i((double)k/n));
}

// mod() of a vector
// This should already be in itpp...
template <class T>
itpp::Vec <T> matlab_mod(const itpp::Vec <T> v,const unsigned int n) {
  itpp::Vec <T> r(length(v));
  for (int t=0;t<length(v);t++) {
    r(t)=matlab_mod(v(t),n);
  }
  return r;
}

// Flatten a 3d complex float vector into a cvec.
// Implements matlab's m(:) operation.
itpp::cvec flatten(const vcf3d & m);
itpp::vec flatten(const vf3d & m);

// Similar to Matlab's diff function. Return v(2:end)-v(1:end-1).
template <class T>
itpp::Vec <T> diff(const itpp::Vec <T> & v) {
  itpp::Vec <T> retval;
  if (v.length()>=2) {
    retval.set_size(v.length()-1,false);
    for (int32 t=0;t<v.length()-1;t++) {
      retval(t)=v(t+1)-v(t);
    }
  }
  return retval;
}

// Logical AND of all the bits.
// Is there another way to do this???
bool and_reduce(
  const itpp::bvec & v
);
bool and_reduce(
  const itpp::ivec & v
);

// Return the last element of a vector
template <class T>
inline T last(const itpp::Vec <T> & v) {
  return v(length(v)-1);
}

// Read an rtl_sdr capture file into a cvec.
void rtl_sdr_to_cvec(const std::string & filename,itpp::cvec & v);

}

#endif

