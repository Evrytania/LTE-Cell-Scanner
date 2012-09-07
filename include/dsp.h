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

#ifndef HAVE_DSP_H
#define HAVE_DSP_H

// Return the average power of a vector.
template <class myType>
double sigpower(const myType v) {
  double r=0;
  for (int t=0;t<length(v);t++) {
    r+=pow(real(v(t)),2)+pow(imag(v(t)),2);
  }
  return (r/length(v));
}

// Wrapers to properly scale fft and ifft output so that
// sigpower(fft(x))==sigpower(x).
#define idft(A) (ifft(A)*sqrt(length(A)))
#define dft(A) (fft(A)/sqrt(length(A)))

// Shift a vector in frequency
//itpp::cvec fshift(const itpp::cvec & v,const double f,const double fs);
//itpp::cvec fshift(const itpp::cvec & v,const double f);
// Shift vector seq up by f Hz assuming that seq was sampled at fs Hz.
inline itpp::cvec fshift(const itpp::cvec &seq,const double f,const double fs) {
  //std::complex <double> k=std::complex<double>(0,pi*f/(fs/2));
  double k=itpp::pi*f/(fs/2);
  const uint32 len=length(seq);
  itpp::cvec r(len);
  std::complex<double>coeff;
  for (uint32 t=0;t<len;t++) {
    coeff.real()=cos(k*t);
    coeff.imag()=sin(k*t);
    //r(t)=seq(t)*exp(k*((double)t));
    r(t)=seq(t)*coeff;
  }
  return r;
}
// Shift vector seq up by f Hz assuming that seq was sampled at 2 Hz.
inline itpp::cvec fshift(const itpp::cvec &seq,const double f) {
  return fshift(seq,f,2);
}
inline void fshift_inplace(itpp::cvec &seq,const double f,const double fs) {
  //std::complex <double> k=std::complex<double>(0,pi*f/(fs/2));
  double k=itpp::pi*f/(fs/2);
  const uint32 len=length(seq);
  std::complex<double>coeff;
  for (uint32 t=0;t<len;t++) {
    coeff.real()=cos(k*t);
    coeff.imag()=sin(k*t);
    //r(t)=seq(t)*exp(k*((double)t));
    seq(t)*=coeff;
  }
}
inline void fshift_inplace(itpp::cvec &seq,const double f) {
  fshift_inplace(seq,f,2);
}

// Shift a vector in time
template <class vectype>
// Cyclically shift vector to the right by n samples.
void tshift(vectype &v,const double n) {
  // Currently, only support integer values of n.
  ASSERT(n==itpp::floor_i(n));
  if (n>=0) {
    vectype v_save=v.right(n);
    for (uint32 t=v.length()-1;t>=n;t--) {
      v[t]=v[t-n];
    }
    for (uint32 t=0;t<n;t++) {
      v[t]=v_save[t];
    }
  } else {
    vectype v_save=v.left(n);
    for (uint32 t=0;t<v.length()-n;t++) {
      v[t]=v[t+n];
    }
    for (uint32 t=0;t<n;t++) {
      v[t+v.length()-n]=v_save[t];
    }
  }
}

// Convert to/from dB and linear and power values
template <class myType>
myType db10(const myType s) {
  return (10*log10(s));
}

template <class myType>
myType db20(const myType s) {
  return (20*log10(s));
}

template <class myType>
myType udb10(const myType s) {
  myType result;

  result.set_length(length(s),false);
  for (int t=0;t<length(s);t++) {
    result(t)=10^(s(t)/10);
  }

  return result;
}
inline double udb10(const double v) {
  return pow(10.0,v/10.0);
}

template <class myType>
myType udb20(const myType s) {
  myType result;

  result.set_length(length(s),false);
  for (int t=0;t<length(s);t++) {
    result(t)=10^(s(t)/20);
  }

  return result;
}

inline double udb20(const double v) {
  return pow(10,v/20);
}

// Currently, this simply returns a vector of complex gaussian noise with
// zero mean and average power 1.0.
inline itpp::cvec blnoise(
  const uint32 & n_samp
) {
  return itpp::randn_c(n_samp);
}

// Similar to Matlab's interp1 function except that only linear interpolation
// is supported.
template <class T>
itpp::Vec <T> interp1(
  const itpp::vec & X,
  const itpp::Vec <T> & Y,
  const itpp::vec & x
) {
  ASSERT(itpp_ext::and_reduce(itpp_ext::diff(X)>0));
  ASSERT(length(X)==length(Y));

  itpp::Vec <T> retval;
  retval.set_size(length(x));

  // No interpolation possible if only one X,Y pair is supplied.
  if (length(X)==1) {
    retval=Y(0);
    return retval;
  }

  for (int32 t=0;t<length(x);t++) {
    uint32 try_l=0;
    uint32 try_r=length(X)-1;
    while (try_r-try_l>1) {
      uint32 try_mid=itpp::round_i((try_r+try_l)/2.0);
      //std::cout << try_l << " " << try_mid << " " << try_r << std::endl;
      if (x(t)>=X(try_mid)) {
        try_l=try_mid;
      } else {
        try_r=try_mid;
      }
    }
    retval(t)=Y(try_l)+(x(t)-X(try_l))*(Y(try_r)-Y(try_l))/(X(try_r)-X(try_l));
  }

  return retval;
}

// Returns the value x such that chi2cdf(x,k) is p.
inline double chi2cdf_inv(
  const double & p,
  const double & k
) {
  return 2*boost::math::gamma_p_inv(k/2,p);
}

// Returns the chi squared cdf evaluated at x for k degrees of freedom.
inline double chi2cdf(
  const double & x,
  const double & k
) {
  return boost::math::gamma_p(k/2,x/2);
}

// Interpolate, using fft's, time domain signal x so that is of length n_y.
itpp::cvec interpft(
  const itpp::cvec & x,
  const uint32 & n_y_pre
);

#endif

