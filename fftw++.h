/* Fast Fourier transform header class for FFTW3 Library
   Version 1.0 Copyright (C) 2004 John C. Bowman

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. */

#ifndef __fftwpp_h__
#define __fftwpp_h__ 1

#define __FFTWPP_H_VERSION__ 1.00

#include <fstream>
#include <iostream>
#include <fftw3.h>

#ifndef __Complex_h__
#include <complex>
typedef std::complex<double> Complex;
#endif

using std::ifstream;
using std::ofstream;
using std::cerr;
using std::endl;

#ifdef __Array_h__
using Array::array1;
using Array::Array1;
using Array::array2;
using Array::Array2;
using Array::array3;
using Array::Array3;

static Array1<Complex> NULL1;  
static Array2<Complex> NULL2;  
static Array3<Complex> NULL3;  
#endif

inline void fftw_export_wisdom(void (*emitter)(char c, ofstream& s),
			       ofstream& s)
{
  fftw_export_wisdom((void (*) (char, void *)) emitter, (void *) &s);
}

inline int fftw_import_wisdom(int (*g)(ifstream& s), ifstream &s)
{
  return fftw_import_wisdom((int (*) (void *)) g, (void *) &s);
}

inline void PutWisdom(char c, ofstream& s) {s.put(c);}
inline int GetWisdom(ifstream& s) {return s.get();}

inline Complex *fftnew(size_t size)
{
  void *mem=fftw_malloc(size*sizeof(Complex));  
  if(size && !mem) cerr << endl << "Memory limits exceeded" << endl;
  return (Complex *) mem;
}

inline void fftdelete(Complex *ptr)
{
  fftw_free(ptr);
}

// Base clase for fft routines
//
class fftw {
protected:
  unsigned int size;
  int sign;
  double norm;

  bool inplace;
  fftw_plan plan;
  
  static unsigned int effort;
  static bool Wise;
  static const char *WisdomName;
  static ifstream ifWisdom;
  static ofstream ofWisdom;
  
  unsigned int Dist(unsigned int n, unsigned int stride, unsigned int dist) {
    return dist ? dist : ((stride == 1) ? n : 1);
  }
  
  unsigned int realsize(unsigned int n, Complex *in, Complex *out) {
    return n/2+(!out || in == out);
  }
  
  // Shift the Fourier origin to (nx/2,0)
  void Shift(Complex *data, unsigned int nx, unsigned int ny) {
    const unsigned int nyp=ny/2+1;
    Complex *pstop=data+nx*nyp;
    int pinc=2*nyp;
    for(Complex *p=data+nyp; p < pstop; p += pinc) {
      //#pragma ivdep
      for(unsigned int j=0; j < nyp; j++) p[j]=-p[j];
    }
  }

  // Shift the Fourier origin to (nx/2,ny/2,0)
  void Shift(Complex *data, unsigned int nx, unsigned int ny,
	     unsigned int nz) {
    const unsigned int nzp=nz/2+1;
    const unsigned int nyzp=ny*nzp;
    const unsigned int pinc=2*nzp;
    Complex *p,*pstop;
    p=pstop=data;
    for(unsigned i=0; i < nx; i++) {
      if(i % 2) p -= nzp;
      else p += nzp;
      pstop += nyzp;
      for(; p < pstop; p += pinc) {
	//#pragma ivdep
	for(unsigned int k=0; k < nzp; k++) p[k]=-p[k];
      }
    }
  }

public:
  fftw(unsigned int size, int sign, unsigned int n=0) : 
    size(size), sign(sign), norm(1.0/(n ? n : size)) {}
  
  virtual ~fftw() {}
  
  virtual fftw_plan Plan(Complex *in, Complex *out)=0;
  
  inline void CheckAlign(Complex *p, const char *s) {
    if((size_t) p % sizeof(Complex) == 0) return;
    cerr << "ERROR: " << s << " array is not " << sizeof(Complex) 
	 << "-byte aligned" << endl;
    exit(1);
  }
  
  void Setup(Complex *in, Complex *out) {
    if(!Wise) LoadWisdom();
    bool alloc=!in;
    if(alloc) in=fftnew(size);
    CheckAlign(in,"constructor input");
    if(out) CheckAlign(out,"constructor output");
    else out=in;
    inplace=(out==in);
    plan=Plan(in,out);
    if(!plan) {
      cerr << "Unable to construct FFTW plan" << endl;
      exit(1);
    }
    
    if(alloc) fftdelete(in);
    SaveWisdom();
  }
  
  void LoadWisdom() {
    ifWisdom.open(WisdomName);
    fftw_import_wisdom(GetWisdom,ifWisdom);
    ifWisdom.close();
    Wise=true;
  }

  void SaveWisdom() {
    ofWisdom.open(WisdomName);
    fftw_export_wisdom(PutWisdom,ofWisdom);
    ofWisdom.close();
  }
  
  virtual void Execute(Complex *in, Complex *out) {
    fftw_execute_dft(plan,(fftw_complex *) in,(fftw_complex *) out);
  }
    
  void Setout(Complex *in, Complex *&out) {
    CheckAlign(in,"input");
    if(out) CheckAlign(out,"output");
    else out=in;
    if(inplace ^ (out == in)) {
      cerr << "ERROR: fft constructor and call must be either both in place or out of place" << endl; 
      exit(1);
    }
  }
  
  void fft(Complex *in, Complex *out=NULL) {
    Setout(in,out);
    Execute(in,out);
  }
    
  virtual void fftNormalized(Complex *in, Complex *out=NULL) {
    Setout(in,out);
    Execute(in,out);
    for(unsigned int i=0; i < size; i++) out[i] *= norm;
  }
  
  void fftNormalized(Complex *in, Complex *out,
		     unsigned int nx, unsigned int m,
		     unsigned int stride, unsigned int dist) {
    if(stride == 1 && dist == nx) fftw::fftNormalized(in,out);
    else if(stride == nx && dist == 1) fftw::fftNormalized(in,out);
    else {
      Setout(in,out);
      Execute(in,out);
      for(unsigned int k=0; k < m; k++) {
	for(unsigned int j=0; j < nx; j++) {
	  out[j*stride+k*dist] *= norm;
	}
      }
    }
  }

};

// Compute the complex Fourier transform of n complex values.
// Before calling fft(), the arrays in and out (which may coincide) must be
// allocated as Complex[n].
//
// Out-of-place usage: 
//
//   fft1d Forward(n,-1,in,out);
//   Forward.fft(in,out);
//
//   fft1d Backward(n,1,out,in);
//   Backward.fft(out,in);
//
//   fft1d Backward(n,1,out,in);
//   Backward.fftNormalized(out,in); // True inverse of Forward.fft(in,out);
//
// In-place usage:
//
//   fft1d Forward(n,-1);
//   Forward.fft(in);
//
//   fft1d Backward(n,1);
//   Backward.fft(in);
//
class fft1d : public fftw {
  unsigned int nx;
public:  
  fft1d(unsigned int nx, int sign, Complex *in=NULL, Complex *out=NULL) 
    : fftw(nx,sign), nx(nx) {Setup(in,out);} 
  
#ifdef __Array_h__
  fft1d(int sign, const array1<Complex>& in, const array1<Complex>& out=NULL1) 
    : fftw(in.Nx(),sign), nx(in.Nx()) {Setup(in,out);} 
  
  fft1d(int sign, const Array1<Complex>& in, const Array1<Complex>& out=NULL1) 
    : fftw(in.Nx(),sign), nx(in.Nx()) {Setup(in,out);} 
#endif  
  
  fftw_plan Plan(Complex *in, Complex *out) {
    return fftw_plan_dft_1d(nx,(fftw_complex *) in, (fftw_complex *) out,
		       sign,effort);
  }
};
  
// Compute the complex Fourier transform of m complex vectors, each of
// length n.
// Before calling fft(), the arrays in and out (which may coincide) must be
// allocated as Complex[m*n].
//
// Out-of-place usage: 
//
//   mfft1d Forward(n,-1,m,stride,dist,in,out);
//   Forward.fft(in,out);
//
// In-place usage:
//
//   mfft1d Forward(n,-1,m,stride,dist);
//   Forward.fft(in);
//
// Notes:
//   stride is the spacing between the elements of each Complex vector;
//   dist is the spacing between the first elements of the vectors;
//
//
class mfft1d : public fftw {
  unsigned int nx;
  unsigned int m;
  unsigned int stride;
  unsigned int dist;
public:  
  mfft1d(unsigned int nx, int sign, unsigned int m=1, unsigned int stride=1,
	 unsigned int dist=0, Complex *in=NULL, Complex *out=NULL) 
    : fftw((nx-1)*stride+(m-1)*Dist(nx,stride,dist)+1,sign,nx),
      nx(nx), m(m), stride(stride), dist(Dist(nx,stride,dist))
  {Setup(in,out);} 
  
  fftw_plan Plan(Complex *in, Complex *out) {
    int n[1]={nx};
    return fftw_plan_many_dft(1,n,m,
			    (fftw_complex *) in,NULL,stride,dist,
			    (fftw_complex *) out,NULL,stride,dist,
			    sign,effort);
  }
  
  void fftNormalized(Complex *in, Complex *out=NULL) {
    fftw::fftNormalized(in,out,nx,m,stride,dist);
  }
};
  
// Compute the complex Fourier transform of n real values, using phase sign -1.
// Before calling fft(), the array in must be allocated as Complex[n/2] and
// the array out must be allocated as Complex[n/2+1]. The arrays in and out
// may coincide, in which case they must both be allocated as Complex[n/2+1].
//
// Out-of-place usage: 
//
//   rcfft1d Forward(n,in,out);
//   Forward.fft(in,out);
//
// In-place usage:
//
//   rcfft1d Forward(n);
//   Forward.fft(out);
// 
// Notes:
//   in contains the n real values stored as a Complex array;
//   out contains the first n/2+1 Complex Fourier values.
//
class rcfft1d : public fftw {
  unsigned int nx;
public:  
  rcfft1d(unsigned int nx, Complex *in=NULL, Complex *out=NULL) 
    : fftw(nx/2+1,-1,nx), nx(nx) {Setup(in,out);} 
  
#ifdef __Array_h__
  rcfft1d(const array1<Complex>& in, const array1<Complex>& out=NULL1)
    : fftw(in.Nx(),-1,2*(in.Nx()-1)), nx(2*(in.Nx()-1)) {Setup(in,out);} 
  
  rcfft1d(const Array1<Complex>& in, const Array1<Complex>& out=NULL1)
    : fftw(in.Nx(),-1,2*(in.Nx()-1)), nx(2*(in.Nx()-1)) {Setup(in,out);} 
#endif  
  
  fftw_plan Plan(Complex *in, Complex *out) {
    return fftw_plan_dft_r2c_1d(nx,(double *) in, (fftw_complex *) out, effort);
  }
  
  void Execute(Complex *in, Complex *out) {
    fftw_execute_dft_r2c(plan,(double *) in,(fftw_complex *) out);
  }
};
  
// Compute the real inverse Fourier transform of the n/2+1 Complex values
// corresponding to the non-negative part of the frequency spectrum, using
// phase sign +1.
// Before calling fft(), the array in must be allocated as Complex[n/2+1]
// and the array out must be allocated as Complex[n/2]. The arrays in and out
// may coincide, in which case they must both be allocated as Complex[n/2+1]. 
//
// Out-of-place usage: 
//
//   crfft1d Backward(n,in,out);
//   Backward.fft(in,out);
//
// In-place usage:
//
//   crfft1d Backward(n);
//   Backward.fft(in);
// 
// Notes:
//   in contains the first n/2+1 Complex Fourier values.
//   out contains the n real values stored as a Complex array;
//
class crfft1d : public fftw {
  unsigned int nx;
public:  
  crfft1d(unsigned int nx, Complex *in=NULL, Complex *out=NULL) 
    : fftw(realsize(nx,in,out),1,nx), nx(nx) {Setup(in,out);} 
  
#ifdef __Array_h__
  crfft1d(const array1<Complex>& in, const array1<Complex>& out=NULL1) 
    : fftw(in.Nx(),1,2*(in.Nx()-1)), nx(2*(in.Nx()-1)) {Setup(in,out);} 
  
  crfft1d(const Array1<Complex>& in, const Array1<Complex>& out=NULL1) 
    : fftw(in.Nx(),1,2*(in.Nx()-1)), nx(2*(in.Nx()-1)) {Setup(in,out);} 
#endif  
  
  fftw_plan Plan(Complex *in, Complex *out) {
    return fftw_plan_dft_c2r_1d(nx,(fftw_complex *) in, (double *) out,effort);
  }
  
  void Execute(Complex *in, Complex *out) {
    fftw_execute_dft_c2r(plan,(fftw_complex *) in,(double *) out);
  }
};
  
// Compute the real Fourier transform of m real vectors, each of length n,
// using phase sign -1. Before calling fft(), the array in must be
// allocated as Complex[m*n/2] and the array out must be allocated as
// Complex[m*(n/2+1)]. The arrays in and out may coincide, in which case
// they must both be allocated as Complex[m*(n/2+1)].
//
// Out-of-place usage: 
//
//   mrcfft1d Forward(n,m,stride,dist,in,out);
//   Forward.fft(in,out);
//
// In-place usage:
//
//   mrcfft1d Forward(n,m,stride,dist);
//   Forward.fft(out);
// 
// Notes:
//   stride is the spacing between the elements of each Complex vector;
//   dist is the spacing between the first elements of the vectors;
//   in contains the n real values stored as a Complex array;
//   out contains the first n/2+1 Complex Fourier values.
//
class mrcfft1d : public fftw {
  unsigned int nx;
  unsigned int m;
  unsigned int stride;
  unsigned int dist;
public:  
  mrcfft1d(unsigned int nx, unsigned int m=1, unsigned int stride=1,
	 unsigned int dist=0, Complex *in=NULL, Complex *out=NULL) 
    : fftw(nx/2*stride+(m-1)*Dist(nx,stride,dist)+1,-1,nx), nx(nx), m(m),
      stride(stride), dist(Dist(nx,stride,dist)) {Setup(in,out);} 
  
  fftw_plan Plan(Complex *in, Complex *out) {
    const int n[1]={nx};
    return fftw_plan_many_dft_r2c(1,n,m,
				(double *) in,NULL,stride,2*dist,
				(fftw_complex *) out,NULL,stride,dist,
				effort);
  }
  
  void Execute(Complex *in, Complex *out) {
    fftw_execute_dft_r2c(plan,(double *) in,(fftw_complex *) out);
  }
  
  void fftNormalized(Complex *in, Complex *out=NULL) {
    fftw::fftNormalized(in,out,(nx/2+1),m,stride,dist);
  }
};
  
// Compute the real inverse Fourier transform of m complex vectors, each of
// length n/2+1, corresponding to the non-negative parts of the frequency
// spectra, using phase sign +1. Before calling fft(), the array in must be
// allocated as Complex[m*(n/2+1)] and the array out must be allocated as
// Complex[m*n/2]. The arrays in and out may coincide, in which case they
// must both be allocated as Complex[m*(n/2+1)].  
//
// Out-of-place usage: 
//
//   mcrfft1d Backward(n,m,stride,dist,in,out);
//   Backward.fft(in,out);
//
// In-place usage:
//
//   mcrfft1d Backward(n,m,stride,dist);
//   Backward.fft(out);
// 
// Notes:
//   stride is the spacing between the elements of each Complex vector;
//   dist is the spacing between the first elements of the vectors;
//   in contains the first n/2+1 Complex Fourier values.
//   out contains the n real values stored as a Complex array;
//
class mcrfft1d : public fftw {
  unsigned int nx;
  unsigned int m;
  unsigned int stride;
  unsigned int dist;
public:
  mcrfft1d(unsigned int nx, unsigned int m=1, unsigned int stride=1,
	 unsigned int dist=0, Complex *in=NULL, Complex *out=NULL) 
    : fftw((realsize(nx,in,out)-1)*stride+(m-1)*Dist(nx,stride,dist)+1,1,nx),
      nx(nx), m(m), stride(stride), dist(Dist(nx,stride,dist)) {Setup(in,out);}
  
  fftw_plan Plan(Complex *in, Complex *out) {
    const int n[1]={nx};
    return fftw_plan_many_dft_c2r(1,n,m,
				(fftw_complex *) in,NULL,stride,dist,
				(double *) out,NULL,stride,2*dist,
				effort);
  }
  
  void Execute(Complex *in, Complex *out) {
    fftw_execute_dft_c2r(plan,(fftw_complex *) in,(double *) out);
  }
  
  void fftNormalized(Complex *in, Complex *out=NULL) {
    fftw::fftNormalized(in,out,(nx/2+1),m,stride,dist);
  }
};
  
// Compute the complex two-dimensional Fourier transform of nx times ny
// complex values. Before calling fft(), the arrays in and out (which may
// coincide) must be allocated as Complex[nx*ny].
//
// Out-of-place usage: 
//
//   fft2d Forward(nx,ny,-1,in,out);
//   Forward.fft(in,out);
//
//   fft2d Backward(nx,ny,1,out,in);
//   Backward.fft(out,in);
//
//   fft2d Backward(nx,ny,1,out,in);
//   Backward.fftNormalized(out,in); // True inverse of Forward.fft(in,out);
//
// In-place usage:
//
//   fft2d Forward(nx,ny,-1);
//   Forward.fft(in);
//
//   fft2d Backward(nx,ny,1);
//   Backward.fft(in);
//
// Note:
//   in[ny*i+j] contains the ny Complex values for each i=0,...,nx-1.
//
class fft2d : public fftw {
  unsigned int nx;
  unsigned int ny;
public:  
  fft2d(unsigned int nx, unsigned int ny, int sign, Complex *in=NULL,
	Complex *out=NULL)
    : fftw(nx*ny,sign), nx(nx), ny(ny) {Setup(in,out);} 
  
#ifdef __Array_h__
  fft2d(int sign, const array2<Complex>& in, const array2<Complex>& out=NULL2) 
    : fftw(in.Nx()*in.Ny(),sign), nx(in.Nx()), ny(in.Ny()) {Setup(in,out);} 
  
  fft2d(int sign, const Array2<Complex>& in, const Array2<Complex>& out=NULL2) 
    : fftw(in.Nx()*in.Ny(),sign), nx(in.Nx()), ny(in.Ny()) {Setup(in,out);} 
#endif  
  
  fftw_plan Plan(Complex *in, Complex *out) {
    return fftw_plan_dft_2d(nx,ny,(fftw_complex *) in, (fftw_complex *) out,
			  sign,effort);
  }
};

// Compute the complex two-dimensional Fourier transform of nx times ny real
// values, using phase sign -1.
// Before calling fft(), the array in must be allocated as Complex[nx*ny/2] and
// the array out must be allocated as Complex[nx*(ny/2+1)]. The arrays in
// and out may coincide, in which case they must both be allocated as
// Complex[nx*(ny/2+1)]. 
//
// Out-of-place usage: 
//
//   rcfft2d Forward(nx,ny,in,out);
//   Forward.fft(in,out);       // Origin of Fourier domain at (nx/2,0)
//   Forward.fft(in,false,out); // Origin of Fourier domain at (0,0)
//
//
// In-place usage:
//
//   rcfft2d Forward(nx,ny);
//   Forward.fft(in);           // Origin of Fourier domain at (nx/2,0)
//   Forward.fft(in,false);     // Origin of Fourier domain at (0,0)
// 
// Notes:
//   in contains the nx*ny real values stored as a Complex array;
//   out contains the upper-half portion (ky >= 0) of the Complex transform.
//
class rcfft2d : public fftw {
  unsigned int nx;
  unsigned int ny;
  bool shift;
public:  
  rcfft2d(unsigned int nx, unsigned int ny, Complex *in=NULL,
	  Complex *out=NULL) 
    : fftw(nx*(ny/2+1),-1,nx*ny), nx(nx), ny(ny), shift(true) {Setup(in,out);} 
  
#ifdef __Array_h__
  rcfft2d(const array2<Complex>& in, const array2<Complex>& out=NULL2)
    : fftw(in.Nx()*in.Ny(),-1,in.Nx()*2*(in.Ny()-1)),
      nx(in.Nx()), ny(2*(in.Ny()-1)), shift(true) {Setup(in,out);} 
  
  rcfft2d(const Array2<Complex>& in, const Array2<Complex>& out=NULL2)
    : fftw(in.Nx()*in.Ny(),-1,in.Nx()*2*(in.Ny()-1)),
      nx(in.Nx()), ny(2*(in.Ny()-1)), shift(true) {Setup(in,out);} 
#endif  
  
  fftw_plan Plan(Complex *in, Complex *out) {
    return fftw_plan_dft_r2c_2d(nx,ny,(double *) in, (fftw_complex *) out,
			      effort);
  }
  
  void fft(Complex *in, bool shift0=true, Complex *out=NULL) {
    shift=shift0;
    fftw::fft(in,out);
  }
    
  void Execute(Complex *in, Complex *out) {
    if(shift) Shift(in,nx,ny);
    fftw_execute_dft_r2c(plan,(double *) in,(fftw_complex *) out);
  }
};
  
// Compute the real two-dimensional inverse Fourier transform of the
// nx*(ny/2+1) Complex values corresponding to the spectral values in the
// half-plane ky >= 0, using phase sign +1.
// Before calling fft(), the array in must be allocated as
// Complex[nx*(ny+1)/2] and the array out must be allocated as
// Complex[nx*ny]. The arrays in and out may coincide, in which case they
// must both be allocated as Complex[nx*(ny/2+1)]. 
//
// Out-of-place usage: 
//
//   crfft2d Backward(nx,ny,in,out);
//   Backward.fft(in,out);       // Origin of Fourier domain at (nx/2,0)
//   Backward.fft(in,false,out); // Origin of Fourier domain at (0,0)
//
//
// In-place usage:
//
//   crfft2d Backward(nx,ny);
//   Backward.fft(in);           // Origin of Fourier domain at (nx/2,0)
//   Backward.fft(in,false);     // Origin of Fourier domain at (0,0)
// 
// Notes:
//   in contains the upper-half portion (ky >= 0) of the Complex transform;
//   out contains the nx*ny real values stored as a Complex array.
//
class crfft2d : public fftw {
  unsigned int nx;
  unsigned int ny;
  bool shift;
public:  
  crfft2d(unsigned int nx, unsigned int ny, Complex *in=NULL,
	  Complex *out=NULL) 
    : fftw(nx*(realsize(ny,in,out)),1,nx*ny), nx(nx), ny(ny), shift(true)
  {Setup(in,out);} 
  
#ifdef __Array_h__
  crfft2d(const array2<Complex>& in, const array2<Complex>& out=NULL2)
    : fftw(in.Nx()*in.Ny(),1,in.Nx()*2*(in.Ny()-1)),
      nx(in.Nx()), ny(2*(in.Ny()-1)), shift(true) {Setup(in,out);} 
  
  crfft2d(const Array2<Complex>& in, const Array2<Complex>& out=NULL2)
    : fftw(in.Nx()*in.Ny(),1,in.Nx()*2*(in.Ny()-1)),
      nx(in.Nx()), ny(2*(in.Ny()-1)), shift(true) {Setup(in,out);} 
#endif  
  
  fftw_plan Plan(Complex *in, Complex *out) {
    return fftw_plan_dft_c2r_2d(nx,ny,(fftw_complex *) in, (double *) out,
			      effort);
  }
  
  void fft(Complex *in, bool shift0=true, Complex *out=NULL) {
    shift=shift0;
    fftw::fft(in,out);
  }
    
  void Execute(Complex *in, Complex *out) {
    fftw_execute_dft_c2r(plan,(fftw_complex *) in,(double *) out);
    if(shift) Shift(out,nx,ny);
  }
};

// Compute the complex three-dimensional Fourier transform of 
// nx times ny times nz complex values. Before calling fft(), the arrays in
// and out (which may coincide) must be allocated as Complex[nx*ny*nz].
//
// Out-of-place usage: 
//
//   fft3d Forward(nx,ny,nz,-1,in,out);
//   Forward.fft(in,out);
//
//   fft3d Backward(nx,ny,nz,1,out,in);
//   Backward.fft(out,in);
//
//   fft3d Backward(nx,ny,nz,1,out,in);
//   Backward.fftNormalized(out,in); // True inverse of Forward.fft(in,out);
//
// In-place usage:
//
//   fft3d Forward(nx,ny,nz,-1);
//   Forward.fft(in);
//
//   fft3d Backward(nx,ny,nz,1);
//   Backward.fft(in);
//
// Note:
//   in[nz*(ny*i+j)+k] contains the (i,j,k)th Complex value,
//   indexed by i=0,...,nx-1, j=0,...,ny-1, and k=0,...,nz-1.
//
class fft3d : public fftw {
  unsigned int nx;
  unsigned int ny;
  unsigned int nz;
public:  
  fft3d(unsigned int nx, unsigned int ny, unsigned int nz,
	int sign, Complex *in=NULL, Complex *out=NULL)
    : fftw(nx*ny*nz,sign), nx(nx), ny(ny), nz(nz) {Setup(in,out);} 
  
#ifdef __Array_h__
  fft3d(int sign, const array3<Complex>& in=NULL3, const array3<Complex>& out=NULL3)
    : fftw(in.Nx()*in.Ny()*in.Nz(),sign),
      nx(in.Nx()), ny(in.Ny()), nz(in.Nz()) {Setup(in,out);} 
  
  fft3d(int sign, const Array3<Complex>& in=NULL3, const Array3<Complex>& out=NULL3)
    : fftw(in.Nx()*in.Ny()*in.Nz(),sign),
      nx(in.Nx()), ny(in.Ny()), nz(in.Nz()) {Setup(in,out);} 
#endif  
  
  fftw_plan Plan(Complex *in, Complex *out) {
    return fftw_plan_dft_3d(nx, ny, nz, (fftw_complex *) in,
			  (fftw_complex *) out, sign, effort);
  }
};

// Compute the complex two-dimensional Fourier transform of
// nx times ny times nz real values, using phase sign -1.
// Before calling fft(), the array in must be allocated as Complex[nx*ny*nz/2]
// and the array out must be allocated as Complex[nx*ny*(nz/2+1)]. The
// arrays in and out may coincide, in which case they must both be allocated as
// Complex[nx*ny*(nz/2+1)]. 
//
// Out-of-place usage: 
//
//   rcfft3d Forward(nx,ny,nz,in,out);
//   Forward.fft(in,out);       // Origin of Fourier domain at (nx/2,ny/2,0)
//   Forward.fft(in,false,out); // Origin of Fourier domain at (0,0)
//
//
// In-place usage:
//
//   rcfft3d Forward(nx,ny,nz);
//   Forward.fft(in);           // Origin of Fourier domain at (nx/2,ny/2,0)
//   Forward.fft(in,false);     // Origin of Fourier domain at (0,0)
// 
// Notes:
//   in contains the nx*ny*nz real values stored as a Complex array;
//   out contains the upper-half portion (kz >= 0) of the Complex transform.
//
class rcfft3d : public fftw {
  unsigned int nx;
  unsigned int ny;
  unsigned int nz;
  bool shift;
public:  
  rcfft3d(unsigned int nx, unsigned int ny, unsigned int nz, Complex *in=NULL,
	  Complex *out=NULL) : fftw(nx*ny*(nz/2+1),-1,nx*ny*nz),
			       nx(nx), ny(ny), nz(nz),
			       shift(true) {Setup(in,out);} 
  
#ifdef __Array_h__
  rcfft3d(const array3<Complex>& in=NULL3, const array3<Complex>& out=NULL3) :
    fftw(in.Nx()*in.Ny()*in.Nz(),-1,in.Nx()*in.Ny()*2*(in.Nz()-1)),
    nx(in.Nx()), ny(in.Ny()), nz(2*(in.Nz()-1)), shift(true) {Setup(in,out);} 
  
  rcfft3d(const Array3<Complex>& in=NULL3, const Array3<Complex>& out=NULL3) :
    fftw(in.Nx()*in.Ny()*in.Nz(),-1,in.Nx()*in.Ny()*2*(in.Nz()-1)),
    nx(in.Nx()), ny(in.Ny()), nz(2*(in.Nz()-1)), shift(true) {Setup(in,out);} 
#endif  
  
  fftw_plan Plan(Complex *in, Complex *out) {
    return fftw_plan_dft_r2c_3d(nx,ny,nz,(double *) in, (fftw_complex *) out,
			      effort);
  }
  
  void fft(Complex *in, bool shift0=true, Complex *out=NULL) {
    shift=shift0;
    fftw::fft(in,out);
  }
    
  void Execute(Complex *in, Complex *out) {
    if(shift) Shift(in,nx,ny,nz);
    fftw_execute_dft_r2c(plan,(double *) in,(fftw_complex *) out);
  }
};
  
// Compute the real two-dimensional inverse Fourier transform of the
// nx*ny*(nz/2+1) Complex values corresponding to the spectral values in the
// half-plane kz >= 0, using phase sign +1.
// Before calling fft(), the array in must be allocated as
// Complex[nx*ny*(nz+1)/2] and the array out must be allocated as
// Complex[nx*ny*nz]. The arrays in and out may coincide, in which case they
// must both be allocated as Complex[nx*ny*(nz/2+1)]. 
//
// Out-of-place usage: 
//
//   crfft3d Backward(nx,ny,nz,in,out);
//   Backward.fft(in,out);       // Origin of Fourier domain at (nx/2,ny/2,0)
//   Backward.fft(in,false,out); // Origin of Fourier domain at (0,0)
//
//
// In-place usage:
//
//   crfft3d Backward(nx,ny,nz);
//   Backward.fft(in);           // Origin of Fourier domain at (nx/2,ny/2,0)
//   Backward.fft(in,false);     // Origin of Fourier domain at (0,0)
// 
// Notes:
//   in contains the upper-half portion (kz >= 0) of the Complex transform;
//   out contains the nx*ny*nz real values stored as a Complex array.
//
class crfft3d : public fftw {
  unsigned int nx;
  unsigned int ny;
  unsigned int nz;
  bool shift;
public:  
  crfft3d(unsigned int nx, unsigned int ny, unsigned int nz, Complex *in=NULL,
	  Complex *out=NULL) : fftw(nx*ny*(realsize(nz,in,out)),1,nx*ny*nz),
			       nx(nx), ny(ny), nz(nz),
			       shift(true) {Setup(in,out);} 
  
#ifdef __Array_h__
  crfft3d(const array3<Complex>& in=NULL3, const array3<Complex>& out=NULL3) :
    fftw(in.Nx()*in.Ny()*in.Nz(),1,in.Nx()*in.Ny()*2*(in.Nz()-1)),
    nx(in.Nx()), ny(in.Ny()), nz(2*(in.Nz()-1)), shift(true) {Setup(in,out);} 
  
  crfft3d(const Array3<Complex>& in=NULL3, const Array3<Complex>& out=NULL3) :
    fftw(in.Nx()*in.Ny()*in.Nz(),1,in.Nx()*in.Ny()*2*(in.Nz()-1)),
    nx(in.Nx()), ny(in.Ny()), nz(2*(in.Nz()-1)), shift(true) {Setup(in,out);} 
#endif  
  
  fftw_plan Plan(Complex *in, Complex *out) {
    return fftw_plan_dft_c2r_3d(nx,ny,nz,(fftw_complex *) in, (double *) out,
			      effort);
  }
  
  void fft(Complex *in, bool shift0=true, Complex *out=NULL) {
    shift=shift0;
    fftw::fft(in,out);
  }
    
  void Execute(Complex *in, Complex *out) {
    fftw_execute_dft_c2r(plan,(fftw_complex *) in,(double *) out);
    if(shift) Shift(out,nx,ny,nz);
  }
};
  
#endif
