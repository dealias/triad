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

#include <fstream>
#include <iostream>
#include <fftw3.h>
#include <complex>


#ifndef COMPLEX
typedef std::complex<double> Complex;
#endif

using std::ifstream;
using std::ofstream;
using std::cerr;
using std::endl;

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
  Complex *mem=(Complex *) fftw_malloc(size*sizeof(Complex));  
  if(size && !mem) cerr << endl << "Memory limits exceeded" << endl;
  return mem;
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
  static const unsigned int effort=FFTW_PATIENT;
  fftw_plan plan;
  
  static bool Wise;
  static const char *WisdomName;
  static ifstream ifWisdom;
  static ofstream ofWisdom;
  
  unsigned int realsize(unsigned int n, Complex *in, Complex *out) {
    return n/2+(!out || in == out);
  }
  
  // Shift the Fourier origin to (nx/2,0)
  void Shift(Complex *data, unsigned int nx, unsigned int nyp) {
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
  fftw(unsigned int size, int sign, unsigned int n=0) : size(size), sign(sign),
				      norm(1.0/(n ? n : size)) {}
  
  virtual void Plan(Complex *in, Complex *out)=0;
  
  void Setup(Complex *in, Complex *out) {
    if(!Wise) LoadWisdom();
    bool alloc=!in;
    if(alloc) in=fftnew(size);
    if(!out) out=in;
    inplace=(out==in);
    
    Plan(in,out);
    
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
    if(!out) out=in;
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
//   fft Forward(n,-1,in,out);
//   Forward.fft(in,out);
//
//   fft Backward(n,1,out,in);
//   Backward.fft(out,in);
//
//   fft Backward(n,1,out,in);
//   Backward.fftNormalized(out,in); // True inverse of Forward.fft(in,out);
//
// In-place usage:
//
//   fft Forward(n,-1);
//   Forward.fft(in);
//
//   fft Backward(n,1);
//   Backward.fft(in);
//
class fft : public fftw {
  unsigned int nx;
public:  
  fft(unsigned int nx, int sign, Complex *in=NULL, Complex *out=NULL) 
    : fftw(nx,sign), nx(nx) {Setup(in,out);} 
  
  void Plan(Complex *in, Complex *out) {
    plan=fftw_plan_dft_1d(nx,(fftw_complex *) in, (fftw_complex *) out,
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
//   mfft Forward(n,-1,m,stride,dist,in,out);
//   Forward.fft(in,out);
//
// In-place usage:
//
//   mfft Forward(n,-1,m,stride,dist);
//   Forward.fft(in);
//
// Notes:
//   stride is the spacing between the elements of each Complex vector;
//   dist is the spacing between the first elements of the vectors;
//
//
class mfft : public fftw {
  unsigned int nx;
  unsigned int m;
  unsigned int stride;
  unsigned int dist;
public:  
  mfft(unsigned int nx, int sign, unsigned int m=1, unsigned int stride=1,
	 unsigned int dist=1, Complex *in=NULL, Complex *out=NULL) 
    : fftw((nx-1)*stride+(m-1)*dist+1,sign,nx), nx(nx), m(m),
      stride(stride), dist(dist) {Setup(in,out);} 
  
  void Plan(Complex *in, Complex *out) {
    int n[1]={nx};
    plan=fftw_plan_many_dft(1,n,m,
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
//   rcfft Forward(n,in,out);
//   Forward.fft(in,out);
//
// In-place usage:
//
//   rcfft Forward(n);
//   Forward.fft(out);
// 
// Notes:
//   in contains the n real values stored as a Complex array;
//   out contains the first n/2+1 Complex Fourier values.
//
class rcfft : public fftw {
  unsigned int nx;
public:  
  rcfft(unsigned int nx, Complex *in=NULL, Complex *out=NULL) 
    : fftw(nx/2+1,-1,nx), nx(nx) {Setup(in,out);} 
  
  void Plan(Complex *in, Complex *out) {
    plan=fftw_plan_dft_r2c_1d(nx,(double *) in, (fftw_complex *) out, effort);
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
//   crfft Backward(n,in,out);
//   Backward.fft(in,out);
//
// In-place usage:
//
//   crfft Backward(n);
//   Backward.fft(in);
// 
// Notes:
//   in contains the first n/2+1 Complex Fourier values.
//   out contains the n real values stored as a Complex array;
//
class crfft : public fftw {
  unsigned int nx;
public:  
  crfft(unsigned int nx, Complex *in=NULL, Complex *out=NULL) 
    : fftw(realsize(nx,in,out),1,nx), nx(nx) {Setup(in,out);} 
  
  void Plan(Complex *in, Complex *out) {
    plan=fftw_plan_dft_c2r_1d(nx,(fftw_complex *) in, (double *) out, effort);
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
//   mrcfft Forward(n,m,stride,dist,in,out);
//   Forward.fft(in,out);
//
// In-place usage:
//
//   mrcfft Forward(n,m,stride,dist);
//   Forward.fft(out);
// 
// Notes:
//   stride is the spacing between the elements of each Complex vector;
//   dist is the spacing between the first elements of the vectors;
//   in contains the n real values stored as a Complex array;
//   out contains the first n/2+1 Complex Fourier values.
//
class mrcfft : public fftw {
  unsigned int nx;
  unsigned int m;
  unsigned int stride;
  unsigned int dist;
public:  
  mrcfft(unsigned int nx, unsigned int m=1, unsigned int stride=1,
	 unsigned int dist=1, Complex *in=NULL, Complex *out=NULL) 
    : fftw(nx/2*stride+(m-1)*dist+1,-1,nx), nx(nx), m(m),
      stride(stride), dist(dist) {Setup(in,out);} 
  
  void Plan(Complex *in, Complex *out) {
    const int n[1]={nx};
    plan=fftw_plan_many_dft_r2c(1,n,m,
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
//   mcrfft Backward(n,m,stride,dist,in,out);
//   Backward.fft(in,out);
//
// In-place usage:
//
//   mcrfft Backward(n,m,stride,dist);
//   Backward.fft(out);
// 
// Notes:
//   stride is the spacing between the elements of each Complex vector;
//   dist is the spacing between the first elements of the vectors;
//   in contains the first n/2+1 Complex Fourier values.
//   out contains the n real values stored as a Complex array;
//
class mcrfft : public fftw {
  unsigned int nx;
  unsigned int m;
  unsigned int stride;
  unsigned int dist;
public:
  mcrfft(unsigned int nx, unsigned int m=1, unsigned int stride=1,
	 unsigned int dist=1, Complex *in=NULL, Complex *out=NULL) 
    : fftw((realsize(nx,in,out)-1)*stride+(m-1)*dist+1,1,nx), nx(nx), m(m),
      stride(stride), dist(dist) {Setup(in,out);}
  
  void Plan(Complex *in, Complex *out) {
    const int n[1]={nx};
    plan=fftw_plan_many_dft_c2r(1,n,m,
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
  
  void Plan(Complex *in, Complex *out) {
    plan=fftw_plan_dft_2d(nx,ny,(fftw_complex *) in, (fftw_complex *) out,
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
  
  void Plan(Complex *in, Complex *out) {
    plan=fftw_plan_dft_r2c_2d(nx,ny,(double *) in, (fftw_complex *) out,
			      effort);
  }
  
  void fft(Complex *in, bool shift0=true, Complex *out=NULL) {
    shift=shift0;
    fftw::fft(in,out);
  }
    
  void Execute(Complex *in, Complex *out) {
    if(shift) Shift(in,nx,ny/2+1);
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
  
  void Plan(Complex *in, Complex *out) {
    plan=fftw_plan_dft_c2r_2d(nx,ny,(fftw_complex *) in, (double *) out,
			      effort);
  }
  
  void fft(Complex *in, bool shift0=true, Complex *out=NULL) {
    shift=shift0;
    fftw::fft(in,out);
  }
    
  void Execute(Complex *in, Complex *out) {
    fftw_execute_dft_c2r(plan,(fftw_complex *) in,(double *) out);
    if(shift) Shift(out,nx,ny/2+1);
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
  
  void Plan(Complex *in, Complex *out) {
    plan=fftw_plan_dft_3d(nx, ny, nz, (fftw_complex *) in,
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
  
  void Plan(Complex *in, Complex *out) {
    plan=fftw_plan_dft_r2c_3d(nx,ny,nz,(double *) in, (fftw_complex *) out,
			      effort);
  }
  
  void fft(Complex *in, bool shift0=true, Complex *out=NULL) {
    shift=shift0;
    fftw::fft(in,out);
  }
    
  void Execute(Complex *in, Complex *out) {
    if(shift) Shift(in,nx,ny,nz/2+1);
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
  
  void Plan(Complex *in, Complex *out) {
    plan=fftw_plan_dft_c2r_3d(nx,ny,nz,(fftw_complex *) in, (double *) out,
			      effort);
  }
  
  void fft(Complex *in, bool shift0=true, Complex *out=NULL) {
    shift=shift0;
    fftw::fft(in,out);
  }
    
  void Execute(Complex *in, Complex *out) {
    fftw_execute_dft_c2r(plan,(fftw_complex *) in,(double *) out);
    if(shift) Shift(out,nx,ny,nz/2+1);
  }
};
  
#endif
