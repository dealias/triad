/* Derived Fast Fourier transform routines
   Version 1.2 Copyright (C) 1997 John C. Bowman (bowman@math.ualberta.ca)

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

#include "options.h"
#include "fft.h"

Complex *wpTable;
unsigned int wpTableSize=0;

static Complex *WTable;
static unsigned int WTableSize=0;

#ifdef _CRAY
const int offset=1;
#else
const int offset=0;
#endif	

inline Complex expim1(Real phase)
{
  double sinyby2=sin(0.5*phase);
  return Complex(-2.0*sinyby2*sinyby2,sin(phase));
}

void fft_init(unsigned int log2n)
{
  unsigned int n=1 << log2n;
	
  wpTable=new(wpTable,log2n) Complex;
  unsigned int mmax=1 << wpTableSize;
  while (n > mmax) {
    mmax <<= 1;
    wpTable[wpTableSize]=expim1(twopi/mmax);
    wpTableSize++;
  }
}
	
void rfft_init(unsigned int log2n)
{
  if(log2n+1 > wpTableSize) fft_init(log2n+1);
  if(log2n == 0) return;
  unsigned int n4=1 << (log2n-1);
	
  WTable=new(WTable,n4) Complex;
  WTableSize=n4;
  WTable[0]=Complex(0.0,0.5);
	
  Complex wp=1.0+wpTable[log2n];
  for(unsigned int i=1; i < n4; i++) WTable[i]=WTable[i-1]*wp;
}

// Return the Fourier transform of n real values.
// Before calling, data must be allocated as Complex[n/2+1].
// On entry: data contains the n real values stored as a Complex array of
// length n/2;
//           log2n contains the base-2 logarithm of n;
//           isign is the sign (+/- 1) of the phase;
//           [scale is a constant by which the results will be multiplied;]
//           [bitreverse is 0 for a true fft of data (default);
//                         +1 for a fft of bit-reversed data (faster);
//                         -1 for a bit-reversed fft of data (faster).]
// On exit:  data contains the n/2+1 Complex Fourier values.

void rcfft(Complex *data, unsigned int log2n, int isign, Real scale,
	   int bitreverse) 
{
  log2n--;
  unsigned int n2=1 << log2n;
  unsigned int n4=n2 >> 1;
  if(WTableSize != n4) rfft_init(log2n);
	
  fft(data,log2n,1,scale,bitreverse);
	
  data[n2]=data[0].re-data[0].im;
  data[0]=data[0].re+data[0].im;
	
  if(isign == 1) {
    //#pragma ivdep	
    for(unsigned int i=1; i < n4; i++) {
      Real ure=data[i].re, uim=data[i].im;
      Real vre=data[n2-i].re, vim=-data[n2-i].im; 
      Real Are=0.5*(ure+vre), Aim=0.5*(uim+vim);
      ure -= vre; uim -= vim;
      vre=WTable[i].re*ure-WTable[i].im*uim;
      vim=WTable[i].re*uim+WTable[i].im*ure;
      data[i].re=Are-vre;
      data[i].im=Aim-vim;
      data[n2-i].re=Are+vre;
      data[n2-i].im=-Aim-vim;
    }
  } else {
    //#pragma ivdep	
    for(unsigned int i=1; i < n4; i++) {
      Real ure=data[i].re, uim=data[i].im;
      Real vre=data[n2-i].re, vim=-data[n2-i].im; 
      Real Are=0.5*(ure+vre), Aim=0.5*(uim+vim);
      ure -= vre; uim -= vim;
      vre=WTable[i].re*ure-WTable[i].im*uim;
      vim=WTable[i].re*uim+WTable[i].im*ure;
      data[i].re=Are-vre;
      data[i].im=vim-Aim;
      data[n2-i].re=Are+vre;
      data[n2-i].im=Aim+vim;
    }
    data[n4].im=-data[n4].im;
  }
}

// Return the real inverse Fourier transform of the n/2+1
// Complex values corresponding to the non-negative part of the frequency
// spectrum. Before calling, data must be allocated as Complex[n/2+1].
// On entry: data contains the n/2+1 Complex Fourier transform values;
//           log2n contains the base-2 logarithm of n;
//           isign is the sign (+/- 1) of the phase;
//           [scale is a constant by which the results will be multiplied;]
//           [bitreverse is 0 for a true fft of data (default);
//                         +1 for a fft of bit-reversed data (faster);
//                         -1 for a bit-reversed fft of data (faster).]
// On exit:  data contains the n real inverse Fourier transform values
// stored as a Complex array of length n/2.
// Note: To compute a true inverse transform, set scale=1.0/n.

void crfft(Complex *data, unsigned int log2n, int isign, Real scale,
	   int bitreverse)
{		 
  log2n--;
  unsigned int n2=1 << log2n;
  unsigned int n4=n2 >> 1;
  if(WTableSize != n4) rfft_init(log2n);
	
  data[0].im=data[0].re-data[n2].re;
  data[0].re += data[n2].re;
	
  if(isign == 1) {
    //#pragma ivdep
    for(unsigned int i=1; i < n4; i++) {
      Real ure=data[i].re, uim=-data[i].im;
      Real vre=data[n2-i].re, vim=data[n2-i].im; 
      Real Are=ure+vre, Aim=uim+vim;
      ure -= vre; uim -= vim;
      vre=2.0*(WTable[i].re*ure+WTable[i].im*uim);
      vim=2.0*(WTable[i].re*uim-WTable[i].im*ure);
      data[i].re=Are-vre;
      data[i].im=Aim-vim;
      data[n2-i].re=Are+vre;
      data[n2-i].im=-Aim-vim;
    }
    data[n4].re *= 2.0;
    data[n4].im *= -2.0;
  } else {
    //#pragma ivdep
    for(unsigned int i=1; i < n4; i++) {
      Real ure=data[i].re, uim=data[i].im;
      Real vre=data[n2-i].re, vim=-data[n2-i].im; 
      Real Are=ure+vre, Aim=uim+vim;
      ure -= vre; uim -= vim;
      vre=2.0*(WTable[i].re*ure+WTable[i].im*uim);
      vim=2.0*(WTable[i].re*uim-WTable[i].im*ure);
      data[i].re=Are-vre;
      data[i].im=Aim-vim;
      data[n2-i].re=Are+vre;
      data[n2-i].im=-Aim-vim;
    }
    data[n4].re *= 2.0;
    data[n4].im *= 2.0;
  }
	
  fft(data,log2n,-1,scale,bitreverse);
}

// Return the Fourier transform of nk real vectors, each of length n.
// Before calling, data must be allocated as Complex[nk*(n/2+1)].
// On entry: data contains the n real values stored as a Complex array of
// length n/2, for each k=0,...,nk-1;
//           log2n contains the base-2 logarithm of n;
//           isign is the sign (+/- 1) of the phase;
//           nk is the number of Complex vectors;
//           [inc1 is the stride between the elements of each Complex vector;]
//           [inc2 is the stride between first elements of the vectors;]
//           [scale is a constant by which the results will be multiplied;]
//           [bitreverse is 0 for a true fft of data (default);
//                         +1 for a fft of bit-reversed data (faster);
//                         -1 for a bit-reversed fft of data (faster).]
// On exit:  data contains the n/2+1 Complex Fourier values.

void mrcfft0(Complex *data, unsigned int log2n, int isign, unsigned int nk,
	     unsigned int inc1, unsigned int inc2, Real scale, int bitreverse)
{		 
  if(log2n == 0) return;
  log2n--;
  if(inc1 == 0) inc1=nk;
  unsigned int kstop=nk*inc2;
	
  unsigned int n2=1 << log2n;
  unsigned int n4=n2 >> 1;
  if(WTableSize != n4) rfft_init(log2n);
	
  mfft(data,log2n,1,nk,inc1,inc2,scale,bitreverse);
	
  Complex *q=data+n2*inc1;
  //#pragma ivdep	
  for(unsigned int k=0; k < kstop; k += inc2) {
    q[k]=data[k].re-data[k].im;
    data[k]=data[k].re+data[k].im;
  }
	
  if(isign == 1) {
    for(unsigned int i=1; i < n4; i++) {
      Complex *p=data+i*inc1, *q=data+(n2-i)*inc1;
      Real Wre=WTable[i].re, Wim=WTable[i].im;
      //#pragma ivdep	
      for(unsigned int k=0; k < kstop; k += inc2) {
	Real ure=p[k].re, uim=p[k].im;
	Real vre=q[k].re, vim=-q[k].im;
	Real Are=0.5*(ure+vre), Aim=0.5*(uim+vim);
	ure -= vre; uim -= vim;
	vre=Wre*ure-Wim*uim;
	vim=Wre*uim+Wim*ure;
	p[k].re=Are-vre;
	p[k].im=Aim-vim;
	q[k].re=Are+vre;
	q[k].im=-Aim-vim;
      }
    }
  } else {
    for(unsigned int i=1; i < n4; i++) {
      Complex *p=data+i*inc1, *q=data+(n2-i)*inc1;
      Real Wre=WTable[i].re, Wim=WTable[i].im;
      //#pragma ivdep	
      for(unsigned int k=0; k < kstop; k += inc2) {
	Real ure=p[k].re, uim=p[k].im;
	Real vre=q[k].re, vim=-q[k].im;
	Real Are=0.5*(ure+vre), Aim=0.5*(uim+vim);
	ure -= vre; uim -= vim;
	vre=Wre*ure-Wim*uim;
	vim=Wre*uim+Wim*ure;
	p[k].re=Are-vre;
	p[k].im=vim-Aim;
	q[k].re=Are+vre;
	q[k].im=Aim+vim;
      }
    }
    Complex *p=data+n4*inc1;
    //#pragma ivdep	
    for(unsigned int k=0; k < kstop; k += inc2) p[k].im=-p[k].im;
  }
}

// Return the real inverse Fourier transform of nk Complex vectors, of
// length n/2+1, corresponding to the non-negative part of the frequency
// spectrum. Before calling, data must be allocated as Complex[nk*(n/2+1)].
// On entry: data contains the n/2+1 Complex Fourier transform values for
// each k=0,...,nk-1, 
//           log2n contains the base-2 logarithm of n;
//           isign is the sign (+/- 1) of the phase;
//           nk is the number of Complex vectors;
//           [inc1 is the stride between the elements of each Complex vector;]
//           [inc2 is the stride between first elements of the vectors;]
//           [scale is a constant by which the results will be multiplied;]
//           [bitreverse is 0 for a true fft of data (default);
//                         +1 for a fft of bit-reversed data (faster);
//                         -1 for a bit-reversed fft of data (faster).]
// On exit:  data contains the n real inverse Fourier transform values
// stored as a Complex array of length n/2.
// Note: To compute a true inverse transform, set scale=1.0/n.

void mcrfft0(Complex *data, unsigned int log2n, int isign, unsigned int nk,
	     unsigned int inc1, unsigned int inc2, Real scale, int bitreverse)
{		 
  if(log2n == 0) return;
  log2n--;
  if(inc1 == 0) inc1=nk;
  unsigned int kstop=nk*inc2;
	
  unsigned int n2=1 << log2n;
  unsigned int n4=n2 >> 1;
  if(WTableSize != n4) rfft_init(log2n);
	
  Complex *q=data+n2*inc1;
  //#pragma ivdep	
  for(unsigned int k=0; k < kstop; k += inc2) {
    data[k].im=data[k].re-q[k].re;
    data[k].re += q[k].re;
  }
	
  if(isign == 1) {
    for(unsigned int i=1; i < n4; i++) {
      Complex *p=data+i*inc1, *q=data+(n2-i)*inc1;
      Real Wre=2.0*WTable[i].re, Wim=2.0*WTable[i].im;
      //#pragma ivdep
      for(unsigned int k=0; k < kstop; k += inc2) {
	Real ure=p[k].re, uim=-p[k].im;
	Real vre=q[k].re, vim=q[k].im;
	Real Are=ure+vre, Aim=uim+vim;
	ure -= vre; uim -= vim;
	vre=Wre*ure+Wim*uim;
	vim=Wre*uim-Wim*ure;
	p[k].re=Are-vre;
	p[k].im=Aim-vim;
	q[k].re=Are+vre;
	q[k].im=-Aim-vim;
      }
    }
    Complex *p=data+n4*inc1;
    //#pragma ivdep	
    for(unsigned int k=0; k < kstop; k += inc2) p[k]=2.0*conj(p[k]);
  } else {
    for(unsigned int i=1; i < n4; i++) {
      Complex *p=data+i*inc1, *q=data+(n2-i)*inc1;
      Real Wre=2.0*WTable[i].re, Wim=2.0*WTable[i].im;
      //#pragma ivdep
      for(unsigned int k=0; k < kstop; k += inc2) {
	Real ure=p[k].re, uim=p[k].im;
	Real vre=q[k].re, vim=-q[k].im;
	Real Are=ure+vre, Aim=uim+vim;
	ure -= vre; uim -= vim;
	vre=Wre*ure+Wim*uim;
	vim=Wre*uim-Wim*ure;
	p[k].re=Are-vre;
	p[k].im=Aim-vim;
	q[k].re=Are+vre;
	q[k].im=-Aim-vim;
      }
    }
    Complex *p=data+n4*inc1;
    //#pragma ivdep
    for(unsigned int k=0; k < kstop; k += inc2) {
      p[k].re *= 2.0;
      p[k].im *= 2.0;
    }
  }
	
  mfft(data,log2n,-1,nk,inc1,inc2,scale,bitreverse);
}

// Return the two-dimensional Fourier transform of a Complex vector.
// Before calling, data must be allocated as Complex[nx*ny].
// On entry: data[ny*i+j] contains the ny Complex values for each i=0,...,nx-1;
//           log2nx contains the base-2 logarithm of nx;
//           log2ny contains the base-2 logarithm of ny;
//           isign is the sign (+/- 1) of the phase;
//           [scale is a constant by which the results will be multiplied;]
//           [bitreverse is 0 for a true fft of data (default);
//                         +1 for a fft of bit-reversed data (faster);
//                         -1 for a bit-reversed fft of data (faster).]
// On exit:  data[ny*i+j] contains the ny Complex Fourier values for
// each i=0,...nx-1.
// Note: To compute a true inverse transform, set scale=1.0/(nx*ny).

void fft2d(Complex *data, unsigned int log2nx, unsigned int log2ny, int isign,
	   Real scale, int bitreverse)
{
  unsigned int nx=1 << log2nx;
  unsigned int ny=1 << log2ny;
	
  mfft(data,log2nx,isign,ny,ny,1,scale,bitreverse);
  mfft(data,log2ny,isign,nx,1,ny,1.0,bitreverse);
}

// Return the two-dimensional Fourier transform of nx*ny real values.
// Before calling, data must be allocated as Complex[nx*(ny/2+1)].
// On entry:  ((Real *) data)[(ny+2)*i+j] must contain
// the (i,j)th real value, indexed by i=0,...,nx-1 and j=0,...,ny-1;
//           log2nx contains the base-2 logarithm of nx;
//           log2ny contains the base-2 logarithm of ny;
//           isign is the sign (+/- 1) of the phase;
//           [scale is a constant by which the results will be multiplied;]
//           [bitreverse is 0 for a true fft of data (default);
//                         +1 for a fft of bit-reversed data (faster);
//                         -1 for a bit-reversed fft of data (faster).]
// On exit: data[(ny/2+1)*i+j] contains the ny/2+1 Complex values for
// each i=0,...,nx-1. The origin of the Fourier domain is located at (nx/2,0).

void rcfft2d(Complex *data, unsigned int log2nx, unsigned int log2ny,
	     int isign, Real scale, int bitreverse)
{
  unsigned int nx=1 << log2nx;
  unsigned int ny=1 << log2ny;
  const unsigned int nyp=ny/2+1;

  Complex *pstop=data+nx*nyp;
  int pinc=2*nyp;
  for(Complex *p=data+nyp; p < pstop; p += pinc) {
    //#pragma ivdep
    for(unsigned int j=0; j < nyp; j++) p[j]=-p[j];
  }
	
  mrcfft(data,log2ny,isign,nx,1,nyp,scale,bitreverse);
  mfft(data,log2nx,isign,nyp,nyp,1);
}

// Return the two-dimensional real inverse Fourier transform of the
// nx*(ny/2+1) spectral values in the half-plane ky >= 0.
// Before calling, data must be allocated as Complex[nx*(ny/2+1)].
// On entry: data[(ny/2+1)*i+j] contains the ny/2+1 Complex values for
// each i=0,...,nx-1, 
//           log2nx contains the base-2 logarithm of nx;
//           log2ny contains the base-2 logarithm of ny;
//           isign is the sign (+/- 1) of the phase;
//           [scale is a constant by which the results will be multiplied;]
//           [bitreverse is 0 for a true fft of data (default);
//                         +1 for a fft of bit-reversed data (faster);
//                         -1 for a bit-reversed fft of data (faster).]
// The values corresponding to ky=0 are assumed to satisfy the reality
// condition (see crfft2d_sym below). The origin of the Fourier domain is
// located at (nx/2,0).
// On exit:  ((Real *) data)[(ny+2)*i+j] contains
// the (i,j)th real value, indexed by i=0,...,nx-1 and j=0,...,ny-1.
// Note: To compute a true inverse transform, set scale=1.0/(nx*ny).

void crfft2d(Complex *data, unsigned int log2nx, unsigned int log2ny,
	     int isign, Real scale, int bitreverse)
{
  unsigned int nx=1 << log2nx;
  unsigned int ny=1 << log2ny;
  const unsigned int nyp=ny/2+1;
	
  mfft(data,log2nx,isign,nyp,nyp,1);
  mcrfft(data,log2ny,isign,nx,1,nyp,scale,bitreverse);
	
  Complex *pstop=data+nx*nyp;
  int pinc=2*nyp;
  for(Complex *p=data+nyp; p < pstop; p += pinc) {
    //#pragma ivdep
    for(unsigned int j=0; j < nyp; j++) p[j]=-p[j];
  }
}

// Call crfft2d but first enforce reality condition by symmetrizing data.

void crfft2d_sym(Complex *data, unsigned int log2nx, unsigned int log2ny,
		 int isign, Real scale, int bitreverse)
{
  const unsigned int nx=1 << log2nx;
  const unsigned int ny=1 << log2ny;
  const unsigned int nx2=nx/2;
  const unsigned int nyp=ny/2+1;
	
  data[nyp*nx2].im=0.0;
  //#pragma ivdep
  for(unsigned int i=1; i < nx2; i++) data[nyp*i]=conj(data[nyp*(nx-i)]);
	
  crfft2d(data,log2nx,log2ny,isign,scale,bitreverse);
}

// Return the two-dimensional Fourier transform of nx*ny real values.
// Before calling, data must be allocated as Complex[nx*(ny/2+1)].
// On entry: ((Real *) data)[2*i+j+2*(nx-1)*(j/2)] must contain 
// the (i,j)th real value, indexed by i=0,...,nx-1 and j=0,...,ny-1;
//           log2nx contains the base-2 logarithm of nx;
//           log2ny contains the base-2 logarithm of ny;
//           isign is the sign (+/- 1) of the phase;
//           [scale is a constant by which the results will be multiplied;]
//           [bitreverse is 0 for a true fft of data (default);
//                         +1 for a fft of bit-reversed data (faster);
//                         -1 for a bit-reversed fft of data (faster).]
// On exit: data[i+nx*j] contains the nx Complex values for
// each j=0,...,ny/2. The origin of the Fourier domain is located at (nx/2,0).

void rcfft2dT(Complex *data, unsigned int log2nx, unsigned int log2ny,
	      int isign, Real scale, int bitreverse)
{
  unsigned int i;
  unsigned int nx=1 << log2nx;
  unsigned int ny=1 << log2ny;
  const unsigned int nyp=ny/2+1;
  const int nx1=nx+offset;
  Complex *p;

  mrcfft0(data,log2ny,isign,nx,nx1,1,scale,bitreverse);
	
  Complex *pstop=data+nx1*nyp;
#ifdef _CRAY
  for(i=1; i < nx; i += 2) {
    //#pragma ivdep
    for(p=data+i; p < pstop; p += nx1) *p=-(*p);
  }
#else	
  for(p=data; p < pstop; p += nx1) {
    for(i=1; i < nx; i += 2) p[i]=-p[i];
  }
#endif
	
  mfft(data,log2nx,isign,nyp,1,nx1);
}

// Return the two-dimensional real inverse Fourier transform of the
// nx*(ny/2+1) spectral values in the half-plane ky >= 0.
// Before calling, data must be allocated as Complex[nx*(ny/2+1)].
// On entry: data[i+nx*j] contains nx Complex values for each j=0,...,ny/2;
//           log2nx contains the base-2 logarithm of nx;
//           log2ny contains the base-2 logarithm of ny;
//           isign is the sign (+/- 1) of the phase;
//           [scale is a constant by which the results will be multiplied;]
//           [bitreverse is 0 for a true fft of data (default);
//                         +1 for a fft of bit-reversed data (faster);
//                         -1 for a bit-reversed fft of data (faster).]
// The values corresponding to ky=0 are assumed to satisfy the reality
// condition (see crfft2dT_sym below). The origin of the Fourier domain is
// located at (nx/2,0).
// On exit: ((Real *) data)[2*i+j+(nx-1)*(j/2)*2] contains 
// the (i,j)th real value, indexed by i=0,...,nx-1 and j=0,...,ny-1.
// Note: To compute a true inverse transform, set scale=1.0/(nx*ny).

void crfft2dT(Complex *data, unsigned int log2nx, unsigned int log2ny,
	      int isign, Real scale, int bitreverse)
{
  unsigned int nx=1 << log2nx;
  unsigned int ny=1 << log2ny;
  const unsigned int nyp=ny/2+1;
  const int nx1=nx+offset;

  mfft(data,log2nx,isign,nyp,1,nx1);
	
  Complex *p;
  Complex *pstop=data+nx1*nyp;
#ifdef _CRAY
  for(unsigned int i=1; i < nx; i += 2) {
    //#pragma ivdep
    for(p=data+i; p < pstop; p += nx1) *p=-(*p);
  }
#else	
  for(p=data; p < pstop; p += nx1) {
    for(unsigned int i=1; i < nx; i += 2) p[i]=-p[i];
  }
#endif	
	
  mcrfft0(data,log2ny,isign,nx,nx1,1,scale,bitreverse);
}

// Call crfft2dT but first enforce reality condition by symmetrizing data.

void crfft2dT_sym(Complex *data, unsigned int log2nx, unsigned int log2ny,
		  int isign, Real scale, int bitreverse)
{
  unsigned int nx=1 << log2nx;
  const unsigned int nx2=nx/2;
  data[nx2].im=0.0;
  //#pragma ivdep
  for(unsigned int i=1; i < nx2; i++) data[i]=conj(data[nx-i]);
	
  crfft2dT(data,log2nx,log2ny,isign,scale,bitreverse);
}

// Return the three-dimensional Fourier transform of a Complex vector.
// Before calling, data must be allocated as Complex[nx*ny*nz].
// On entry: data[nz*(ny*i+j)+k] contains the (i,j,k)th Complex value,
// indexed by i=0,...,nx-1, j=0,...,ny-1, and k=0,...,nz-1;
//           log2nx contains the base-2 logarithm of nx;
//           log2ny contains the base-2 logarithm of ny;
//           log2nz contains the base-2 logarithm of nz;
//           isign is the sign (+/- 1) of the phase;
//           [scale is a constant by which the results will be multiplied;]
//           [bitreverse is 0 for a true fft of data (default);
//                         +1 for a fft of bit-reversed data (faster);
//                         -1 for a bit-reversed fft of data (faster).]
// On exit:  data[nz*(ny*i+j)+k] contains the (i,j,k)th Complex Fourier
// values indexed by i=0,...,nx-1, j=0,...,ny-1, and k=0,...,nz-1.
// Note: To compute a true inverse transform, set scale=1.0/(nx*ny*nz).

void fft3d(Complex *data, unsigned int log2nx, unsigned int log2ny,
	   unsigned int log2nz, int isign, Real scale, int bitreverse)
{
  unsigned int nx=1 << log2nx;
  unsigned int ny=1 << log2ny;
  unsigned int nz=1 << log2nz;
  int nyz=ny*nz;
	
  mfft(data,log2nx,isign,nyz,nyz,1,scale,bitreverse);
  for(unsigned int i=0; i < nx; i++)
    mfft(data+i*nyz,log2ny,isign,nz,nz,1,1.0,bitreverse);
  mfft(data,log2nz,isign,nx*ny,1,nz,1.0,bitreverse);
}

// Return the three-dimensional Fourier transform of nx*ny*nz real values.
// Before calling, data must be allocated as Complex[nx*ny*(nz/2+1)].
// On entry:  ((Real *) data)[(nz+2)*(ny*i+j)+k] must contain
// the (i,j,k)th real value, indexed by i=0,...,nx-1, j=0,...,ny-1, and
// k=0,...,nz-1;
//           log2nx contains the base-2 logarithm of nx;
//           log2ny contains the base-2 logarithm of ny;
//           log2nz contains the base-2 logarithm of nz;
//           isign is the sign (+/- 1) of the phase;
//           [scale is a constant by which the results will be multiplied;]
//           [bitreverse is 0 for a true fft of data (default);
//                         +1 for a fft of bit-reversed data (faster);
//                         -1 for a bit-reversed fft of data (faster).]
// On exit: data[(nz/2+1)*(ny*i+j)+k] contains the nz/2+1 Complex values for
// each i=0,...,nx-1 and j=0,...,ny-1. The origin of the Fourier domain is
// located at (nx/2,ny/2,0).

void rcfft3d(Complex *data, unsigned int log2nx, unsigned int log2ny,
	     unsigned int log2nz, int isign, Real scale, int bitreverse)
{
  unsigned int i;
  unsigned int nx=1 << log2nx;
  unsigned int ny=1 << log2ny;
  unsigned int nz=1 << log2nz;
  const unsigned int nzp=nz/2+1;
  Complex *p,*pstop;

  mrcfft(data,log2nz,isign,nx*ny,1,nzp,scale,bitreverse);
	
  int nyzp=ny*nzp;
  int pinc=2*nzp;
  p=pstop=data;
  for(i=0; i < nx; i++) {
    if(i % 2) p -= nzp;
    else p += nzp;
    pstop += nyzp;
    for(; p < pstop; p += pinc) {
      //#pragma ivdep
      for(unsigned int k=0; k < nzp; k++) p[k]=-p[k];
    }
  }
	
  mfft(data,log2nx,isign,nyzp,nyzp,1);
	
  for(i=0; i < nx; i++)
    mfft(data+i*nyzp,log2ny,isign,nzp,nzp,1);
}

// Return the three-dimensional real inverse Fourier transform of the
// nx*ny*(nz/2+1) spectral values in the half-space kz >= 0.
// Before calling, data must be allocated as Complex[nx*ny*(nz/2+1)].
// On entry: data[(nz/2+1)*(ny*i+j)+k] must contain the nz/2+1 Complex
// values for each i=0,...,nx-1 and j=0,...,ny-1; 
//           log2nx contains the base-2 logarithm of nx;
//           log2ny contains the base-2 logarithm of ny;
//           log2nz contains the base-2 logarithm of nz;
//           isign is the sign (+/- 1) of the phase;
//           [scale is a constant by which the results will be multiplied;]
//           [bitreverse is 0 for a true fft of data (default);
//                         +1 for a fft of bit-reversed data (faster);
//                         -1 for a bit-reversed fft of data (faster).]
// The values corresponding to kz=0 are assumed to satisfy the reality
// condition (see crfft3d_sym below). The origin of the Fourier domain is
// located at (nx/2,ny/2,0).
// On exit:  ((Real *) data)[(nz+2)*(ny*i+j)+k] contains the (i,j,k)th real 
// value, indexed by i=0,...,nx-1, j=0,...,ny-1, and k=0,...,nz-1.
// Note: To compute a true inverse transform, set scale=1.0/(nx*ny*nz).

void crfft3d(Complex *data, unsigned int log2nx, unsigned int log2ny,
	     unsigned int log2nz, int isign, Real scale, int bitreverse)
{
  unsigned int i;
  unsigned int nx=1 << log2nx;
  unsigned int ny=1 << log2ny;
  unsigned int nz=1 << log2nz;
  const unsigned int nzp=nz/2+1;
  Complex *p,*pstop;

  const int unsigned nyzp=ny*nzp;
	
  for(i=0; i < nx; i++)
    mfft(data+i*nyzp,log2ny,isign,nzp,nzp,1);
	
  mfft(data,log2nx,isign,nyzp,nyzp,1);
	
  int pinc=2*nzp;
  p=pstop=data;
  for(i=0; i < nx; i++) {
    if(i % 2) p -= nzp;
    else p += nzp;
    pstop += nyzp;
    for(; p < pstop; p += pinc) {
      //#pragma ivdep
      for(unsigned int k=0; k < nzp; k++) p[k]=-p[k];
    }
  }
	
  mcrfft(data,log2nz,isign,nx*ny,1,nzp,scale,bitreverse);
}

// Call crfft3d but first enforce reality condition by symmetrizing data.

void crfft3d_sym(Complex *data, unsigned int log2nx, unsigned int log2ny,
		 unsigned int log2nz, int isign, Real scale, int bitreverse)
{
  const unsigned int nx=1 << log2nx;
  const unsigned int ny=1 << log2ny;
  const unsigned int nz=1 << log2nz;
  const unsigned int nx2=nx/2;
  const unsigned int nzp=nz/2+1;
  const unsigned int nyzp=ny*nzp;
	
  data[nyzp*nx2].im=0.0;
  for(unsigned int i=1; i < nx2; i++) {
    data[nyzp*i]=conj(data[nyzp*(nx-i)]);
    //#pragma ivdep
    for(unsigned int j=1; j < ny; j++)
      data[nzp*(ny*i+j)]=conj(data[nzp*(ny*(nx-i)+(ny-j))]);
  }
	
  crfft3d(data,log2nx,log2ny,log2nz,isign,scale,bitreverse);
}

