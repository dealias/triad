#include "options.h"
#include "fft.h"

Complex *wpTable;
unsigned int wpTableSize=0;

static Complex *WTable;
static unsigned int WTableSize=0;

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
// length n/2.
//           log2n contains the base-2 logarithm of n.
//           isign is +1 for a forward transform, -1 for an inverse transform.
//           bitreverse is 0 for a true fft, 1 for a faster bit-reversed fft.
// On exit:  data contains the n/2+1 Complex Fourier values.

void rcfft(Complex *data, unsigned int log2n, int isign, int bitreverse)
{		 
	log2n--;
	unsigned int n2=1 << log2n;
	unsigned int n4=n2 >> 1;
	if(WTableSize != n4) rfft_init(log2n);
	
	unsigned int log4n=log2n/2;
	// Special case for powers of 4.
	if(2*log4n == log2n) fft4(data,log4n,1);
	else fft(data,log2n,1,bitreverse);
	
	data[n2]=data[0].re-data[0].im;
	data[0]=data[0].re+data[0].im;
	if(isign == 1) {
#pragma ivdep	
		for(unsigned int i=1; i < n4; i++) {
			Complex u=data[i], v=conj(data[n2-i]);
			Complex A=0.5*(u+v), B=WTable[i]*(u-v);
			data[i]=A-B;
			data[n2-i]=conj(A+B);
		}
	} else {
#pragma ivdep	
		for(unsigned int i=1; i < n4; i++) {
			Complex u=data[i], v=conj(data[n2-i]);
			Complex A=0.5*(u+v), B=WTable[i]*(u-v);
			data[i]=conj(A-B);
			data[n2-i]=A+B;
		}
		data[n4]=conj(data[n4]);
	}
}

// Return the real inverse Fourier transform of the n/2+1
// Complex values corresponding to the non-negative part of the frequency
// spectrum. Before calling, data must be allocated as Complex[n/2+1].
// On entry: data contains the n/2+1 Complex Fourier transform values.
//           log2n contains the base-2 logarithm of n.
//           isign is +1 for a forward transform, -1 for an inverse transform.
//           bitreverse is 0 for a true fft, 1 for a faster bit-reversed fft.
// On exit:  data contains the n real inverse Fourier transform values
// stored as a Complex array of length n/2.
// Note: The final result must be divided by n.

void crfft(Complex *data, unsigned int log2n, int isign, int bitreverse)
{		 
	log2n--;
	unsigned int n2=1 << log2n;
	unsigned int n4=n2 >> 1;
	if(WTableSize != n4) rfft_init(log2n);
	
	data[0].im=data[0].re-data[n2].re;
	data[0].re += data[n2].re;
	
	if(isign == 1) {
#pragma ivdep
		for(unsigned int i=1; i < n4; i++) {
			Complex u=conj(data[i]), v=data[n2-i];
			Complex A=u+v, B=2.0*conj(WTable[i])*(u-v);
			data[i]=A-B;
			data[n2-i]=conj(A+B);
		}
		data[n4]=2.0*conj(data[n4]);
	} else {
#pragma ivdep
		for(unsigned int i=1; i < n4; i++) {
			Complex u=data[i], v=conj(data[n2-i]);
			Complex A=u+v, B=2.0*conj(WTable[i])*(u-v);
			data[i]=A-B;
			data[n2-i]=conj(A+B);
		}
		data[n4] *= 2.0;
	}
	
	unsigned int log4n=log2n/2;
	// Special case for powers of 4.
	if(2*log4n == log2n) fft4(data,log4n,-1);
	else fft(data,log2n,-1,bitreverse);
}

// Return the real inverse Fourier transform of nk Complex vectors, of
// length n/2+1, corresponding to the non-negative part of the frequency
// spectrum. Before calling, data must be allocated as Complex[nk*(n/2+1)].
// On entry: data contains the n/2+1 Complex Fourier transform values for
// each k=0,...,nk-1. 
//           log2n contains the base-2 logarithm of n.
//           isign is +1 for a forward transform, -1 for an inverse transform.
//           inc1 is the stride between the elements of each Complex vector.
//           inc2 is the stride between first elements of the vectors.
//           bitreverse is 0 for a true fft, 1 for a faster bit-reversed fft.
// On exit:  data contains the n real inverse Fourier transform values
// stored as a Complex array of length n/2.
// Note: The final result must be divided by n.

void mcrfft(Complex *data, unsigned int log2n, int isign, unsigned int nk,
			unsigned int inc1, unsigned int inc2, int bitreverse)
{		 
	log2n--;
	if(inc1 == 0) inc1=nk;
	unsigned int kstop=nk*inc2;
	
	unsigned int n2=1 << log2n;
	unsigned int n4=n2 >> 1;
	if(WTableSize != n4) rfft_init(log2n);
	
	Complex *q=data+n2*inc1;
	for(unsigned int k=0; k < kstop; k += inc2) {
		data[k].im=data[k].re-q[k].re;
		data[k].re += q[k].re;
	}
	
	if(isign == 1) {
#pragma ivdep
		for(unsigned int i=1; i < n4; i++) {
			Complex *p=data+i*inc1, *q=data+(n2-i)*inc1;
			Complex W=conj(WTable[i]);
			for(unsigned int k=0; k < kstop; k += inc2) {
				Complex u=conj(p[k]), v=q[k];
				Complex A=u+v, B=2.0*W*(u-v);
				p[k]=A-B;
				q[k]=conj(A+B);
			}
		}
		Complex *p=data+n4*inc1;
		for(unsigned int k=0; k < kstop; k += inc2) p[k]=2.0*conj(p[k]);
	} else {
#pragma ivdep
		for(unsigned int i=1; i < n4; i++) {
			Complex *p=data+i*inc1, *q=data+(n2-i)*inc1;
			Complex W=conj(WTable[i]);
			for(unsigned int k=0; k < kstop; k += inc2) {
				Complex u=p[k], v=conj(q[k]);
				Complex A=u+v, B=2.0*W*(u-v);
				p[k]=A-B;
				q[k]=conj(A+B);
			}
		}
		Complex *p=data+n4*inc1;
		for(unsigned int k=0; k < kstop; k += inc2) p[k] *= 2.0;
	}
	
	mfft(data,log2n,-1,nk,inc1,inc2,bitreverse);
}

// Return the two-dimensional real inverse Fourier transform of the
// nx*(ny/2+1) spectral values taken from the positive frequency half-plane.
// Before calling, data must be allocated as Complex[nx*(ny/2+1)].
// On entry: data[i+nx*j] contains the nx Complex values for
// each j=0,...,ny/2. 
//           log2nx contains the base-2 logarithm of nx.
//           log2ny contains the base-2 logarithm of ny.
// On exit: ((Real *) data)[2*i+j+(nx-1)*(j/2)*2] contains 
// the (i,j)th real value, indexed by i=0,...,nx-1 and j=0,...,ny-1.
// Note: The final result must be divided by nx*ny.

void crfft2dT(Complex *data, unsigned int log2nx, unsigned int log2ny,
			 int isign)
{
	unsigned int i,j;
	unsigned int nx=1 << log2nx;
	unsigned int ny=1 << log2ny;
	const unsigned int nyp=ny/2+1;

	for(j=0; j < nyp; j++) {
		Complex *p=data+nx*j;
		for(i=1; i < nx; i += 2) p[i] *= -1.0;
	}

	mfft(data,log2nx,isign,nyp,1,nx);
	
	for(j=0; j < nyp; j++) {
		Complex *p=data+nx*j;
		for(i=1; i < nx; i += 2) p[i] *= -1.0;
	}
	
	for(j=1; j < nyp; j += 2) {
		Complex *p=data+nx*j;
		for(i=0; i < nx; i++) p[i] *= -1.0;
	}

	mcrfft(data,log2ny,isign,nx);
}

// Return the two-dimensional real inverse Fourier transform of the
// nx*(ny/2+1) spectral values taken from the positive frequency half-plane.
// Before calling, data must be allocated as Complex[nx*(ny/2+1)].
// On entry: data[(ny/2+1)*i+j] contains the (ny/2+1) Complex values for
// each i=0,...,nx-1. 
//           log2nx contains the base-2 logarithm of nx.
//           log2ny contains the base-2 logarithm of ny.
// On exit:  ((Real *) data)[(ny+2)*i+j] contains
// the (i,j)th real value, indexed by i=0,...,nx-1 and j=0,...,ny-1.
// Note: The final result must be divided by nx*ny.

void crfft2d(Complex *data, unsigned int log2nx, unsigned int log2ny,
			 int isign)
{
	unsigned int i,j;
	unsigned int nx=1 << log2nx;
	unsigned int ny=1 << log2ny;
	const unsigned int nyp=ny/2+1;

	for(i=1; i < nx; i += 2) {
		Complex *p=data+i*nyp;
		for(j=0; j < nyp; j++) p[j] *= -1.0;
	}

	mfft(data,log2nx,-1,nyp);
	
	for(i=1; i < nx; i += 2) {
		Complex *p=data+i*nyp;
		for(j=0; j < nyp; j++) p[j] *= -1.0;
	}
	
	for(i=0; i < nx; i++) {
		Complex *p=data+i*nyp;
		for(j=1; j < nyp; j += 2) p[j] *= -1.0;
	}
	
	mcrfft(data,log2ny,isign,nx,1,nyp);
}
