#include "options.h"
#include "fft.h"

Complex *wpTable;
unsigned int wpTableSize=0;

static Complex *WTable;
static unsigned int WTableSize=0;

#if _CRAY
const int offset=1;
#else
const int offset=0;
#endif	

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
		data[n4].im=-data[n4].im;
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

// Return the Fourier transform of nk real vectors, each of length n.
// Before calling, data must be allocated as Complex[nk*(n/2+1)].
// On entry: data contains the n real values stored as a Complex array of
// length n/2, for each k=0,...,nk-1.
//           log2n contains the base-2 logarithm of n.
//           isign is +1 for a forward transform, -1 for an inverse transform.
//           inc1 is the stride between the elements of each Complex vector.
//           inc2 is the stride between first elements of the vectors.
//           bitreverse is 0 for a true fft, 1 for a faster bit-reversed fft.
// On exit:  data contains the n/2+1 Complex Fourier values.

void mrcfft(Complex *data, unsigned int log2n, int isign, unsigned int nk,
			unsigned int inc1, unsigned int inc2, int bitreverse)
{		 
	log2n--;
	if(inc1 == 0) inc1=nk;
	unsigned int kstop=nk*inc2;
	
	unsigned int n2=1 << log2n;
	unsigned int n4=n2 >> 1;
	if(WTableSize != n4) rfft_init(log2n);
	
	mfft(data,log2n,1,nk,inc1,inc2,bitreverse);
	
	Complex *q=data+n2*inc1;
#pragma ivdep	
	for(unsigned int k=0; k < kstop; k += inc2) {
		q[k]=data[k].re-data[k].im;
		data[k]=data[k].re+data[k].im;
	}
	
	if(isign == 1) {
		for(unsigned int i=1; i < n4; i++) {
			Complex *p=data+i*inc1, *q=data+(n2-i)*inc1;
			Complex W=WTable[i];
#pragma ivdep	
			for(unsigned int k=0; k < kstop; k += inc2) {
				Complex u=p[k], v=conj(q[k]);
				Complex A=0.5*(u+v), B=W*(u-v);
				p[k]=A-B;
				q[k]=conj(A+B);
			}
		}
	} else {
		for(unsigned int i=1; i < n4; i++) {
			Complex *p=data+i*inc1, *q=data+(n2-i)*inc1;
			Complex W=WTable[i];
#pragma ivdep	
			for(unsigned int k=0; k < kstop; k += inc2) {
				Complex u=p[k], v=conj(q[k]);
				Complex A=0.5*(u+v), B=W*(u-v);
				p[k]=conj(A-B);
				q[k]=A+B;
			}
		}
		Complex *p=data+n4*inc1;
#pragma ivdep	
		for(unsigned int k=0; k < kstop; k += inc2) p[k].im=-p[k].im;
	}
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
#pragma ivdep	
	for(unsigned int k=0; k < kstop; k += inc2) {
		data[k].im=data[k].re-q[k].re;
		data[k].re += q[k].re;
	}
	
	if(isign == 1) {
		for(unsigned int i=1; i < n4; i++) {
			Complex *p=data+i*inc1, *q=data+(n2-i)*inc1;
			Complex W=2.0*conj(WTable[i]);
#pragma ivdep
			for(unsigned int k=0; k < kstop; k += inc2) {
				Complex u=conj(p[k]), v=q[k];
				Complex A=u+v, B=W*(u-v);
				p[k]=A-B;
				q[k]=conj(A+B);
			}
		}
		Complex *p=data+n4*inc1;
#pragma ivdep	
		for(unsigned int k=0; k < kstop; k += inc2) p[k]=2.0*conj(p[k]);
	} else {
		for(unsigned int i=1; i < n4; i++) {
			Complex *p=data+i*inc1, *q=data+(n2-i)*inc1;
			Complex W=2.0*conj(WTable[i]);
#pragma ivdep
			for(unsigned int k=0; k < kstop; k += inc2) {
				Complex u=p[k], v=conj(q[k]);
				Complex A=u+v, B=W*(u-v);
				p[k]=A-B;
				q[k]=conj(A+B);
			}
		}
		Complex *p=data+n4*inc1;
#pragma ivdep
		for(unsigned int k=0; k < kstop; k += inc2) p[k] *= 2.0;
	}
	
	mfft(data,log2n,-1,nk,inc1,inc2,bitreverse);
}

// Return the two-dimensional Fourier transform of nx*ny real values.
// Before calling, data must be allocated as Complex[nx*(ny/2+1)].
// On entry: ((Real *) data)[2*i+j+2*(nx-1)*(j/2)] must contain 
// the (i,j)th real value, indexed by i=0,...,nx-1 and j=0,...,ny-1.
//           log2nx contains the base-2 logarithm of nx.
//           log2ny contains the base-2 logarithm of ny.
//           isign is +1 for a forward transform, -1 for an inverse transform.
// On exit: data[i+nx*j] contains the nx Complex values for
// each j=0,...,ny/2.
// The origin of the Fourier domain is located at (nx/2,0).

void rcfft2dT(Complex *data, unsigned int log2nx, unsigned int log2ny,
			  int isign)
{
	unsigned int i;
	unsigned int nx=1 << log2nx;
	unsigned int ny=1 << log2ny;
	const unsigned int nyp=ny/2+1;
	const int nx1=nx+offset;
	Complex *p;

	mrcfft(data,log2ny,isign,nx,nx1,1);
	
	Complex *pstop=data+nx1*nyp;
#if _CRAY
	for(i=1; i < nx; i += 2) {
#pragma ivdep
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
// nx*(ny/2+1) spectral values taken from the frequency half-plane.
// Before calling, data must be allocated as Complex[nx*(ny/2+1)].
// On entry: data[i+nx*j] contains nx/2+1 Complex values for j=0
// and nx values for j=1,...,ny/2.
// The origin of the Fourier domain is located at (nx/2,0).
//           log2nx contains the base-2 logarithm of nx.
//           log2ny contains the base-2 logarithm of ny.
//           isign is +1 for a forward transform, -1 for an inverse transform.
// On exit: ((Real *) data)[2*i+j+(nx-1)*(j/2)*2] contains 
// the (i,j)th real value, indexed by i=0,...,nx-1 and j=0,...,ny-1.
// Note: The final result must be divided by nx*ny.

void crfft2dT(Complex *data, unsigned int log2nx, unsigned int log2ny,
			  int isign)
{
	unsigned int i;
	unsigned int nx=1 << log2nx;
	unsigned int ny=1 << log2ny;
	const unsigned int nyp=ny/2+1;
	const unsigned int nx2=nx/2;
	const int nx1=nx+offset;
	Complex *p;

	// Enforce reality condition
	data[0].im=0.0;
#pragma ivdep
	for(i=1; i < nx2; i++) data[i]=conj(data[nx-i]);
	
	mfft(data,log2nx,isign,nyp,1,nx1);

	Complex *pstop=data+nx1*nyp;
#if _CRAY
	for(i=1; i < nx; i += 2) {
#pragma ivdep
		for(p=data+i; p < pstop; p += nx1) *p=-(*p);
	}
#else	
	for(p=data; p < pstop; p += nx1) {
		for(i=1; i < nx; i += 2) p[i]=-p[i];
	}
#endif	
	
	mcrfft(data,log2ny,isign,nx,nx1,1);
}

// Return the two-dimensional Fourier transform of nx*ny real values.
// Before calling, data must be allocated as Complex[nx*(ny/2+1)].
// On entry:  ((Real *) data)[(ny+2)*i+j] must contain
// the (i,j)th real value, indexed by i=0,...,nx-1 and j=0,...,ny-1.
//           log2nx contains the base-2 logarithm of nx.
//           log2ny contains the base-2 logarithm of ny.
//           isign is +1 for a forward transform, -1 for an inverse transform.
// On exit: data[(ny/2+1)*i+j] contains the ny/2+1 Complex values for
// each i=0,...,nx-1. 
// The origin of the Fourier domain is located at (nx/2,0).

void rcfft2d(Complex *data, unsigned int log2nx, unsigned int log2ny,
			 int isign)
{
	unsigned int nx=1 << log2nx;
	unsigned int ny=1 << log2ny;
	const unsigned int nyp=ny/2+1;

	mrcfft(data,log2ny,isign,nx,1,nyp);

	Complex *pstop=data+nx*nyp;
	int pinc=2*nyp;
	for(Complex *p=data+nyp; p < pstop; p += pinc) {
#pragma ivdep
		for(unsigned int j=0; j < nyp; j++) p[j]=-p[j];
	}
	
	mfft(data,log2nx,isign,nyp);
}

// Return the two-dimensional real inverse Fourier transform of the
// nx*(ny/2+1) spectral values taken from the non-negative frequency
// half-plane. 
// Before calling, data must be allocated as Complex[nx*(ny/2+1)].
// On entry: data[(ny/2+1)*i+j] contains the ny/2+1 Complex values for
// each i=0,...,nx-1. 
// The origin of the Fourier domain is located at (nx/2,0).
//           log2nx contains the base-2 logarithm of nx.
//           log2ny contains the base-2 logarithm of ny.
//           isign is +1 for a forward transform, -1 for an inverse transform.
// On exit:  ((Real *) data)[(ny+2)*i+j] contains
// the (i,j)th real value, indexed by i=0,...,nx-1 and j=0,...,ny-1.
// Note: The final result must be divided by nx*ny.

void crfft2d(Complex *data, unsigned int log2nx, unsigned int log2ny,
			 int isign)
{
	unsigned int nx=1 << log2nx;
	unsigned int ny=1 << log2ny;
	const unsigned int nyp=ny/2+1;
	const unsigned int nx2=nx/2;

	// Enforce reality condition
	data[0].im=0.0;
#pragma ivdep
	for(unsigned int i=1; i < nx2; i++) data[nyp*i]=conj(data[nyp*(nx-i)]);
	
	mfft(data,log2nx,isign,nyp);
	
	Complex *pstop=data+nx*nyp;
	int pinc=2*nyp;
	for(Complex *p=data+nyp; p < pstop; p += pinc) {
#pragma ivdep
		for(unsigned int j=0; j < nyp; j++) p[j]=-p[j];
	}
	
	mcrfft(data,log2ny,isign,nx,1,nyp);
}

// Return the three-dimensional Fourier transform of nx*ny*nz real values.
// Before calling, data must be allocated as Complex[nx*ny*(nz/2+1)].
// On entry:  ((Real *) data)[(nz+2)*(ny*i+j)+k] must contain
// the (i,j,k)th real value, indexed by i=0,...,nx-1, j=0,...,ny-1, and
// k=0,...,nz-1.
//           log2nx contains the base-2 logarithm of nx.
//           log2ny contains the base-2 logarithm of ny.
//           log2nz contains the base-2 logarithm of nz.
//           isign is +1 for a forward transform, -1 for an inverse transform.
// On exit: data[(nz/2+1)*(ny*i+j)+k] contains the nz/2+1 Complex values for
// each i=0,...,nx-1 and j=0,...,ny-1. 
// The origin of the Fourier domain is located at (nx/2,ny/2,0).

void rcfft3d(Complex *data, unsigned int log2nx, unsigned int log2ny,
			 unsigned int log2nz, int isign)
{
	unsigned int i;
	unsigned int nx=1 << log2nx;
	unsigned int ny=1 << log2ny;
	unsigned int nz=1 << log2nz;
	const unsigned int nzp=nz/2+1;
	Complex *p,*pstop;

	mrcfft(data,log2nz,isign,nx*ny,1,nzp);
	
	int nyzp=ny*nzp;
	int pinc=2*nzp;
	p=pstop=data;
	for(i=0; i < nx; i++) {
		if(i % 2) p -= nzp;
		else p += nzp;
		pstop += nyzp;
		for(; p < pstop; p += pinc) {
#pragma ivdep
			for(unsigned int k=0; k < nzp; k++) p[k]=-p[k];
		}
	}
	
	mfft(data,log2nx,isign,nyzp);
	
	for(i=0; i < nx; i++)
		mfft(data+i*nyzp,log2ny,isign,nzp);
}

// On exit:  ((Real *) data)[(ny+2)*i+j] contains
// the (i,j)th real value, indexed by i=0,...,nx-1 and j=0,...,ny-1.
// Note: The final result must be divided by nx*ny.

// Return the three-dimensional real inverse Fourier transform of the
// nx*ny*(nz/2+1) spectral values taken from the non-negative frequency
// half-plane. 
// Before calling, data must be allocated as Complex[nx*ny*(nz/2+1)].
// On entry: data[(nz/2+1)*(ny*i+j)+k] must contain the nz/2+1 Complex
// values for each i=0,...,nx-1 and j=0,...,ny-1. 
// The origin of the Fourier domain is located at (nx/2,ny/2,0).
//           log2nx contains the base-2 logarithm of nx.
//           log2ny contains the base-2 logarithm of ny.
//           log2nz contains the base-2 logarithm of nz.
//           isign is +1 for a forward transform, -1 for an inverse transform.
// On exit:  ((Real *) data)[(nz+2)*(ny*i+j)+k] contains the (i,j,k)th real
// value, indexed by i=0,...,nx-1, j=0,...,ny-1, and k=0,...,nz-1.
// Note: The final result must be divided by nx*ny*nz.

void crfft3d(Complex *data, unsigned int log2nx, unsigned int log2ny,
			 unsigned int log2nz, int isign)
{
	unsigned int i;
	unsigned int nx=1 << log2nx;
	unsigned int ny=1 << log2ny;
	unsigned int nz=1 << log2nz;
	const unsigned int nzp=nz/2+1;
	const unsigned int nx2=nx/2;
	Complex *p,*pstop;

	int nyzp=ny*nzp;
	
	// Enforce reality condition
	data[0].im=0.0;
	for(i=1; i < nx2; i++) {
		data[nyzp*i]=conj(data[nyzp*(nx-i)]);
#pragma ivdep
		for(unsigned int j=1; j < ny; j++)
			data[nzp*(ny*i+j)]=conj(data[nzp*(ny*(nx-i)+(ny-j))]);
	}
	
	for(i=0; i < nx; i++)
		mfft(data+i*nyzp,log2ny,isign,nzp);
		
	mfft(data,log2nx,isign,nyzp);
	
	int pinc=2*nzp;
	p=pstop=data;
	for(i=0; i < nx; i++) {
		if(i % 2) p -= nzp;
		else p += nzp;
		pstop += nyzp;
		for(; p < pstop; p += pinc) {
#pragma ivdep
			for(unsigned int k=0; k < nzp; k++) p[k]=-p[k];
		}
	}
		
	mcrfft(data,log2nz,isign,nx*ny,1,nzp);
}

// Return the Fourier transform of a Complex vector of length n=power of 4.
// Before calling, data must be allocated as Complex[n].
// On entry: data contains the n Complex values.
//           log4n contains the base-4 logarithm of n.
//           isign is +1 for a forward transform, -1 for an inverse transform.
// On exit:  data contains the n Complex Fourier values.
// Note: When computing an inverse transform, the result must be divided by n.

void fft4(Complex *data, unsigned int log4n, int isign)
{
	unsigned int m,k,j;
	static unsigned int lastsize=0;
	static Complex *phase;

	m=1 << log4n;
	
	if(m != lastsize) {
		if(m > lastsize) phase=new(phase,m*m) Complex;
		lastsize=m;
		Complex factor=expi(twopi/(m*m));
		Complex kphase=1.0;
		for(k=0; k < m; k++) {
			Complex *p=phase+m*k;
			p[0]=1.0;
			for(j=1; j < m; j++) p[j]=p[j-1]*kphase;
			kphase *= factor;
		}
	}
	
	mfft(data,log4n,isign,m);
	
	if(isign == 1) {
		for(j=0; j < m; j++) {
			int mj=m*j;
			Complex *p=data+mj;
			Complex *c=phase+mj;
			Complex *q=data+j;
			p[j] *= c[j];
			for(k=0; k < j; k++, q += m) {
				Complex temp = p[k];
				p[k] = (*q)*c[k];
				(*q)=temp*c[k];
			}
		}
	} else {
		for(j=0; j < m; j++) {
			int mj=m*j;
			Complex *p=data+mj;
			Complex *c=phase+mj;
			Complex *q=data+j;
			p[j] *= conj(c[j]);
			for(k=0; k < j; k++, q += m) {
				Complex temp = p[k];
				p[k] = (*q)*conj(c[k]);
				(*q)=temp*conj(c[k]);
			}
		}
	}
	
	mfft(data,log4n,isign,m);
}

// Return the two-dimensional Fourier transform of a Complex vector.
// Before calling, data must be allocated as Complex[nx*ny].
// On entry: data[ny*i+j] contains the ny Complex values for each i=0,...,nx-1.
//           log2nx contains the base-2 logarithm of nx.
//           log2ny contains the base-2 logarithm of ny.
//           isign is +1 for a forward transform, -1 for an inverse transform.
// On exit:  data[ny*i+j] contains the ny Complex Fourier values for
// each i=0,...nx-1.
// Note: When computing an inverse transform, the result must be divided
// by nx*ny.

void fft2d(Complex *data, unsigned int log2nx, unsigned int log2ny, int isign)
{
	unsigned int nx=1 << log2nx;
	unsigned int ny=1 << log2ny;
	
	mfft(data,log2nx,isign,ny);
	mfft(data,log2ny,isign,nx,1,ny);
}

// Return the three-dimensional Fourier transform of a Complex vector.
// Before calling, data must be allocated as Complex[nx*ny*nz].
// On entry: data[nz*(ny*i+j)+k] contains the (i,j,k)th Complex value,
// indexed by i=0,...,nx-1, j=0,...,ny-1, and k=0,...,nz-1.
//           log2nx contains the base-2 logarithm of nx.
//           log2ny contains the base-2 logarithm of ny.
//           log2nz contains the base-2 logarithm of nz.
//           isign is +1 for a forward transform, -1 for an inverse transform.
// On exit:  data[nz*(ny*i+j)+k] contains the (i,j,k)th Complex Fourier
// values indexed by i=0,...,nx-1, j=0,...,ny-1, and k=0,...,nz-1.
// Note: When computing an inverse transform, the result must be divided
// by nx*ny*nz.

void fft3d(Complex *data, unsigned int log2nx, unsigned int log2ny,
		   unsigned int log2nz, int isign)
{
	unsigned int nx=1 << log2nx;
	unsigned int ny=1 << log2ny;
	unsigned int nz=1 << log2nz;
	
	mfft(data,log2nx,isign,ny*nz);
	int nyz=ny*nz;
	for(unsigned int i=0; i < nx; i++)
		mfft(data+i*nyz,log2ny,isign,nz);
	mfft(data,log2nz,isign,nx*ny,1,nz);
}
