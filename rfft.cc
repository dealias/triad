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
	if(log2n > wpTableSize) fft_init(log2n);
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
//           bitreverse is 0 for a true fft, 1 for a faster bit-reversed fft.
// On exit:  data contains the n/2+1 Complex Fourier values.

void rfft(Complex *data, unsigned int log2n, int bitreverse)
{		 
	unsigned int i;
	if(log2n > wpTableSize) fft_init(log2n+1);
	
	unsigned int log4n=log2n/2;
	if(2*log4n == log2n) fft4(data,log4n,1); // Special case for powers of 4.
	else fft(data,log2n,1,bitreverse);
	
	unsigned int n2=1 << log2n;
	unsigned int n4=n2 >> 1;
	if(WTableSize != n4) rfft_init(log2n);
	
	data[n2]=data[0].re-data[0].im;
	data[0]=data[0].re+data[0].im;
#pragma ivdep	
	for(i=1; i < n4; i++) {
		Complex u=data[i], v=conj(data[n2-i]);
		Complex A=0.5*(u+v), B=WTable[i]*(u-v);
		data[i]=A-B;
		data[n2-i]=conj(A+B);
	}
}


// Return the real inverse Fourier transform of the n/2+1
// Complex values corresponding to the non-negative part of the frequency
// spectrum. Before calling, data must be allocated as Complex[n/2+1].
// On entry: data contains the n/2+1 Complex Fourier transform values.
//           log2n contains the base-2 logarithm of n.
//           bitreverse is 0 for a true fft, 1 for a faster bit-reversed fft.
// On exit:  data contains the n real inverse Fourier transform values
// stored as a Complex array of length n/2.

void rfft_inv(Complex *data, unsigned int log2n, int bitreverse)
{		 
	unsigned int i;
	if(log2n > wpTableSize) fft_init(log2n+1);
	
	unsigned int n2=1 << log2n;
	unsigned int n4=n2 >> 1;
	if(WTableSize != n4) rfft_init(log2n);
	
	data[0].im=0.5*(data[0].re-data[n2].re);
	data[0].re=0.5*(data[0].re+data[n2].re);
#pragma ivdep	
	for(i=1; i < n4; i++) {
		Complex u=data[i], v=conj(data[n2-i]);
		Complex A=0.5*(u+v), B=conj(WTable[i])*(u-v);
		data[i]=A-B;
		data[n2-i]=conj(A+B);
	}
	
	unsigned int log4n=log2n/2;
	if(2*log4n == log2n) fft4(data,log4n,-1); // Special case for powers of 4.
	else fft(data,log2n,-1,bitreverse);
}
