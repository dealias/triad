#include "utils.h"

Complex *wpTablep,*wpTablen;
static Complex *WTablep,*WTablen;
static unsigned int TableSize=0,WTablepSize=0,WTablenSize=0;

void fft_br(Complex *data, unsigned int log2n);
void fft_brinv(Complex *data, unsigned int log2n);
	
void rfft_br(Complex *data, unsigned int log2n)
{		 
	unsigned int i;
	
	fft_br(data,log2n);
	
	unsigned int n2=1 << log2n;
	unsigned int n4=n2/2;
	if(n4 != WTablepSize) {
		Complex wp=1.0+wpTablep[log2n];
		WTablep=new(WTablep,n4) Complex;
		WTablepSize=n4;
		WTablep[0]=Complex(0.0,0.5);
		for(i=1; i < n4; i++) WTablep[i]=WTablep[i-1]*wp;
	}
	
	data[0]=Complex(1.0,1.0)*data[0].re+Complex(1.0,-1.0)*data[0].im;
#pragma ivdep	
	for(i=1; i < n4; i++) {
		Complex u=data[i], v=conj(data[n2-i]);
		Complex A=0.5*(u+v), B=WTablep[i]*(u-v);
		data[i]=A-B;
		data[n2-i]=conj(A+B);
	}
}

void rfft_brinv(Complex *data, unsigned int log2n)
{		 
	unsigned int i;
	unsigned int n2=1 << log2n;
	unsigned int n4=n2/2;
	if(n4 != WTablenSize) {
		Complex wp=1.0+wpTablen[log2n];
		WTablen=new(WTablen,n4) Complex;
		WTablenSize=n4;
		WTablen[0]=Complex(0.0,-0.5);
		for(i=1; i < n4; i++) WTablen[i]=WTablen[i-1]*wp;
	}
	
	data[0]=Complex(0.5,0.5)*data[0].re+Complex(0.5,-0.5)*data[0].im;
#pragma ivdep	
	for(i=1; i < n4; i++) {
		Complex u=data[i], v=conj(data[n2-i]);
		Complex A=0.5*(u+v), B=WTablen[i]*(u-v);
		data[i]=A-B;
		data[n2-i]=conj(A+B);
	}
	
	fft_brinv(data,log2n);
}

// Compute H = F (*) G, where F and G contain the non-negative Fourier
// components of real functions f and g, respectively. Dealiasing via
// zero-padding is implemented automatically.
//
// Arrays F[n/2], G[n/2] must be distinct, with n=2^log2n.
// Input F[i] (0 <= i < m <= n/3), g[i] (0 <= i < n/2).
// Output H[i] = F (*) G  (0 <= i < m), f[i], g[i] (0 <= i < n/2).
//
// Array H[n/2] can coincide with either F or G, in which case the output H
// subsumes f or g, respectively.

void convolve0(Complex *H, Complex *F, Complex *g, unsigned int m, unsigned
			   int log2n)
{
	unsigned int n=1 << log2n;
	unsigned int i, n2=n/2;

	for(i=m; i < n/2; i++) F[i]=0.0;
	rfft_brinv(F,log2n-1);
	
	Real ninv=2.0/n2;
#pragma ivdep	
	for(i=0; i < n2; i++) {
		H[i].re=F[i].re*g[i].re*ninv;
		H[i].im=F[i].im*g[i].im*ninv;
	}
	
	rfft_br(H,log2n-1);
	H[0].im=0.0;
}

// Compute H = F (*) G, where F and G contain the non-negative Fourier
// components of real functions f and g, respectively. Dealiasing via
// zero-padding is implemented automatically.
//
// Arrays F[n/2], G[n/2] must be distinct, with n=2^log2n.
// Input F[i], G[i] (0 <= i < m), where m <= n/3.
// Output H[i] = F (*) G  (0 <= i < m), F[i]=f[i], G[i]=g[i] (0 <= i < n/2).
//
// Array H[n/2] can coincide with either F or G, in which case the output H
// subsumes f or g, respectively.

void convolve(Complex *H, Complex *F, Complex *G, unsigned int m, unsigned
			  int log2n)
{
	unsigned int n=1 << log2n;
	unsigned int i,mmax;

//	if(m > n/3) msg(ERROR, "Insufficient room for dealiasing");
	
	if(log2n > TableSize) {
		wpTablep=new(wpTablep,log2n) Complex;
		wpTablen=new(wpTablen,log2n) Complex;
		
		mmax=1 << TableSize;
		while (n > mmax) {
			mmax <<= 1;
			Real arg=twopi/mmax;
			wpTablep[TableSize]=expim1(arg);
			wpTablen[TableSize]=expim1(-arg);
			TableSize++;
		}
	}
	
#pragma ivdep	
	for(i=m; i < n/2; i++) G[i]=0.0;
	rfft_brinv(G,log2n-1);
	
	convolve0(H,F,G,m,log2n);
}	

// Compute H = F (*) G, where F and G contain the non-negative Fourier
// components of real functions f and g, respectively, via direct convolution
// instead of a Fast Fourier Transform technique.
//
// Input F[i], G[i] (0 <= i < m).
// Output H[i] = F (*) G  (0 <= i < m), F and G unchanged.
//
// Array H[m] must be distinct from F and G.

void convolve_direct(Complex *H, Complex *F, Complex *G, unsigned int m)
{
	unsigned int i,j;
	for(i=0; i < m; i++) {
		H[i]=0.0;
		for(j=0; j <= i; j++) H[i] += F[j]*G[i-j];
		for(j=i+1; j < m; j++) H[i] += F[j]*conj(G[j-i]);
		for(j=1; j < m-i; j++) H[i] += conj(F[j])*G[i+j];
	}
}	

