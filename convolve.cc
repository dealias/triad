#include "options.h"
#include "fft.h"

const int bitreverse=1; // Use faster bit-reversed FFT's if available

// Compute H = F (*) G, where F and G contain the non-negative Fourier
// components of real functions f and g, respectively. Dealiasing via
// zero-padding is implemented automatically.
//
// Arrays F[n/2+1], G[n/2+1] must be distinct, with n=2^log2n.
// Input F[i] (0 <= i < m <= n/3), g[i] (0 <= i < n/2).
// For a 1D convolution, 3*m must be less than n+2.
// Output H[i] = F (*) G  (0 <= i < m), f[i], g[i] (0 <= i < n/2).
//
// Array H[n/2+1] can coincide with either F or G, in which case the output H
// subsumes f or g, respectively.

void convolve0(Complex *H, Complex *F, Complex *g, unsigned int m, unsigned
			   int log2n)
{
	const unsigned int n=1 << log2n, n2=n/2;
	unsigned int i;

#pragma ivdep	
	for(i=m; i < n2+1; i++) F[i]=0.0;
	crfft(F,log2n,-1,1.0,-bitreverse);
	
#pragma ivdep	
	for(i=0; i < n2; i++) {
		H[i].re=F[i].re*g[i].re;
		H[i].im=F[i].im*g[i].im;
	}
	
	rcfft(H,log2n,1,1.0/n,bitreverse);
}

// Compute H = F (*) G, where F and G contain the non-negative Fourier
// components of real functions f and g, respectively. Dealiasing via
// zero-padding is implemented automatically.
//
// Arrays F[n/2+1], G[n/2+1] must be distinct, with n=2^log2n.
// Input F[i], G[i] (0 <= i < m), where m <= n/3.
// For a 1D convolution, 3*m must be less than n+2.
// Output H[i] = F (*) G  (0 <= i < m), F[i]=f[i], G[i]=g[i] (0 <= i < n/2).
//
// Array H[n/2+1] can coincide with either F or G, in which case the output H
// subsumes f or g, respectively.

void convolve(Complex *H, Complex *F, Complex *G, unsigned int m, unsigned
			  int log2n)
{
	unsigned int n=1 << log2n;
	unsigned int i;

#pragma ivdep	
	for(i=m; i < n/2+1; i++) G[i]=0.0;
	crfft(G,log2n,-1,1.0,-bitreverse);
	
	convolve0(H,F,G,m,log2n);
}	

// Compute H = F (*) G, where F and G contain the non-negative Fourier
// components of real functions f and g, respectively, via direct convolution
// instead of a Fast Fourier Transform technique.
//
// Input F[i], G[i] (0 <= i < m).
// Output H[i] = F (*) G  (0 <= i < m), F and G unchanged.
//
// Array H[m] must be distinct from F[m] and G[m].

void convolve_direct(Complex *H, Complex *F, Complex *G, unsigned int m)
{
	unsigned int i,j;
	for(i=0; i < m; i++) {
		H[i]=0.0;
#pragma ivdep	
		for(j=0; j <= i; j++) H[i] += F[j]*G[i-j];
#pragma ivdep	
		for(j=i+1; j < m; j++) H[i] += F[j]*conj(G[j-i]);
#pragma ivdep	
		for(j=1; j < m-i; j++) H[i] += conj(F[j])*G[i+j];
	}
}	

