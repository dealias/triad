#include "utils.h"

static Complex *wpTablep;
static Complex *wpTablen;
static unsigned int TableSize=0;

void fft_br(Complex *data, unsigned int log2n)
{
	unsigned int n=1 << log2n;
	unsigned int mmax=1;
	Complex *pstop=data+n,*wp=wpTablep;
	
 	while (mmax < n) {
		unsigned int istep=mmax << 1;
		Real c=1.0, s=0.0;
		Real wpre=wp->re, wpim=wp->im;
		wp++;
		for (unsigned int m=0; m < mmax; m++) {
			Complex *p,*q;
#pragma ivdep			
			for (p=data+m; p < pstop; p += istep) {
				q=p+mmax;
				Real tempre=c*q->re-s*q->im;
				Real tempim=c*q->im+s*q->re;
				q->re=p->re-tempre;
				p->re += tempre;
				q->im=p->im-tempim;
				p->im += tempim;
			}
			Real wtemp=c;
			c += c*wpre-s*wpim;
			s += s*wpre+wtemp*wpim;
		}
		mmax=istep;
	}
}

void fft_brinv(Complex *data, unsigned int log2n)
{
	unsigned int n=1 << log2n;
	unsigned int istep=n;
	Complex *pstop=data+n,*wp=wpTablen+log2n-1;

 	while (istep > 1) {
		unsigned int mmax=istep >> 1;
		Real c=1.0, s=0.0;
		Real wpre=wp->re, wpim=wp->im;
		wp--;
		for (unsigned int m=0; m < mmax; m++) {
			Complex *p,*q;
#pragma ivdep			
			for (p=data+m; p < pstop; p += istep) {
				q=p+mmax;
				Real tempre=p->re-q->re;
				p->re += q->re;
				Real tempim=p->im-q->im;
				p->im += q->im;
			    q->re=c*tempre-s*tempim;
				q->im=c*tempim+s*tempre;
			}
			Real wtemp=c;
			c += c*wpre-s*wpim;
			s += s*wpre+wtemp*wpim;
		}
		istep=mmax;
	}
}

inline void rfft_br(Complex *data, unsigned int log2n)
{		 
	fft_br(data,log2n);
	
	unsigned int n2=1 << log2n;
	unsigned int n4=n2/2;
	Complex wp=1.0+wpTablep[log2n], W=Complex(0.0,0.5)*wp;
	
	data[0]=Complex(1.0,1.0)*data[0].re+Complex(1.0,-1.0)*data[0].im;
	for(unsigned int i=1; i < n4; i++) {
		Complex u=data[i], v=conj(data[n2-i]);
		Complex A=0.5*(u+v), B=W*(u-v);
		data[i]=A-B;
		data[n2-i]=conj(A+B);
		W *= wp;
	}
}

inline void rfft_brinv(Complex *data, unsigned int log2n)
{		 
	unsigned int n2=1 << log2n;
	unsigned int n4=n2/2;
	Complex wp=1.0+wpTablen[log2n], W=Complex(0.0,0.5)*wp;
	
	data[0]=Complex(0.5,0.5)*data[0].re+Complex(0.5,-0.5)*data[0].im;
	for(unsigned int i=1; i < n4; i++) {
		Complex u=data[i], v=conj(data[n2-i]);
		Complex A=0.5*(u+v), B=W*(u-v);
		data[i]=A+B;
		data[n2-i]=conj(A-B);
		W *= wp;
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
// Array H[m] can coincide with either F or G, in which case the output H
// subsumes f or g, respectively.

void convolve0(Complex *H, Complex *F, Complex *g, unsigned int m, unsigned
			   int log2n)
{
	unsigned int n=1 << log2n;
	unsigned int i, n2=n/2;

	for(i=m; i < n/2; i++) F[i] = 0.0;
	rfft_brinv(F,log2n-1);
	
	Real ninv=2.0/n2;
	for(i=0; i < n2; i++) {
		H[i].re=F[i].re*g[i].re*ninv;
		H[i].im=F[i].im*g[i].im*ninv;
	}
	
	rfft_br(H,log2n-1);
}

// Compute H = F (*) G, where F and G contain the non-negative Fourier
// components of real functions f and g, respectively. Dealiasing via
// zero-padding is implemented automatically.
//
// Arrays F[n/2], G[n/2] must be distinct, with n=2^log2n.
// Input F[i], G[i] (0 <= i < m), where m <= n/3.
// Output H[i] = F (*) G  (0 <= i < m), F[i]=f[i], G[i]=g[i] (0 <= i < n/2).
//
// Array H[m] can coincide with either F or G, in which case the output H
// subsumes f or g, respectively.

void convolve(Complex *H, Complex *F, Complex *G, unsigned int m, unsigned
			  int log2n)
{
	unsigned int n=1 << log2n;
	unsigned int i,mmax;

	if(m > n/3) msg(ERROR, "Insufficient room for dealiasing");
	
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
	
	for(i=m; i < n/2; i++) G[i] = 0.0;
	rfft_brinv(G,log2n-1);
	
	convolve0(H,F,G,m,log2n);
}	


