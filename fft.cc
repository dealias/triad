#include "utils.h"

extern Complex *wpTablep,*wpTablen;

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
				Real tempim=p->im-q->im;
				p->re += q->re;
			    q->re=c*tempre-s*tempim;
				p->im += q->im;
				q->im=c*tempim+s*tempre;
			}
			Real wtemp=c;
			c += c*wpre-s*wpim;
			s += s*wpre+wtemp*wpim;
		}
		istep=mmax;
	}
}

