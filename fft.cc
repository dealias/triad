#include "utils.h"

void fft(Complex *data, unsigned int log2n, int isign)
{
	unsigned int mmax,m,j,istep,i,n;
	static Complex *wpTable[2];
	static unsigned int TableSize[2]={0,0};

	n=1 << log2n;

	j=0;
	for (i=0; i < n-1; i++) {
		if (j > i) {Complex *p=data+i, *q=data+j;
					Real temp;
					temp=q->re; q->re=p->re; p->re=temp;
					temp=p->im; p->im=q->im; q->im=temp;
				}
		
		m=n >> 1;
		while (j >= m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	
	int index=(isign == 1) ? 1 : 0;
	unsigned int TableSizei=TableSize[index];
	if(log2n > TableSizei) wpTable[index]=new(wpTable[index],log2n) Complex;
	mmax=1 << TableSizei;
	Real pmtwopi=isign*twopi;
	while (mmax < n) {
		mmax <<= 1;
		wpTable[index][TableSizei++]=expim1(pmtwopi/mmax);
	}
	
	mmax=1; 
	Complex *wp=wpTable[index], *pstop=data+n;
 	while (mmax < n) {
		istep=mmax << 1;
		Real c=1.0, s=0.0;
		Real wpre=wp->re, wpim=wp->im;
		wp++;
		for (m=0; m < mmax; m++) {
			Complex *p,*q;
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

