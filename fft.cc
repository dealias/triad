#include "utils.h"

void fft1(Complex data[], unsigned int log2n, int isign)
{
	unsigned long mmax,m,j,istep,i,n;
	static Complex *wpTable;
	static unsigned long TableSize=0;

	n=1 << log2n;
	j=0;
	for (i=0; i < n; i++) {
		if (j > i) {Complex *p=data+i, *q=data+j;
					Real temp;
					temp=q->re; q->re=p->re; p->re=temp;
					temp=p->im; p->im=q->im; q->im=temp;
				}
		
		m=n >> 1;
		while (m >= 1 && j >= m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	
	if(log2n > TableSize) wpTable=new(wpTable,log2n) Complex;
	mmax=1 << TableSize;
	Real pmtwopi=isign*twopi;
	while (n > mmax) {
		mmax <<= 1;
		wpTable[TableSize++]=expim1(pmtwopi/mmax);
	}
	
	mmax=1; 
	Complex *wp=wpTable, *pstop=data+n;
 	while (n > mmax) {
		istep=mmax << 1;
		Real wre=1.0, wim=0.0;
		Real wpre=wp->re, wpim=wp->im;
		wp++;
		for (m=0; m < mmax; m++) {
			Complex *p,*q;
			for (p=data+m; p < pstop; p += istep) {
				q=p+mmax;
				Real tempre=wre*q->re-wim*q->im;
				Real tempim=wre*q->im+wim*q->re;
				q->re=p->re-tempre;
				p->re += tempre;
				q->im=p->im-tempim;
				p->im += tempim;
			}
			Real wtemp=wre;
			wre += wre*wpre-wim*wpim;
			wim += wim*wpre+wtemp*wpim;
		}
		mmax=istep;
	}
}

