#include "utils.h"

static Complex *wpTablep,*wpTablen;
static Complex *WTablep,*WTablen;
static unsigned int TableSize=0, WTablepSize=0, WTablenSize=0;

void fft_init(unsigned int log2n)
{
	unsigned int n=1 << log2n;
	
	wpTablep=new(wpTablep,log2n) Complex;
	wpTablen=new(wpTablen,log2n) Complex;
		
	unsigned int mmax=1 << TableSize;
	while (n > mmax) {
		mmax <<= 1;
		Real arg=twopi/mmax;
		wpTablep[TableSize]=expim1(arg);
		wpTablen[TableSize]=expim1(-arg);
		TableSize++;
	}
}
	
void fft_br(Complex *data, unsigned int log2n)
{
	if(log2n > TableSize) fft_init(log2n);
		
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
	if(log2n > TableSize) fft_init(log2n);
	
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

void rfft_br(Complex *data, unsigned int log2n)
{		 
	unsigned int i;
	
	if(log2n > TableSize) fft_init(log2n+1);
	fft_br(data,log2n);
	
	unsigned int n2=1 << log2n;
	unsigned int n4=n2/2;
	if(WTablepSize != n4) {
		Complex wp=1.0+wpTablep[log2n];
		WTablep=new(WTablep,n4) Complex;
		WTablepSize=n4;
		WTablep[0]=Complex(0.0,0.5);
		for(i=1; i < n4; i++) WTablep[i]=WTablep[i-1]*wp;
	}
	
	Real temp=data[0].re+data[0].im;
	data[0].im=data[0].re-data[0].im;
	data[0].re=temp;
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
	if(log2n > TableSize) fft_init(log2n+1);
	
	unsigned int n2=1 << log2n;
	unsigned int n4=n2/2;
	if(WTablenSize != n4) {
		Complex wp=1.0+wpTablen[log2n];
		WTablen=new(WTablen,n4) Complex;
		WTablenSize=n4;
		WTablen[0]=Complex(0.0,-0.5);
		for(i=1; i < n4; i++) WTablen[i]=WTablen[i-1]*wp;
	}
	
	Real temp=0.5*(data[0].re+data[0].im);
	data[0].im=0.5*(data[0].re-data[0].im);
	data[0].re=temp;
#pragma ivdep	
	for(i=1; i < n4; i++) {
		Complex u=data[i], v=conj(data[n2-i]);
		Complex A=0.5*(u+v), B=WTablen[i]*(u-v);
		data[i]=A-B;
		data[n2-i]=conj(A+B);
	}
	
	fft_brinv(data,log2n);
}

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
