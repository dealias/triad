#include "utils.h"

static Complex *wpTable;
static Complex *WTable;
static unsigned int TableSize=0, WTableSize=0;

void fft(Complex *data, unsigned int log2n, int isign);

void fft_init(unsigned int log2n)
{
	unsigned int n=1 << log2n;
	
	wpTable=new(wpTable,log2n) Complex;
	unsigned int mmax=1 << TableSize;
	while (n > mmax) {
		mmax <<= 1;
		wpTable[TableSize]=expim1(twopi/mmax);
		TableSize++;
	}
}
	
void rfft_init(unsigned int log2n)
{
	unsigned int n4=1 << (log2n-1);
	
	WTable=new(WTable,n4) Complex;
	WTableSize=n4;
	WTable[0]=Complex(0.0,0.5);
	
	Complex wp=1.0+wpTable[log2n];
	for(unsigned int i=1; i < n4; i++) WTable[i]=WTable[i-1]*wp;
}

void fft_br(Complex *data, unsigned int log2n)
{
	if(log2n > TableSize) fft_init(log2n);
		
	unsigned int n=1 << log2n;
	Complex *pstop=data+n,*wp=wpTable+2;
	
 	if (n > 1) {
#pragma ivdep			
		for(Complex *p=data; p < pstop; p += 2) {
			Complex *q=p+1;
			Real tempre=p->re-q->re;
			Real tempim=p->im-q->im;
			p->re += q->re;
			p->im += q->im;
			q->re=tempre;
			q->im=tempim;
		}
	}
	
 	if(n > 2) {
#pragma ivdep			
		for(Complex *p=data; p < pstop; p += 4) {
			Complex *p1=p+1,*p2=p1+1,*p3=p2+1;
			Real tempre=p->re;
			Real tempim=p->im;
			Real temp2re=p1->re;
			Real temp2im=p1->im-p3->re;
			p->re += p2->re;
			p->im += p2->im;
			p1->re -= p3->im;
			p1->im += p3->re;
			p2->re=tempre-p2->re;
			p2->im=tempim-p2->im;
			p3->re=temp2re+p3->im;
			p3->im=temp2im;
		}
	}
	
	unsigned int mmax=4;
 	while (mmax < n) {
		unsigned int istep=mmax << 1;
		Complex *p;
#pragma ivdep			
		for(p=data; p < pstop; p += istep) {
			Complex *q=p+mmax;
			Real tempre=p->re-q->re;
			Real tempim=p->im-q->im;
			p->re += q->re;
			p->im += q->im;
			q->re=tempre;
			q->im=tempim;
		}

		Real c=1.0+wp->re, s=wp->im;
		for(unsigned int m=1; m < mmax-1; m += 2) {
			Real c2=c+c*wp->re-s*wp->im;
			Real s2=s+s*wp->re+c*wp->im;
#pragma ivdep			
			for(p=data+m; p < pstop; p += istep) {
				Complex *p1=p+1,*q=p+mmax,*q1=q+1;
				Real tempre=c*q->re-s*q->im;
				Real tempim=c*q->im+s*q->re;
				Real temp2re=c2*q1->re-s2*q1->im;
				Real temp2im=c2*q1->im+s2*q1->re;
				q->re=p->re-tempre;
				q->im=p->im-tempim;
				q1->re=p1->re-temp2re;
				q1->im=p1->im-temp2im;
				p->re += tempre;
				p->im += tempim;
				p1->re += temp2re;
				p1->im += temp2im;
			}
			c=c2+c2*wp->re-s2*wp->im;
			s=s2+s2*wp->re+c2*wp->im;
		}
#pragma ivdep			
			for(p=data+mmax-1; p < pstop; p += istep) {
				Complex *q=p+mmax;
				Real tempre=c*q->re-s*q->im;
				Real tempim=c*q->im+s*q->re;
				q->re=p->re-tempre;
				q->im=p->im-tempim;
				p->re += tempre;
				p->im += tempim;
			}
		mmax=istep;
		wp++;
	}
}

void fft_brinv(Complex *data, unsigned int log2n)
{
	if(log2n > TableSize) fft_init(log2n);
	
	unsigned int n=1 << log2n;
	unsigned int istep=n;
	Complex *pstop=data+n,*wp=wpTable+log2n-1;

 	while (istep > 4) {
		unsigned int mmax=istep >> 1;
		Complex *p;
#pragma ivdep			
		for(p=data; p < pstop; p += istep) {
			Complex *q=p+mmax;
			Real tempre=p->re-q->re;
			Real tempim=p->im-q->im;
			p->re += q->re;
			p->im += q->im;
			q->re=tempre;
			q->im=tempim;
		}
		
		Real c=1.0+wp->re, s=-wp->im;
		for(unsigned int m=1; m < mmax-1; m += 2) {
			Real c2=c+c*wp->re+s*wp->im;
			Real s2=s+s*wp->re-c*wp->im;
#pragma ivdep			
			for(p=data+m; p < pstop; p += istep) {
				Complex *p1=p+1,*q=p+mmax,*q1=q+1;
				Real tempre=p->re-q->re;
				Real tempim=p->im-q->im;
				Real temp2re=p1->re-q1->re;
				Real temp2im=p1->im-q1->im;
				p->re += q->re;
				p->im += q->im;
 				p1->re += q1->re;
				p1->im += q1->im;
			    q->re=c*tempre-s*tempim;
				q->im=c*tempim+s*tempre;
			    q1->re=c2*temp2re-s2*temp2im;
				q1->im=c2*temp2im+s2*temp2re;
			}
			c=c2+c2*wp->re+s2*wp->im;
			s=s2+s2*wp->re-c2*wp->im;
		}
#pragma ivdep			
		for(p=data+mmax-1; p < pstop; p += istep) {
			Complex *q=p+mmax;
			Real tempre=p->re-q->re;
			Real tempim=p->im-q->im;
			p->re += q->re;
			p->im += q->im;
			q->re=c*tempre-s*tempim;
			q->im=c*tempim+s*tempre;
		}
		istep=mmax;
		wp--;
	}
	
 	if (n > 2) {
#pragma ivdep			
		for(Complex *p=data; p < pstop; p += 4) {
			Complex *p1=p+1,*p2=p1+1,*p3=p2+1;
			Real tempre=p->re;
			Real tempim=p->im;
			Real temp2re=p3->re-p1->re;
			Real temp2im=p1->im;
			p->re += p2->re;
			p->im += p2->im;
			p1->re += p3->re;
			p1->im += p3->im;
			p2->re=tempre-p2->re;
			p2->im=tempim-p2->im;
			p3->re=temp2im-p3->im;
			p3->im=temp2re;
		}
	}
	
	if (n > 1) {
#pragma ivdep			
		for(Complex *p=data; p < pstop; p += 2) {
			Complex *q=p+1;
			Real tempre=p->re-q->re;
			Real tempim=p->im-q->im;
			p->re += q->re;
			p->im += q->im;
			q->re=tempre;
			q->im=tempim;
		}
	}
}

void rfft_br(Complex *data, unsigned int log2n)
{		 
	unsigned int i;
	
	if(log2n > TableSize) fft_init(log2n+1);
	fft_br(data,log2n);
//	fft(data,log2n,1);
	
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

void rfft_brinv(Complex *data, unsigned int log2n)
{		 
	unsigned int i;
	if(log2n > TableSize) fft_init(log2n+1);
	
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
	
	fft_brinv(data,log2n);
//	fft(data,log2n,-1);
}

void fft(Complex *data, unsigned int log2n, int isign)
{
	unsigned int m,j,istep,i,n;

	n=1 << log2n;

	j=0;
	for(i=0; i < n-1; i++) {
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
	
	if(log2n > TableSize) fft_init(log2n);
	
	Complex *pstop=data+n, *wp=wpTable+1;
	
	if(n > 1) {
		for(Complex *p=data; p < pstop; p += 2) {
			Complex *q=p+1;
			Real tempre=p->re-q->re;
			Real tempim=p->im-q->im;
			p->re += q->re;
			p->im += q->im;
			q->re=tempre;
			q->im=tempim;
		}
	}
	
	unsigned int mmax=2; 
 	while (mmax < n) {
		istep=mmax << 1;
		Complex *p;
		for(p=data; p < pstop; p += istep) {
			Complex *q=p+mmax;
			Real tempre=p->re-q->re;
			Real tempim=p->im-q->im;
			p->re += q->re;
			p->im += q->im;
			q->re=tempre;
			q->im=tempim;
		}
		Real wpre=wp->re, wpim=wp->im*isign;
		wp++;
		Real c=1.0+wpre, s=wpim;
		for(m=1; m < mmax-1; m++) {
			for(p=data+m; p < pstop; p += istep) {
				Complex *q=p+mmax;
				Real tempre=c*q->re-s*q->im;
				Real tempim=c*q->im+s*q->re;
				q->re=p->re-tempre;
				q->im=p->im-tempim;
				p->re += tempre;
				p->im += tempim;
			}
			Real wtemp=c;
			c += c*wpre-s*wpim;
			s += s*wpre+wtemp*wpim;
		}
		for(p=data+mmax-1; p < pstop; p += istep) {
			Complex *q=p+mmax;
			Real tempre=c*q->re-s*q->im;
			Real tempim=c*q->im+s*q->re;
			q->re=p->re-tempre;
			q->im=p->im-tempim;
			p->re += tempre;
			p->im += tempim;
		}
		mmax=istep;
	}
}
