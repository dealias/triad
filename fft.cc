#include "utils.h"

static Complex *wpTablep,*wpTablen;
static Complex *WTablep,*WTablen;
static unsigned int TableSize=0, WTablepSize=0, WTablenSize=0;
void fft(Complex *data, unsigned int log2n, int isign);

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
	Complex *pstop=data+n,*wp=wpTablep+2;
	
 	if (n > 1) {
#pragma ivdep			
		for(Complex *p=data; p < pstop; p += 2) {
			Complex *p1=p+1;
			Real tempre=p->re-p1->re;
			Real tempim=p->im-p1->im;
			p->re += p1->re;
			p->im += p1->im;
			p1->re=tempre;
			p1->im=tempim;
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
	Complex *pstop=data+n,*wp=wpTablen+log2n-1;

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
		
		Real c=1.0+wp->re, s=wp->im;
		for(unsigned int m=1; m < mmax-1; m += 2) {
			Real c2=c+c*wp->re-s*wp->im;
			Real s2=s+s*wp->re+c*wp->im;
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
			c=c2+c2*wp->re-s2*wp->im;
			s=s2+s2*wp->re+c2*wp->im;
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
			Complex *p1=p+1;
			Real tempre=p->re-p1->re;
			Real tempim=p->im-p1->im;
			p->re += p1->re;
			p->im += p1->im;
			p1->re=tempre;
			p1->im=tempim;
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
	unsigned int n4=n2/2;
	if(WTablepSize != n4) {
		Complex wp=1.0+wpTablep[log2n];
		WTablep=new(WTablep,n4) Complex;
		WTablepSize=n4;
		WTablep[0]=Complex(0.0,0.5);
		for(i=1; i < n4; i++) WTablep[i]=WTablep[i-1]*wp;
	}
	
	data[n2]=data[0].re-data[0].im;
	data[0]=data[0].re+data[0].im;
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
	
	data[0].im=0.5*(data[0].re-data[n2].re);
	data[0].re=0.5*(data[0].re+data[n2].re);
#pragma ivdep	
	for(i=1; i < n4; i++) {
		Complex u=data[i], v=conj(data[n2-i]);
		Complex A=0.5*(u+v), B=WTablen[i]*(u-v);
		data[i]=A-B;
		data[n2-i]=conj(A+B);
	}
	
	fft_brinv(data,log2n);
//	fft(data,log2n,-1);
}

void fft(Complex *data, unsigned int log2n, int isign)
{
	unsigned int mmax,m,j,istep,i,n;
	static Complex *wpTable[2];
	static unsigned int TableSize[2]={0,0};

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
		for(m=0; m < mmax; m++) {
			Complex *p,*q;
			for(p=data+m; p < pstop; p += istep) {
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
