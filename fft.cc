#include "options.h"
#include "fft.h"

extern Complex *wpTable;
extern unsigned int wpTableSize;

// Return the bit-reversed forward Fourier transform of a Complex vector.
// Before calling, data must be allocated as Complex[n].
// On entry: data contains the n Complex values.
//           log2n contains the base-2 logarithm of n.
// On exit:  data contains the n Complex Fourier values.

void fft_br(Complex *data, unsigned int log2n)
{
	if(log2n > wpTableSize) fft_init(log2n);
	
	unsigned int n=1 << log2n;
	Complex *pstop=data+n,*wp=wpTable+2;
	
	if(n > 1) {
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

// Return the bit-reversed inverse Fourier transform of a Complex vector.
// Before calling, data must be allocated as Complex[n].
// On entry: data contains the n Complex values.
//           log2n contains the base-2 logarithm of n.
// On exit:  data contains the n Complex bit-reversed inverse Fourier values.
// Note: When computing an inverse transform, the result must be divided by n.

void fft_brinv(Complex *data, unsigned int log2n)
{
	if(log2n > wpTableSize) fft_init(log2n);
	
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

// Return the Fourier transform of a Complex vector.
// Before calling, data must be allocated as Complex[n].
// On entry: data contains the n Complex values.
//           log2n contains the base-2 logarithm of n.
//           isign is +1 for a forward transform, -1 for an inverse transform.
//           bitreverse is 0 for a true fft, 1 for a faster bit-reversed fft.
// On exit:  data contains the n Complex Fourier values.
// Note: When computing an inverse transform, the result must be divided by n.

void fft(Complex *data, unsigned int log2n, int isign, int bitreverse)
{
	if(bitreverse) {
		if(isign == 1) fft_br(data,log2n);
		else fft_brinv(data,log2n);
		return;
	}
	
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
	
	if(log2n > wpTableSize) fft_init(log2n);
	
	Complex *pstop=data+n, *wp=wpTable+1;
	
	if(n > 1) {
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
	
	unsigned int mmax=2; 
 	while (mmax < n) {
		istep=mmax << 1;
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
		Real wpre=wp->re, wpim=wp->im*isign;
		wp++;
		Real c=1.0+wpre, s=wpim;
		for(m=1; m < mmax-1; m++) {
#pragma ivdep			
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
	}
}

// Return the Fourier transform of nk Complex vector's.
// Before calling, data must be allocated as Complex[nk*n].
// On entry: data[k+nk*j] contains the n Complex values for each k=0,...,nk.
//           log2n contains the base-2 logarithm of n.
//           isign is +1 for a forward transform, -1 for an inverse transform.
// On exit:  data contains the n Complex Fourier values for each k=0,...nk.
// Note: When computing an inverse transform, the result must be divided by n.

void mfft(Complex *data, unsigned int log2n, int isign, unsigned int nk)
{
	unsigned int i,j,k,n,m,istep;
	static unsigned int tempsize=0;
	static Complex *temp;

	n=1 << log2n;

	if(nk > tempsize) {
		tempsize=nk;
		temp=new(temp,tempsize) Complex;
	}
	
	j=0;
	for(i=0; i < n-1; i++) {
		if (j > i) {
			Complex *p=data+i*nk, *q=data+j*nk;
#pragma ivdep			
			for(k=0; k < nk; k++) {
				temp[k].re=q[k].re; q[k].re=p[k].re; p[k].re=temp[k].re;
				temp[k].im=p[k].im; p[k].im=q[k].im; q[k].im=temp[k].im;
			}
		}
		
		m=n >> 1;
		while (j >= m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	
	if(log2n > wpTableSize) fft_init(log2n);
	
	Complex *pstop=data+n*nk, *wp=wpTable+1;
	
	if(n > 1) {
		for(Complex *p=data; p < pstop; p += 2*nk) {
			Complex *q=p+nk;
#pragma ivdep			
			for(k=0; k < nk; k++) {
				temp[k].re=p[k].re-q[k].re;
				temp[k].im=p[k].im-q[k].im;
				p[k].re += q[k].re;
				p[k].im += q[k].im;
				q[k].re=temp[k].re;
				q[k].im=temp[k].im;
			}
		}
	}
	
	unsigned int mmax=2; 
 	while (mmax < n) {
		istep=mmax << 1;
		int istepnk=istep*nk;
		int mmaxnk=mmax*nk;
		Complex *p;
		for(p=data; p < pstop; p += istepnk) {
			Complex *q=p+mmaxnk;
#pragma ivdep			
			for(k=0; k < nk; k++) {
				temp[k].re=p[k].re-q[k].re;
				temp[k].im=p[k].im-q[k].im;
				p[k].re += q[k].re;
				p[k].im += q[k].im;
				q[k].re=temp[k].re;
				q[k].im=temp[k].im;
			}
		}
		Real wpre=wp->re, wpim=wp->im*isign;
		wp++;
		Real c=1.0+wpre, s=wpim;
		for(m=1; m < mmax-1; m++) {
			for(p=data+m*nk; p < pstop; p += istepnk) {
				Complex *q=p+mmaxnk;
#pragma ivdep			
				for(k=0; k < nk; k++) {
					temp[k].re=c*q[k].re-s*q[k].im;
					temp[k].im=c*q[k].im+s*q[k].re;
					q[k].re=p[k].re-temp[k].re;
					q[k].im=p[k].im-temp[k].im;
					p[k].re += temp[k].re;
					p[k].im += temp[k].im;
				}
			}
			Real wtemp=c;
			c += c*wpre-s*wpim;
			s += s*wpre+wtemp*wpim;
		}
		for(p=data+mmaxnk-nk; p < pstop; p += istepnk) {
			Complex *q=p+mmaxnk;
#pragma ivdep			
			for(k=0; k < nk; k++) {
				temp[k].re=c*q[k].re-s*q[k].im;
				temp[k].im=c*q[k].im+s*q[k].re;
				q[k].re=p[k].re-temp[k].re;
				q[k].im=p[k].im-temp[k].im;
				p[k].re += temp[k].re;
				p[k].im += temp[k].im;
			}
		}
		mmax=istep;
	}
}

// Return the Fourier transform of a Complex vector of length n=power of 4.
// Before calling, data must be allocated as Complex[n].
// On entry: data contains the n Complex values.
//           log4n contains the base-4 logarithm of n.
//           isign is +1 for a forward transform, -1 for an inverse transform.
// On exit:  data contains the n Complex Fourier values.
// Note: When computing an inverse transform, the result must be divided by n.

void fft4(Complex *data, unsigned int log4n, int isign)
{
	unsigned int m,k,j;
	static unsigned int lastsize=0;
	static Complex *phase;

	m=1 << log4n;
	
	if(m != lastsize) {
		if(m > lastsize) phase=new(phase,m*m) Complex;
		lastsize=m;
		Complex factor=exp(twopi*I/(m*m));
		Complex kphase=1.0;
		for(k=0; k < m; k++) {
			Complex *p=phase+m*k;
			p[0]=1.0;
			for(j=1; j < m; j++) p[j]=p[j-1]*kphase;
			kphase *= factor;
		}
	}
	
	mfft(data,log4n,isign,m);
	
	if(isign == 1) {
		for(j=0; j < m; j++) {
			int mj=m*j;
			Complex *p=data+mj;
			Complex *c=phase+mj;
			Complex *q=data+j;
			p[j] *= c[j];
			for(k=0; k < j; k++, q += m) {
				Complex temp = p[k];
				p[k] = (*q)*c[k];
				(*q)=temp*c[k];
			}
		}
	} else {
		for(j=0; j < m; j++) {
			int mj=m*j;
			Complex *p=data+mj;
			Complex *c=phase+mj;
			Complex *q=data+j;
			p[j] *= conj(c[j]);
			for(k=0; k < j; k++, q += m) {
				Complex temp = p[k];
				p[k] = (*q)*conj(c[k]);
				(*q)=temp*conj(c[k]);
			}
		}
	}
	
	mfft(data,log4n,isign,m);
}

