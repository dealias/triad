/* Efficient Base Fast Fourier Transform routines written in C++.
   Version 1.0 Copyright (C) 1997 John C. Bowman (bowman@ipp-garching.mpg.de)
   
   This file provides two efficient portable fft base routines:
   
   fft:    Single One-Dimensional Fast Fourier Transform
   mfft:   Multiple One-Dimensional Fast Fourier Transform
   
   To obtain slightly improved performance, these routines may be overlayed
   by the following interfaces to existing architecture-specific libraries:
   
   efft.cc (libessl)
   cfft.cc (libsci on Cray MVP systems like the J90)
   tfft.cc (libsci on Cray MPP systems like the T3E)

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. */

#include "options.h"
#include "fft.h"

static void fft_brinv(Complex *data, unsigned int log2n, int isign, 
					  Real scale=1.0);
static void mfft_brinv(Complex *data, unsigned int log2n, int isign,
					   unsigned int nk, unsigned int inc1=0,
					   unsigned int inc2=1, Real scale=1.0);
	
extern Complex *wpTable;
extern unsigned int wpTableSize;

// Return the Fourier transform of a Complex vector.
// Before calling, data must be allocated as Complex[n].
// On entry: data contains the n Complex values;
//           log2n contains the base-2 logarithm of n;
//           isign is the sign (+/- 1) of the phase;
//           [scale is a constant by which the results will be multiplied;]
//           [bitreverse is 0 for a true fft of data (default);
//                         +1 for a fft of bit-reversed data (faster);
//                         -1 for a bit-reversed fft of data (faster).]
// On exit:  data contains the n Complex Fourier values.
// Note: To compute a true inverse transform, set scale=1.0/n.

void fft(Complex *data, unsigned int log2n, int isign, Real scale,
		 int bitreverse) 
{
	if(bitreverse == -1) {fft_brinv(data,log2n,isign,scale); return;}
	
	unsigned int n=1 << log2n;

	Complex *pstop=data+n;
	
	if(scale == 1.0) {
		if(!bitreverse) {
			unsigned int j=0;
			for(unsigned int i=0; i < n-1; i++) {
				if(j > i) {
					Complex *p=data+i, *q=data+j;
					Real tempre=p->re;
					Real tempim=p->im; 
					p->re=q->re;
					p->im=q->im;
					q->re=tempre;
					q->im=tempim;
				}
				
				unsigned int m=n >> 1;
				while(j >= m) {
					j -= m;
					m >>= 1;
				}
				j += m;
			}
		}
	} else {
		if(bitreverse) {
#pragma ivdep			
			for(Complex *p=data; p < pstop; p ++) {
				p->re *= scale;
				p->im *= scale;
			}
		} else {
			unsigned int j=0;
			for(unsigned int i=0; i < n-1; i++) {
				if(j >= i) {
					Complex *p=data+i;
					if(j > i) {
						Real tempre=p->re;
						Real tempim=p->im;
						Complex *q=data+j;
						p->re=q->re*scale;
						p->im=q->im*scale;
						q->re=tempre*scale;
						q->im=tempim*scale;
					} else {
						p->re *= scale;
						p->im *= scale;
					}						
				}

				unsigned int m=n >> 1;
				while(j >= m) {
					j -= m;
					m >>= 1;
				}
				j += m;
			}
			Complex *p=data+n-1;
			p->re *= scale;
			p->im *= scale;
		}
	}
	
	if(log2n > wpTableSize) fft_init(log2n);
	
	Complex *wp=wpTable+2;
	
	if(n > 1) {
#pragma ivdep			
		for(Complex *p=data; p < pstop; p += 2) {
			Complex *q=p+1;
			Real tempre=p->re;
			Real tempim=p->im;
			p->re += q->re;
			p->im += q->im;
			q->re=tempre-q->re;
			q->im=tempim-q->im;
		}
	}
	
	if(n > 2) {
		if(isign == 1) {
#pragma ivdep			
			for(Complex *p=data; p < pstop; p += 4) {
				Complex *p1=p+1,*q=p+2,*q1=p+3;
				Real tempre=p->re;
				Real tempim=p->im;
				Real temp2re=p1->re;
				Real temp2im=p1->im-q1->re;
				p->re += q->re;
				p->im += q->im;
				p1->re -= q1->im;
				p1->im += q1->re;
				q->re=tempre-q->re;
				q->im=tempim-q->im;
				q1->re=temp2re+q1->im;
				q1->im=temp2im;
			}
		} else {
#pragma ivdep			
			for(Complex *p=data; p < pstop; p += 4) {
				Complex *p1=p+1,*q=p+2,*q1=p+3;
				Real tempre=p->re;
				Real tempim=p->im;
				Real temp2re=p1->re;
				Real temp2im=p1->im+q1->re;
				p->re += q->re;
				p->im += q->im;
				p1->re += q1->im;
				p1->im -= q1->re;
				q->re=tempre-q->re;
				q->im=tempim-q->im;
				q1->re=temp2re-q1->im;
				q1->im=temp2im;
			}
		}
	}
	
	unsigned int mmax=4;
 	while(mmax < n) {
		unsigned int istep=mmax << 1;
		Complex *p;
#pragma ivdep			
		for(p=data; p < pstop; p += istep) {
			Complex *q=p+mmax;
			Real tempre=p->re;
			Real tempim=p->im;
			p->re += q->re;
			p->im += q->im;
			q->re=tempre-q->re;
			q->im=tempim-q->im;
		}

		Real wpre=wp->re;
		Real wpim=wp->im*isign;
		Real c=1.0+wpre, s=wpim;
		for(unsigned int m=1; m < mmax-1; m += 2) {
			Real c2=c+c*wpre-s*wpim;
			Real s2=s+s*wpre+c*wpim;
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
			c=c2+c2*wpre-s2*wpim;
			s=s2+s2*wpre+c2*wpim;
		}
#pragma ivdep			
		for(p=data+mmax-1; p < pstop; p += istep) {
			Complex *q=p+mmax;
			Real tempre=c*q->re-s*q->im;
			Real tempim=c*q->im+s*q->re;
			Real qre=p->re-tempre;
			Real qim=p->im-tempim;
			p->re += tempre;
			p->im += tempim;
			q->re=qre;
			q->im=qim;
		}
		mmax=istep;
		wp++;
	}
}

// Return the Fourier transform of nk Complex vectors.
// Before calling, data must be allocated as Complex[nk*n].
// On entry: data contains the n Complex values for each k=0,...,nk-1;
//           log2n contains the base-2 logarithm of n;
//           isign is the sign (+/- 1) of the phase;
//           [inc1 is the stride between the elements of each Complex vector;]
//           [inc2 is the stride between the first elements of the vectors;]
//           [scale is a constant by which the results will be multiplied;]
//           [bitreverse is 0 for a true fft of data (default);
//                         +1 for a fft of bit-reversed data (faster);
//                         -1 for a bit-reversed fft of data (faster).]
// On exit:  data contains the n Complex Fourier values for each k=0,...,nk-1.
// Note: To compute a true inverse transform, set scale=1.0/n.

void mfft(Complex *data, unsigned int log2n, int isign, unsigned int nk,
		  unsigned int inc1, unsigned int inc2, Real scale, int bitreverse)
{
#if !_CRAYMVP
	if(inc1 == 1) {
		Complex *pstop=data+nk*inc2;
		for(Complex *p=data; p < pstop; p += inc2)
			fft(p,log2n,isign,scale,bitreverse);
		return;
	}
#endif	
	
	if(bitreverse == -1) {
		mfft_brinv(data,log2n,isign,nk,inc1,inc2,scale);
		return;
	}

	unsigned int kstop=nk*inc2;
	
	unsigned int n=1 << log2n;
	Complex *pstop=data+n*inc1; 

	if(scale == 1.0) {
		if(!bitreverse) {
			unsigned int j=0;
			for(unsigned int i=0; i < n-1; i++) {
				if(j > i) {
					Complex *p=data+i*inc1, *q=data+j*inc1;
#pragma ivdep			
					for(unsigned int k=0; k < kstop; k += inc2) {
						Real tempre=p[k].re;
						Real tempim=p[k].im;
						p[k].re=q[k].re;
						p[k].im=q[k].im;
						q[k].re=tempre;
						q[k].im=tempim;
					}
				}
				
				unsigned int m=n >> 1;
				while(j >= m) {
					j -= m;
					m >>= 1;
				}
				j += m;
			}
		}
	} else {
			if(bitreverse) {
				for(Complex *p=data; p < pstop; p += inc1) {
#pragma ivdep					
					for(unsigned int k=0; k < kstop; k += inc2) {
						p[k].re *= scale;
						p[k].im *= scale;
					}
				}
			} else {
			unsigned int j=0;
			for(unsigned int i=0; i < n-1; i++) {
				if(j >= i) {
					Complex *p=data+i*inc1;
					if(j > i) {
						Complex *q=data+j*inc1;
#pragma ivdep			
						for(unsigned int k=0; k < kstop; k += inc2) {
							Real tempre=p[k].re;
							Real tempim=p[k].im;
							p[k].re=q[k].re*scale;
							p[k].im=q[k].im*scale;
							q[k].re=tempre*scale;
							q[k].im=tempim*scale;
						}
					} else {
#pragma ivdep			
						for(unsigned int k=0; k < kstop; k += inc2) {
							p[k].re *= scale;
							p[k].im *= scale;
						}
					}
				}
				
				unsigned int m=n >> 1;
				while(j >= m) {
					j -= m;
					m >>= 1;
				}
				j += m;
			}
			Complex *p=data+(n-1)*inc1;
#pragma ivdep			
			for(unsigned int k=0; k < kstop; k += inc2) {
				p[k].re *= scale;
				p[k].im *= scale;
			}
		}
	}
	
	if(log2n > wpTableSize) fft_init(log2n);
	
	Complex *wp=wpTable+2;

	if(n > 1) {
		int twoinc1=2*inc1;
		for(Complex *p=data; p < pstop; p += twoinc1) {
			Complex *q=p+inc1;
#pragma ivdep			
			for(unsigned int k=0; k < kstop; k += inc2) {
				Real tempre=p[k].re;
				Real tempim=p[k].im;
				p[k].re += q[k].re;
				p[k].im += q[k].im;
				q[k].re=tempre-q[k].re;
				q[k].im=tempim-q[k].im;
			}
		}
	}
	
 	if(n > 2) {
		int fourinc1=4*inc1;
		if(isign == 1) {
			for(Complex *p=data; p < pstop; p += fourinc1) {
				Complex *p1=p+inc1,*q=p1+inc1,*q1=q+inc1;
#pragma ivdep			
				for(unsigned int k=0; k < kstop; k += inc2) {
					Real tempre=p[k].re;
					Real tempim=p[k].im;
					Real temp2re=p1[k].re;
					Real temp2im=p1[k].im-q1[k].re;
					p[k].re += q[k].re;
					p[k].im += q[k].im;
					p1[k].re -= q1[k].im;
					p1[k].im += q1[k].re;
					q[k].re=tempre-q[k].re;
					q[k].im=tempim-q[k].im;
					q1[k].re=temp2re+q1[k].im;
					q1[k].im=temp2im;
				}
			}
		} else {
			for(Complex *p=data; p < pstop; p += fourinc1) {
				Complex *p1=p+inc1,*q=p1+inc1,*q1=q+inc1;
#pragma ivdep			
				for(unsigned int k=0; k < kstop; k += inc2) {
					Real tempre=p[k].re;
					Real tempim=p[k].im;
					Real temp2re=p1[k].re;
					Real temp2im=p1[k].im+q1[k].re;
					p[k].re += q[k].re;
					p[k].im += q[k].im;
					p1[k].re += q1[k].im;
					p1[k].im -= q1[k].re;
					q[k].re=tempre-q[k].re;
					q[k].im=tempim-q[k].im;
					q1[k].re=temp2re-q1[k].im;
					q1[k].im=temp2im;
				}
			}
		}
	}
	
	unsigned int mmax=4;
 	while(mmax < n) {
		unsigned int istep=mmax << 1;
		int istepnk=istep*inc1;
		int mmaxnk=mmax*inc1;
		Complex *p;
		for(p=data; p < pstop; p += istepnk) {
			Complex *q=p+mmaxnk;
#pragma ivdep			
			for(unsigned int k=0; k < kstop; k += inc2) {
				Real tempre=p[k].re;
				Real tempim=p[k].im;
				p[k].re += q[k].re;
				p[k].im += q[k].im;
				q[k].re=tempre-q[k].re;
				q[k].im=tempim-q[k].im;
			}
		}

		Real wpre=wp->re;
		Real wpim=wp->im*isign;
		Real c=1.0+wpre, s=wpim;
		for(unsigned int m=1; m < mmax-1; m += 2) {
			Real c2=c+c*wpre-s*wpim;
			Real s2=s+s*wpre+c*wpim;
			for(p=data+m*inc1; p < pstop; p += istepnk) {
				Complex *p1=p+inc1,*q=p+mmaxnk,*q1=q+inc1;
#pragma ivdep			
				for(unsigned int k=0; k < kstop; k += inc2) {
					Real tempre=c*q[k].re-s*q[k].im;
					Real tempim=c*q[k].im+s*q[k].re;
					Real temp2re=c2*q1[k].re-s2*q1[k].im;
					Real temp2im=c2*q1[k].im+s2*q1[k].re;
					q[k].re=p[k].re-tempre;
					q[k].im=p[k].im-tempim;
					q1[k].re=p1[k].re-temp2re;
					q1[k].im=p1[k].im-temp2im;
					p[k].re += tempre;
					p[k].im += tempim;
					p1[k].re += temp2re;
					p1[k].im += temp2im;
				}
			}
			c=c2+c2*wpre-s2*wpim;
			s=s2+s2*wpre+c2*wpim;
		}
		for(p=data+mmaxnk-inc1; p < pstop; p += istepnk) {
			Complex *q=p+mmaxnk;
#pragma ivdep			
			for(unsigned int k=0; k < kstop; k += inc2) {
				Real tempre=c*q[k].re-s*q[k].im;
				Real tempim=c*q[k].im+s*q[k].re;
				Real qre=p[k].re-tempre;
				Real qim=p[k].im-tempim;
				p[k].re += tempre;
				p[k].im += tempim;
				q[k].re=qre;
				q[k].im=qim;
			}
		}
		mmax=istep;
		wp++;
	}
}
	
// Return the bit-reversed inverse Fourier transform of a Complex vector.
// Before calling, data must be allocated as Complex[n].
// On entry: data contains the n Complex values;
//           log2n contains the base-2 logarithm of n;
//           isign is the sign (+/- 1) of the phase;
//           [scale is a constant by which the results will be multiplied.]
// On exit:  data contains the n Complex bit-reversed inverse Fourier values.
// Note: To compute a true inverse transform, set scale=1.0/n.

static void fft_brinv(Complex *data, unsigned int log2n, int isign, Real scale)
{
	if(log2n > wpTableSize) fft_init(log2n);
	
	unsigned int n=1 << log2n;
	unsigned int istep=n;
	Complex *pstop=data+n, *wp=wpTable+log2n-1;

 	while(istep > 4) {
		unsigned int mmax=istep >> 1;
		Complex *p;
#pragma ivdep			
		for(p=data; p < pstop; p += istep) {
			Complex *q=p+mmax;
			Real tempre=p->re;
			Real tempim=p->im;
			p->re += q->re;
			p->im += q->im;
			q->re=tempre-q->re;
			q->im=tempim-q->im;
		}
		
		Real wpre=wp->re;
		Real wpim=wp->im*isign;
		Real c=1.0+wpre, s=wpim;
		for(unsigned int m=1; m < mmax-1; m += 2) {
			Real c2=c+c*wpre-s*wpim;
			Real s2=s+s*wpre+c*wpim;
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
			c=c2+c2*wpre-s2*wpim;
			s=s2+s2*wpre+c2*wpim;
		}
#pragma ivdep			
		for(p=data+mmax-1; p < pstop; p += istep) {
			Complex *q=p+mmax;
			Real tempre=p->re-q->re;
			Real tempim=p->im-q->im;
			Real qre=c*tempre-s*tempim;
			tempim=c*tempim+s*tempre;
			p->re += q->re;
			p->im += q->im;
			q->re=qre;
			q->im=tempim;
		}
		istep=mmax;
		wp--;
	}
	
 	if(n > 2) {
		if(isign == 1) {
#pragma ivdep			
			for(Complex *p=data; p < pstop; p += 4) {
				Complex *p1=p+1,*q=p+2,*q1=p+3;
				Real tempre=p->re;
				Real tempim=p->im;
				Real temp2re=p1->re-q1->re;
				Real temp2im=p1->im;
				p->re += q->re;
				p->im += q->im;
				p1->re += q1->re;
				p1->im += q1->im;
				q->re=tempre-q->re;
				q->im=tempim-q->im;
				q1->re=q1->im-temp2im;
				q1->im=temp2re;
			}
		} else {
#pragma ivdep			
			for(Complex *p=data; p < pstop; p += 4) {
				Complex *p1=p+1,*q=p+2,*q1=p+3;
				Real tempre=p->re;
				Real tempim=p->im;
				Real temp2re=q1->re-p1->re;
				Real temp2im=p1->im;
				p->re += q->re;
				p->im += q->im;
				p1->re += q1->re;
				p1->im += q1->im;
				q->re=tempre-q->re;
				q->im=tempim-q->im;
				q1->re=temp2im-q1->im;
				q1->im=temp2re;
			}
		}
	}
	
	if(scale == 1.0) {
		if(n > 1) {
#pragma ivdep			
			for(Complex *p=data; p < pstop; p += 2) {
				Complex *q=p+1;
				Real tempre=p->re;
				Real tempim=p->im;
				p->re += q->re;
				p->im += q->im;
				q->re=tempre-q->re;
				q->im=tempim-q->im;
			}
		}
	} else {
		if(n > 1) {
#pragma ivdep			
			for(Complex *p=data; p < pstop; p += 2) {
				Complex *q=p+1;
				Real tempre=p->re;
				Real tempim=p->im;
				p->re=(p->re+q->re)*scale;
				p->im=(p->im+q->im)*scale;
				q->re=(tempre-q->re)*scale;
				q->im=(tempim-q->im)*scale;
			}
		} else {
			data->re *= scale;
			data->im *= scale;
		}
	}
}

// Return the bit-reversed inverse Fourier transform of nk Complex vectors.
// Before calling, data must be allocated as Complex[nk*n].
// On entry: data contains the n Complex values for each k=0,...,nk-1;
//           log2n contains the base-2 logarithm of n;
//           isign is the sign (+/- 1) of the phase;
//           [inc1 is the stride between the elements of each Complex vector;]
//           [inc2 is the stride between the first elements of the vectors;]
//           [scale is a constant by which the results will be multiplied.]
// On exit:  data contains the n Complex bit-reversed inverse Fourier values 
// for each k=0,...,nk-1.
// Note: To compute a true inverse transform, set scale=1.0/n.

static void mfft_brinv(Complex *data, unsigned int log2n, int isign,
					   unsigned int nk, unsigned int inc1, unsigned int inc2,
					   Real scale)
{
	if(inc1 == 0) inc1=nk;
	unsigned int kstop=nk*inc2;
	
	if(log2n > wpTableSize) fft_init(log2n);
	
	unsigned int n=1 << log2n;
	unsigned int istep=n;
	Complex *pstop=data+n*inc1, *wp=wpTable+log2n-1;

 	while(istep > 4) {
		unsigned int mmax=istep >> 1;
		int istepnk=istep*inc1;
		int mmaxnk=mmax*inc1;
		Complex *p;
		for(p=data; p < pstop; p += istepnk) {
			Complex *q=p+mmaxnk;
#pragma ivdep			
			for(unsigned int k=0; k < kstop; k += inc2) {
				Real tempre=p[k].re;
				Real tempim=p[k].im;
				p[k].re += q[k].re;
				p[k].im += q[k].im;
				q[k].re=tempre-q[k].re;
				q[k].im=tempim-q[k].im;
			}
		}
		
		Real wpre=wp->re;
		Real wpim=wp->im*isign;
		Real c=1.0+wpre, s=wpim;
		for(unsigned int m=1; m < mmax-1; m += 2) {
			Real c2=c+c*wpre-s*wpim;
			Real s2=s+s*wpre+c*wpim;
			for(p=data+m*inc1; p < pstop; p += istepnk) {
				Complex *p1=p+inc1,*q=p+mmaxnk,*q1=q+inc1;
#pragma ivdep			
				for(unsigned int k=0; k < kstop; k += inc2) {
					Real tempre=p[k].re-q[k].re;
					Real tempim=p[k].im-q[k].im;
					Real temp2re=p1[k].re-q1[k].re;
					Real temp2im=p1[k].im-q1[k].im;
					p[k].re += q[k].re;
					p[k].im += q[k].im;
					p1[k].re += q1[k].re;
					p1[k].im += q1[k].im;
					q[k].re=c*tempre-s*tempim;
					q[k].im=c*tempim+s*tempre;
					q1[k].re=c2*temp2re-s2*temp2im;
					q1[k].im=c2*temp2im+s2*temp2re;
				}
			}
			c=c2+c2*wpre-s2*wpim;
			s=s2+s2*wpre+c2*wpim;
		}
		for(p=data+mmaxnk-inc1; p < pstop; p += istepnk) {
			Complex *q=p+mmaxnk;
#pragma ivdep			
			for(unsigned int k=0; k < kstop; k += inc2) {
				Real tempre=p[k].re-q[k].re;
				Real tempim=p[k].im-q[k].im;
				Real qre=c*tempre-s*tempim;
				tempim=c*tempim+s*tempre;
				p[k].re += q[k].re;
				p[k].im += q[k].im;
				q[k].re=qre;
				q[k].im=tempim;
			}
		}
		istep=mmax;
		wp--;
	}
	
 	if(n > 2) {
		int fourinc1=4*inc1;
		if(isign == 1) {
			for(Complex *p=data; p < pstop; p += fourinc1) {
				Complex *p1=p+inc1,*q=p1+inc1,*q1=q+inc1;
#pragma ivdep			
				for(unsigned int k=0; k < kstop; k += inc2) {
					Real tempre=p[k].re;
					Real tempim=p[k].im;
					Real temp2re=p1[k].re-q1[k].re;
					Real temp2im=p1[k].im;
					p[k].re += q[k].re;
					p[k].im += q[k].im;
					p1[k].re += q1[k].re;
					p1[k].im += q1[k].im;
					q[k].re=tempre-q[k].re;
					q[k].im=tempim-q[k].im;
					q1[k].re=q1[k].im-temp2im;
					q1[k].im=temp2re;
				}
			}
		} else {
			for(Complex *p=data; p < pstop; p += fourinc1) {
				Complex *p1=p+inc1,*q=p1+inc1,*q1=q+inc1;
#pragma ivdep			
				for(unsigned int k=0; k < kstop; k += inc2) {
					Real tempre=p[k].re;
					Real tempim=p[k].im;
					Real temp2re=q1[k].re-p1[k].re;
					Real temp2im=p1[k].im;
					p[k].re += q[k].re;
					p[k].im += q[k].im;
					p1[k].re += q1[k].re;
					p1[k].im += q1[k].im;
					q[k].re=tempre-q[k].re;
					q[k].im=tempim-q[k].im;
					q1[k].re=temp2im-q1[k].im;
					q1[k].im=temp2re;
				}
			}
		}
	}
	
	int twoinc1=2*inc1;
	
	if(scale == 1.0) {
		if(n > 1) {
			for(Complex *p=data; p < pstop; p += twoinc1) {
				Complex *q=p+inc1;
#pragma ivdep			
				for(unsigned int k=0; k < kstop; k += inc2) {
					Real tempre=p[k].re;
					Real tempim=p[k].im;
					p[k].re += q[k].re;
					p[k].im += q[k].im;
					q[k].re=tempre-q[k].re;
					q[k].im=tempim-q[k].im;
				}
			}
		}
	} else {
		if(n > 1) {
			for(Complex *p=data; p < pstop; p += twoinc1) {
				Complex *q=p+inc1;
#pragma ivdep			
				for(unsigned int k=0; k < kstop; k += inc2) {
					Real tempre=p[k].re;
					Real tempim=p[k].im;
					p[k].re=(p[k].re+q[k].re)*scale;
					p[k].im=(p[k].im+q[k].im)*scale;
					q[k].re=(tempre-q[k].re)*scale;
					q[k].im=(tempim-q[k].im)*scale;
				}
			}
		} else {
#pragma ivdep			
			for(unsigned int k=0; k < kstop; k += inc2) {
				data[k].re *= scale;
				data[k].im *= scale;
			}
		}
	}
}
