#include "utils.h"

extern "C" void dcft(const int& init,
					 Complex *x, const int& inc1x, const int& inc2x,
					 Complex *y, const int& inc1y, const int& inc2y,
					 const int& n, const int& m, const int& isign,
					 const double& scale, double *aux1, const int& naux1,
					 double *aux2, const int& naux2); 

static int zero=0, one=1;
static double scale=1.0;

void mfft(Complex *data, unsigned int log2n, int isign, unsigned int nk,
		  unsigned int inc1, unsigned int inc2, int)
{
	static int naux,nlast=0, nklast=0, inc1last=0, inc2last=0;	
	static double *aux1,*aux2;
	unsigned int n=1 << log2n;
	isign = -isign;
	
	if(n != nlast || nk != nklast || inc1 != inc1last || inc2 != inc2last) {
		nlast=n;
		nklast=nk;
		inc1last=inc1;
		inc2last=inc2;
		if(n <= 2048) naux=20000;
		else naux=20000+2.28*n;
		if(inc2 == 1 || n >= 252) naux += (2*n+256)*((nk > 64) ? nk : 64);
	
		aux1=new(aux1,naux) Real;
		aux2=new(aux2,naux) Real;
		dcft(one,data,inc1,inc2,data,inc1,inc2,n,nk,isign,scale,
			 aux1,naux,aux2,naux);
	}
	dcft(zero,data,inc1,inc2,data,inc1,inc2,n,nk,isign,scale,
		 aux1,naux,aux2,naux);
}

void fft(Complex *data, unsigned int log2n, int isign, int)
{
	mfft(data,log2n,isign,1,1,1,1);
}
