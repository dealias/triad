#include <sys/machd.h>

#include "options.h"
#include "fft.h"

extern "C" void CFTFAX(const int& n, int* ifax, Real *trigs);
						
extern "C" void CFFTMLT(Real *ar, Real *ai, Real *work, Real *trigs,
						int *ifax , const int& inc, const int& jump,
						const int& n, const int& lot, const int& isign);

void mfft(Complex *data, unsigned int log2n, int isign, unsigned int nk,
		  unsigned int inc1, unsigned int inc2, int)
{
	unsigned int n=1 << log2n;
	static unsigned int nlast=0;
	static Real *trigs;
	static int ifax[19];
	static Real *work;
	
	if(n != nlast) {
		nlast=n;
		work=new(work,4*n*nk) Real;
		trigs=new(trigs,2*n) Real;
		CFTFAX(n,ifax,trigs);
	}
	
	int inc=2*inc1;
	int jump=2*inc2; // Should be odd. JCB
	CFFTMLT(&data[0].re,&data[0].im,work,trigs,ifax,inc,jump,n,nk,isign);
}

extern "C" void CCFFT(const int& isign, const int& n, Real& scale,
					 Complex *x, Complex *y, Real *table, Real *work, 
					 int* isys);
		 
void fft(Complex *data, unsigned int log2n, int isign, int)
{
	static double *table,*work;
	unsigned int n=1 << log2n;
	static int nlast=0, isys[1]={0};
	const int zero=0;

	Real scale=1.0;
	if(n != nlast) {
		nlast=n;
		table=new(table,100+8*n) Real;
		work=new(work,8*n) Real;
		CCFFT(zero,n,scale,data,data,table,work,isys);
	}
	CCFFT(isign,n,scale,data,data,table,work,isys);
}



