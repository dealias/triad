#include "utils.h"

extern "C" void CFFT(const int& isign, const int& n, double& scale,
					 Complex *x, const int& incx, Complex *y, const int& incy,
					 double *table, const int& ntable, double *work, 
					 const int& nwork);
		 
static double *table,*work;
static int nlast=0,ntable=0,nwork=0,one=1;
static double scale=1.0;

inline void fft0(Complex *data, unsigned int log2n, int isign) {
	unsigned int n=1 << log2n;
	if(n > nlast) {
		nlast=n;
		nwork=12*n;
		work=new(work,nwork) Real;
		ntable=100+8*n;
		table=new(table,ntable) Real;
	}
	CFFT(isign,n,scale,data,one,data,one,table,ntable,work,nwork);
}

void fft_br(Complex *data, unsigned int log2n) {
	fft0(data,log2n,1);
}

void fft_brinv(Complex *data, unsigned int log2n) {
	fft0(data,log2n,-1);
}
