#include "utils.h"

extern "C" void CCFFT(const int& isign, const int& n, double& scale,
					 Complex *x, Complex *y, double *table, double *work, 
					 int* isys);
		 
extern "C" void SCFFT(const int& isign, const int& n, double& scale,
					  Complex *x, Complex *y, double *table, double *work, 
					  int* isys);

extern "C" void CSFFT(const int& isign, const int& n, double& scale,
					  Complex *x, Complex *y, double *table, double *work, 
					  int* isys);

static double *table,*work;
static double *table2,*work2;
static int minus1=-1, zero=0, one=1;
static int nlast=0, nlast2=0, isys[1]={0};

inline void ccfft(Complex *data, unsigned int n, int isign, double scale) {
	if(n != nlast) {
		nlast=n;
		table=new(table,100+8*n) Real;
		work=new(work,8*n) Real;
		CCFFT(zero,n,scale,data,data,table,work,isys);
	}
	CCFFT(isign,n,scale,data,data,table,work,isys);
}

void fft_br(Complex *data, unsigned int log2n) {
	unsigned int n=1 << log2n;
	ccfft(data,n,1,1.0);
}

void fft_brinv(Complex *data, unsigned int log2n) {
	unsigned int n=1 << log2n;
	ccfft(data,n,-1,1.0/n);
}

void rfft_br(Complex *data, unsigned int log2n) {
	unsigned int n=1 << (log2n+1);
	double scale=1.0;
	if(n != nlast2) {
		nlast2=n;
		table2=new(table2,100+4*n) Real;
		work2=new(work2,4+4*n) Real;
		SCFFT(zero,n,scale,data,data,table2,work2,isys);
	}
	SCFFT(one,n,scale,data,data,table2,work2,isys);
}

void rfft_brinv(Complex *data, unsigned int log2n) {
	unsigned int n=1 << (log2n+1);
	double scale=0.5;
	if(n != nlast2) {
		nlast2=n;
		table2=new(table2,100+4*n) Real;
		work2=new(work2,4+4*n) Real;
		CSFFT(zero,n,scale,data,data,table2,work2,isys);
	}
	CSFFT(minus1,n,scale,data,data,table2,work2,isys);
}

