#include <sys/machd.h>

#include "options.h"
#include "fft.h"

static double *table,*work;
static int zero=0;
static int nlast=0, isys[1]={0};

extern "C" void CFTFAX(const int& n, int* ifax, Real *trigs);
						
extern "C" void CFFTMLT(Real *ar, Real *ai, Real *work, Real *trigs,
						int *ifax , const int& inc, const int& jump,
						const int& n, const int& lot, const int& isign);

void fft4(Complex *data, unsigned int log4n, int isign)
{
	unsigned int m,k,j;
	static unsigned int lastsize=0;
	static Complex *phase, *work;
	static Real *trigs;
	static int ifax[19];

	m=1 << log4n;
	int inc=2*m;
	int jump=2;
	
	if(m != lastsize) {
		if(m > lastsize) {
			phase=new(phase,m*m) Complex;
			work=new(work,2*m*m) Complex;
			trigs=new(trigs,2*m) Real;
		}
		lastsize=m;
		Complex factor=exp(twopi*I/(m*m));
		Complex kphase=1.0;
		for(k=0; k < m; k++) {
			Complex *p=phase+m*k;
			p[0]=1.0;
			for(j=1; j < m; j++) p[j]=p[j-1]*kphase;
			kphase *= factor;
		}
		CFTFAX(m,ifax,trigs);
	}

#if 0	
//	int P=min(m/64,NCPU); // Request this many processors.
	int P=1;
	int lot=m/P;
	int remainder=m-P*lot;
	
//#pragma _CRI taskloop private(j) value(data,work,lot,trigs,ifax,jump,inc,m,isign)
	{
		for(int j=0; j < P; j++) {
			int off=j*lot;
			off=0;
			lot=m;
			CFFTMLT(&data[off].re,&data[off].im,work+2*off,trigs,ifax,
					inc,jump,m,lot,isign);
		}
	}
	if(remainder) {
		int off=P*lot;
		CFFTMLT(&data[off].re,&data[off].im,work+2*off,trigs,ifax,
				inc,jump,m,remainder,isign);
	}
#endif	

	CFFTMLT(&data[0].re,&data[0].im,&work[0].re,trigs,ifax,inc,jump,m,m,isign);
	
	for(int i=1; i < m; i++) {
		Complex *p=data+i*m;
		Complex *q=data+i;
		int m1=m+1;
#pragma ivdep			
		for(int k=0; k < m-i; k++) {
			int offset=m1*k;
			Complex temp=p[offset];
			p[offset]=q[offset];
			q[offset]=temp;
		}
	}
	 
	if(isign == 1) {
#pragma ivdep			
		for(i=0; i < m*m; i++) data[i] *= phase[i];
	} else {
#pragma ivdep			
		for(i=0; i < m*m; i++) data[i] *= conj(phase[i]);
	}
	
//#pragma _CRI taskloop private(j) value(m,data,work,phase)
	CFFTMLT(&data[0].re,&data[0].im,&work[0].re,trigs,ifax,inc,jump,m,m,isign);
}

extern "C" void CCFFT(const int& isign, const int& n, double& scale,
					 Complex *x, Complex *y, double *table, double *work, 
					 int* isys);
		 
extern "C" void SCFFT(const int& isign, const int& n, double& scale,
					  Complex *x, Complex *y, double *table, double *work, 
					  int* isys);

extern "C" void CSFFT(const int& isign, const int& n, double& scale,
					  Complex *x, Complex *y, double *table, double *work, 
					  int* isys);

inline void ccfft(Complex *data, unsigned int n, int isign, double scale)
{
	if(n != nlast) {
		nlast=n;
		table=new(table,100+8*n) Real;
		work=new(work,8*n) Real;
		CCFFT(zero,n,scale,data,data,table,work,isys);
	}
	CCFFT(isign,n,scale,data,data,table,work,isys);
}

void fft_br(Complex *data, unsigned int log2n)
{
	unsigned int n=1 << log2n;
	ccfft(data,n,1,1.0);
}

void fft_brinv(Complex *data, unsigned int log2n)
{
	unsigned int n=1 << log2n;
	ccfft(data,n,-1,1.0/n);
}

