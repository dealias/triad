#include <sys/machd.h>

#include "options.h"
#include "fft.h"

static double *table,*work;
static int zero=0;
static int nlast=0, isys[1]={0};

extern "C" void CFTFAX(const int& n, int* ifax, double *trigs);
						
extern "C" void CFFTMLT(double *ar, double *ai, double *work, double *trigs,
						int *ifax , const int& inc, const int& jump,
						const int& n, const int& lot, const int& isign);

// Note: fft4 does not work under optimization level 1 or 2.
// Compile only with optimzation level 0 or 3.

void fft4(Complex *data, unsigned int log4n, int isign)
{
	unsigned int m,k,j;
	static unsigned int lastsize=0;
	static Complex *phase, *work;
	static double *trigs;
	static int ifax[19];

	m=1 << log4n;
	int inc=2*m;
	int jump=2;
	
	if(m != lastsize) {
		if(m > lastsize) {
			phase=new(phase,m*m) Complex;
			work=new(work,2*m*m) Complex;
			trigs=new(trigs,2*m) double;
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

	int P=m/64; // Request this many processors.
    if(P > NCPU) P=NCPU;
	if(P < 1) P=1;
	int lot=m/P;
	int last=lot+m-P*lot;
	{
#pragma _CRI taskloop private(j) value(P,data,work,trigs,ifax,inc,jump,m, \
lot,isign,last)
		for(int j=0; j < P; j++) {
			int off=j*lot;
			CFFTMLT(&data[off].re,&data[off].im,&work[2*off].re,trigs,ifax,
					inc,jump,m,(j == P-1) ? lot : last,isign);
		}
	}
	
	int m1=m+1;
	if(isign == 1) {
		for(int i=1; i < m; i++) {
			Complex *p=data+i*m;
			Complex *q=data+i;
			Complex *c=phase+i;
#pragma ivdep			
			for(int k=0; k < m-i; k++) {
				int offset=m1*k;
				Complex factor=c[offset];
				Complex temp=p[offset]*factor;
				p[offset]=q[offset]*factor;
				q[offset]=temp;
			}
		}
#pragma ivdep			
		for(i=0; i < m*m; i += m1) data[i] *= phase[i];
	} else {
		for(int i=1; i < m; i++) {
			Complex *p=data+i*m;
			Complex *q=data+i;
			Complex *c=phase+i;
#pragma ivdep			
			for(int k=0; k < m-i; k++) {
				int offset=m1*k;
				Complex factor=c[offset];
				Complex temp=p[offset]*conj(factor);
				p[offset]=q[offset]*conj(factor);
				q[offset]=temp;
			}
		}
#pragma ivdep			
		for(i=0; i < m*m; i += m1) data[i] *= conj(phase[i]);
	}
	
	{
#pragma _CRI taskloop private(j) value(P,data,work,trigs,ifax,inc,jump,m, \
lot,isign,last)
		for(int j=0; j < P; j++) {
			int off=j*lot;
			CFFTMLT(&data[off].re,&data[off].im,&work[2*off].re,trigs,ifax,
					inc,jump,m,(j == P-1) ? lot : last,isign);
		}
	}
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
		table=new(table,100+8*n) double;
		work=new(work,8*n) double;
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

