#include "options.h"
#include "fft.h"

extern "C" void CFTFAX(const int& n, int* ifax, Real *trig);
						
extern "C" void CFFTMLT(Real *ar, Real *ai, Real *work, Real *trigs,
						int *ifax , const int& inc, const int& jump,
						const int& n, const int& lot, const int& isign);

void fft4(Complex *data, unsigned int log4n, int isign)
{
	unsigned int m,k,j;
	static unsigned int phasesize=0, lastsize=0;
	static int lastsign=0;
	static Complex *phase;
	static Real *work,*trigs;
	int ifax[19];

	m=1 << log4n;
	
	if(m != lastsize) {
		if(m > lastsize) {
			phase=new(phase,m*m) Complex;
			work=new(work,4*m*m) Real;
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
	
	CFFTMLT(&data[0].re,&data[0].im,work,trigs,ifax,m,1,m,m,isign);
	
	for(k=0; k < m; k++) {
		Complex *p=data+m*k, *q=data+k;
#pragma ivdep			
        for(j=0; j < k; j++, q += m) {
			Complex temp=p[j];
			p[j]=*q;
			*q=temp;
		}
	}

	if(isign == 1) for(j=0; j < m; j++) {
		Complex *p=data+m*j;
		Complex *phasej=phase+m*j;
#pragma ivdep			
		for(k=0; k < m; k++) p[k] *= phasej[k];
	} else for(j=0; j < m; j++) {
		Complex *p=data+m*j;
		Complex *phasej=phase+m*j;
#pragma ivdep			
		for(k=0; k < m; k++) p[k] *= conj(phasej[k]);
	}
	
	CFFTMLT(&data[0].re,&data[0].im,work,trigs,ifax,m,1,m,m,isign);
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

static double *table,*work;
static int zero=0;
static int nlast=0, isys[1]={0};

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

