#include <sys/machd.h>

#include "options.h"
#include "fft.h"

static double *table,*work;
static int zero=0, two=2;
static int nlast=0, isys[1]={0};

extern "C" void CFTFAX(const int& n, int* ifax, Real *trigs);
						
extern "C" void CFFTMLT(Real *ar, Real *ai, Real *work, Real *trigs,
						int *ifax , const int& inc, const int& jump,
						const int& n, const int& lot, const int& isign);

extern "C" void FLOWMARK(char *);

void fft4(Complex *data, unsigned int log4n, int isign)
{
	unsigned int m,k,j;
	static unsigned int lastsize=0;
	static Complex *phase, *work;
	static Real *trigs;
	static int ifax[19];

	m=1 << log4n;
	int twom=2*m;
	
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
	
//#pragma _CRI taskloop private(j) value(data,work,lot,trigs,ifax,twom,two,m,isign)
	{
		for(int j=0; j < P; j++) {
			int off=j*lot;
			off=0;
			lot=m;
			CFFTMLT(&data[off].re,&data[off].im,work+2*off,trigs,ifax,
					twom,two,m,lot,isign);
		}
	}
	if(remainder) {
		int off=P*lot;
		CFFTMLT(&data[off].re,&data[off].im,work+2*off,trigs,ifax,
				twom,two,m,remainder,isign);
	}
#endif	
	FLOWMARK("inner1");
	CFFTMLT(&data[0].re,&data[0].im,&work[0].re,trigs,ifax,twom,two,m,m,isign);
	FLOWMARK("");
	 
#if 0	
	FLOWMARK("inner2");
	if(isign == 1) {
		for(int j=0; j < m; j++) {
			int mj=m*j;
			Complex *temp=work+mj;
			Complex *p=data+mj;
			Complex *c=phase+mj;
			Complex *q=data+j;
			p[j] *= c[j];
#pragma ivdep			
			for(int k=0; k < j; k++) {
				temp[k] = p[k];
				p[k] = q[m*k]*c[k];
				q[m*k]=temp[k]*c[k];
			}
		}
	} else {
		for(int j=0; j < m; j++) {
			int mj=m*j;
			Complex *temp=work+mj;
			Complex *p=data+mj;
			Complex *c=phase+mj;
			Complex *q=data+j;
			p[j] *= conj(c[j]);
#pragma ivdep			
			for(int k=0; k < j; k++) {
				temp[k] = p[k];
				p[k] = q[m*k]*conj(c[k]);
				q[m*k]=temp[k]*conj(c[k]);
			}
		}
	}
    FLOWMARK("");
#endif	

	FLOWMARK("inner2");
	for(k=0; k < m; k++) {
		Complex *p=data+m*k, *q=data+k;
#pragma ivdep			
        for(j=0; j < k; j++, q += m) {
			Complex temp=p[j];
			p[j]=*q;
			*q=temp;
		}
	}
    FLOWMARK("");

	FLOWMARK("inner3");
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
    FLOWMARK("");
	
//#pragma _CRI taskloop private(j) value(m,data,work,phase)
	
    FLOWMARK("inner4");
	CFFTMLT(&data[0].re,&data[0].im,&work[0].re,trigs,ifax,twom,two,m,m,isign);
	FLOWMARK("");
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

