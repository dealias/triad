 #include "options.h"
#include "fft.h"

extern "C" void CFTFAX(const int& n, int* ifax, Real *trigs);
						
extern "C" void CFFTMLT(Real *ar, Real *ai, Real *work, Real *trigs,
						int *ifax , const int& inc, const int& jump,
						const int& n, const int& lot, const int& isign);

typedef int *pint;
typedef Real *pReal;

void mfft(Complex *data, unsigned int log2n, int isign, unsigned int nk,
		  unsigned int inc1, unsigned int inc2, int)
{
	static int TableSize=0;
	static unsigned int *nTable=NULL, *nkTable=NULL;
	static int **ifax=NULL;
	static Real **trigs=NULL,**work=NULL;
	int j;
	
	unsigned int n=1 << log2n;
	
	for(j=0; j < TableSize; j++) if(n == nTable[j] && nk == nkTable[j]) break;
	
    if(j == TableSize) {
		TableSize++;
		nTable=new(nTable,TableSize) unsigned int;
		nkTable=new(nkTable,TableSize) unsigned int;
		ifax=new(ifax,TableSize) pint;
		trigs=new(trigs,TableSize) pReal;
		work=new(work,TableSize) pReal;
		
		nTable[j]=n;
		nkTable[j]=nk;
		ifax[j]=new int[19];
		trigs[j]=new Real[2*n];
		work[j]=new Real[4*n*nk];
		CFTFAX(n,ifax[j],trigs[j]);
	}
	
	int inc=2*inc1;
	int jump=2*inc2;
	CFFTMLT(&data[0].re,&data[0].im,work[j],trigs[j],ifax[j],inc,jump,n,nk,
			isign);
}

extern "C" void CCFFT(const int& isign, const int& n, Real& scale,
					 Complex *x, Complex *y, Real *table, Real *work, 
					 const int& isys);
		 
void fft(Complex *data, unsigned int log2n, int isign, int)
{
	static int TableSize=0;
	static unsigned int *nTable=NULL;
	static Real **table=NULL,**work=NULL;
	const int isys=0;
	const int zero=0;
	const Real scale=1.0;
	int j;
	
	unsigned int n=1 << log2n;
	
	for(j=0; j < TableSize; j++) if(n == nTable[j]) break;
	
    if(j == TableSize) {
		TableSize++;
		nTable=new(nTable,TableSize) unsigned int;
		table=new(table,TableSize) pReal;
		work=new(work,TableSize) pReal;
		nTable[j]=n;
		table[j]=new Real[100+8*n];
		work[j]=new Real[8*n];
		CCFFT(zero,n,scale,data,data,table[j],work[j],isys);
	}
	
	CCFFT(isign,n,scale,data,data,table[j],work[j],isys);
}



