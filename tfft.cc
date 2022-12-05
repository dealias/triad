#include "options.h"
#include "fft.h"

#ifdef _CRAYMPP 

extern "C" void CCFFT(const int& isign, const int& n, const Real& scale,
		      Complex *x, Complex *y, Real *table, Real *work, 
		      const int& isys);
		 
typedef int *pint;
typedef Real *pReal;

void fft(Complex *data, size_t log2n, int isign, Real scale, int)
{
  static int TableSize=0;
  static size_t *nTable=NULL;
  static Real **table=NULL,**work=NULL;
  const int isys=0;
  const int zero=0;
  int j;
	
  size_t n=1 << log2n;
	
  for(j=0; j < TableSize; j++) if(n == nTable[j]) break;
	
  if(j == TableSize) {
    TableSize++;
    nTable=new(nTable,TableSize) size_t;
    table=new(table,TableSize) pReal;
    work=new(work,TableSize) pReal;
    nTable[j]=n;
    table[j]=new Real[2*n];
    work[j]=new Real[4*n];
    CCFFT(zero,n,scale,data,data,table[j],work[j],isys);
  }
	
  CCFFT(isign,n,scale,data,data,table[j],work[j],isys);
}

// Return the Fourier transform of nk Complex vector's.
// Non-vectorizing version.
// Before calling, data must be allocated as Complex[nk*n].
// On entry: data contains the n Complex values for each k=0,...,nk-1.
//           log2n contains the base-2 logarithm of n.
//           isign is +1 for a forward transform, -1 for an inverse transform.
//           inc1 is the stride between the elements of each Complex vector.
//           inc2 is the stride between first elements of the vectors.
// On exit:  data contains the n Complex Fourier values for each k=0,...,nk-1.
// Note: When computing an inverse transform, the result must be divided by n.

void mfft(Complex *data, size_t log2n, int isign, size_t nk,
	  size_t inc1, size_t inc2, Real scale, int)
{
  static size_t data2size=0;
  static Complex *data2=NULL;
	
  if(inc1 == 1) {
    Complex *pstop=data+nk*inc2;
    for(Complex *p=data; p < pstop; p += inc2) fft(p,log2n,isign,scale);
    return;
  }
	
  if(inc1 == 0) inc1=nk;
	
  size_t n=1 << log2n;
  if(n > data2size) data2=new(data2,data2size=n) Complex;

  size_t kstop=nk*inc2;
  for(size_t k=0; k < kstop; k += inc2) {
    Complex *p,*p0=data+k;
    Complex *q,*qstop=data2+n;
    for(p=p0, q=data2; q < qstop; p += inc1) *(q++)=*p;
    fft(data2,log2n,isign,scale);
    for(p=p0, q=data2; q < qstop; p += inc1) *p=*(q++);
  }
}

void mrcfft(Complex *data, size_t log2n, int isign, size_t nk,
	    size_t inc1, size_t inc2, Real scale, int bitreverse)
{		 
  assert(inc1 == 1);
  mrcfft0(data,log2n,isign,nk,inc1,inc2,scale,bitreverse);
}

void mcrfft(Complex *data, size_t log2n, int isign, size_t nk,
	    size_t inc1, size_t inc2, Real scale, int bitreverse)
{		 
  assert(inc1 == 1);
  mcrfft0(data,log2n,isign,nk,inc1,inc2,scale,bitreverse);
}

#endif
