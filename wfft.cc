#include "options.h"
#include <fstream>
#include <fftw.h>
#include <rfftw.h>
#include "fft.h"

inline void fftw_export_wisdom(void (*emitter)(char c, ofstream& s),
			       ofstream& s)
{
  fftw_export_wisdom((void (*) (char, void *)) emitter, (void *) &s);
}

inline fftw_status fftw_import_wisdom(char (*g)(ifstream& s), ifstream &s)
{
  return fftw_import_wisdom((int (*) (void *)) g, (void *) &s);
}

inline void put_wisdom(char c, ofstream& s)
{
  s.put(c);
}

inline char get_wisdom(ifstream& s)
{
  return s.get();
}

static ifstream ifwisdom;
static ofstream ofwisdom;

static char *wisdom_name="wisdom.txt";

static fftw_plan *plan=NULL, *planinv=NULL;
static unsigned int nplan=0;

static rfftwnd_plan *crplan=NULL;
static unsigned int ncrplan=0;

static rfftwnd_plan *rcplan=NULL;
static unsigned int nrcplan=0;

// Create plans for size n
static int init_plan(unsigned int n, int inc1)
{
  nplan++;
  plan=new(plan,nplan) fftw_plan;
  planinv=new(planinv,nplan) fftw_plan;
	
  ifwisdom.open(wisdom_name);
  fftw_import_wisdom(get_wisdom,ifwisdom);
  ifwisdom.close();
	
  int options=FFTW_MEASURE | FFTW_USE_WISDOM | FFTW_IN_PLACE;
  plan[nplan-1]=fftw_create_plan(n, (fftw_direction) 1, options);
  planinv[nplan-1]=fftw_create_plan(n, (fftw_direction) -1, options);
	
  ofwisdom.open(wisdom_name);
  fftw_export_wisdom(put_wisdom,ofwisdom);
  ofwisdom.close();
	
  return nplan;
}

// Create complex_to_real plans for size n
static int init_crplan(unsigned int n)
{
  ncrplan++;
  crplan=new(crplan,ncrplan) rfftwnd_plan;
	
  ifwisdom.open(wisdom_name);
  fftw_import_wisdom(get_wisdom,ifwisdom);
  ifwisdom.close();
	
  int options=FFTW_MEASURE | FFTW_USE_WISDOM | FFTW_IN_PLACE;
  int size[1]={n};
  crplan[ncrplan-1]=rfftwnd_create_plan(1,size,FFTW_COMPLEX_TO_REAL,options);
	
  ofwisdom.open(wisdom_name);
  fftw_export_wisdom(put_wisdom,ofwisdom);
  ofwisdom.close();
	
  return ncrplan;
}

// Create real_to_complex plans for size n
static int init_rcplan(unsigned int n)
{
  nrcplan++;
  rcplan=new(rcplan,nrcplan) rfftwnd_plan;
	
  ifwisdom.open(wisdom_name);
  fftw_import_wisdom(get_wisdom,ifwisdom);
  ifwisdom.close();
	
  int options=FFTW_MEASURE | FFTW_USE_WISDOM | FFTW_IN_PLACE;
  int size[1]={n};
  rcplan[nrcplan-1]=rfftwnd_create_plan(1,size,FFTW_REAL_TO_COMPLEX,options);
	
  ofwisdom.open(wisdom_name);
  fftw_export_wisdom(put_wisdom,ofwisdom);
  ofwisdom.close();
	
  return nrcplan;
}

void scalefft(Complex *data, unsigned int n, unsigned int nk,
	      unsigned int inc1, unsigned int inc2, Real scale)
{
  unsigned int kstop=nk*inc2;
  if(inc1 == 1) {
    for(unsigned int k=0; k < kstop; k += inc2) {
      Complex *p0=data+k, *pstop=p0+n;
      //#pragma ivdep			
      for(Complex *p=p0; p < pstop; p++) {
	p->re *= scale;
	p->im *= scale;
      }
    }
  } else {
    Complex *pstop=data+n*inc1;
    for(Complex *p=data; p < pstop; p += inc1) {
      //#pragma ivdep			
      for(unsigned int k=0; k < kstop; k += inc2) {
	p[k].re *= scale;
	p[k].im *= scale;
      }
    }
  }
}

void signscalefft(Complex *data, unsigned int n, int isign, unsigned int nk,
		  unsigned int inc1, unsigned int inc2, Real scale)
{
  if(scale == 1.0) {
    if(isign != 1.0) {
      unsigned int kstop=nk*inc2;
      if(inc1 == 1) {
	for(unsigned int k=0; k < kstop; k += inc2) {
	  Complex *p0=data+k, *pstop=p0+n;
	  //#pragma ivdep			
	  for(Complex *p=p0; p < pstop; p++) {
	    p->im=-p->im;
	  }
	}
      } else {
	Complex *pstop=data+n*inc1;
	for(Complex *p=data; p < pstop; p += inc1) {
	  //#pragma ivdep			
	  for(unsigned int k=0; k < kstop; k += inc2) {
	    p[k].im=-p[k].im;
	  }
	}
      }
    }
  } else {
    if(isign == 1.0) scalefft(data,n,nk,inc1,inc2,scale);
    else {
      unsigned int kstop=nk*inc2;
      if(inc1 == 1) {
	for(unsigned int k=0; k < kstop; k += inc2) {
	  Complex *p0=data+k, *pstop=p0+n;
	  //#pragma ivdep			
	  for(Complex *p=p0; p < pstop; p++) {
	    p->re *= scale;
	    p->im *= -scale;
	  }
	}
      } else {
	Complex *pstop=data+n*inc1;
	for(Complex *p=data; p < pstop; p += inc1) {
	  //#pragma ivdep			
	  for(unsigned int k=0; k < kstop; k += inc2) {
	    p[k].re *= scale;
	    p[k].im *= -scale;
	  }
	}
      }
    }
  }
}

void mfft(Complex *data, unsigned int log2n, int isign, unsigned int nk,
	  unsigned int inc1, unsigned int inc2, Real scale, int)
{
  static unsigned int TableSize=0;
  static unsigned int *Table=NULL;
	
  unsigned int n=1 << log2n;
	
  if(inc1 == 0) inc1=nk;
	
  if(n > TableSize) {
    unsigned int nold=TableSize;
    TableSize=n;
    Table=new(Table,TableSize) unsigned int;
    for(unsigned int i=nold; i < TableSize; i++) Table[i]=0;
  }
	
  if(Table[n-1] == 0) Table[n-1]=init_plan(n,inc1);
  int j=Table[n-1]-1;
	
  fftw((isign == 1) ? plan[j] : planinv[j],nk,
       (fftw_complex *) data,inc1,inc2,NULL,1,1);
	
  if(scale != 1.0) scalefft(data,n,nk,inc1,inc2,scale);
}

void fft(Complex *data, unsigned int log2n, int isign, Real scale, int)
{
  mfft(data,log2n,isign,1,1,1,scale);
}

// Return the real inverse Fourier transform of nk Complex vectors, of
// length n/2+1, corresponding to the non-negative part of the frequency
// spectrum. Before calling, data must be allocated as Complex[nk*(n/2+1)].
// On entry: data contains the n/2+1 Complex Fourier transform values for
// each k=0,...,nk-1, 
//           log2n contains the base-2 logarithm of n;
//           isign is the sign (+/- 1) of the phase;
//           nk is the number of Complex vectors;
//           [inc1 is the stride between the elements of each Complex vector;]
//           [inc2 is the stride between first elements of the vectors;]
//           [scale is a constant by which the results will be multiplied;]
//           [bitreverse is 0 for a true fft of data (default);
//                         +1 for a fft of bit-reversed data (faster);
//                         -1 for a bit-reversed fft of data (faster).]
// On exit:  data contains the nk*n real inverse Fourier transform values
// stored as Complex arrays of length n/2 for each k=0,...,nk-1.
//
// Note: To compute a true inverse transform, set scale=1.0/n.

void mcrfft(Complex *data, unsigned int log2n, int isign, unsigned int nk,
	    unsigned int inc1, unsigned int inc2, Real scale, int)
{		 
  static unsigned int TableSize=0;
  static unsigned int *Table=NULL;
	
  unsigned int n=1 << log2n;
  if(log2n == 0) return;
  if(inc1 == 0) inc1=nk;
	
  signscalefft(data,n/2+1,isign,nk,inc1,inc2,scale);

  if(n > TableSize) {
    unsigned int nold=TableSize;
    TableSize=n;
    Table=new(Table,TableSize) unsigned int;
    for(unsigned int i=nold; i < TableSize; i++) Table[i]=0;
  }
	
  if(Table[n-1] == 0) Table[n-1]=init_crplan(n);
  int j=Table[n-1]-1;
	
  rfftwnd_complex_to_real(crplan[j],nk,(fftw_complex *) data,inc1,inc2,
			  NULL,1,1);
}

// Return the Fourier transform of nk real vectors, each of length n.
// Before calling, data must be allocated as Complex[nk*(n/2+1)].
// On entry: data contains the nk*n real values stored as Complex
// arrays of length n/2 for each k=0,...,nk-1;
//           log2n contains the base-2 logarithm of n;
//           isign is the sign (+/- 1) of the phase;
//           nk is the number of Complex vectors;
//           [inc1 is the stride between the elements of each Complex vector;]
//           [inc2 is the stride between first elements of the vectors;]
//           [scale is a constant by which the results will be multiplied;]
//           [bitreverse is 0 for a true fft of data (default);
//                         +1 for a fft of bit-reversed data (faster);
//                         -1 for a bit-reversed fft of data (faster).]
// On exit:  data contains the nk*(n/2+1) Complex Fourier values.

void mrcfft(Complex *data, unsigned int log2n, int isign, unsigned int nk,
	    unsigned int inc1, unsigned int inc2, Real scale, int)
{		 
  static unsigned int TableSize=0;
  static unsigned int *Table=NULL;
	
  unsigned int n=1 << log2n;
  if(log2n == 0) return;
  if(inc1 == 0) inc1=nk;
	
  if(n > TableSize) {
    unsigned int nold=TableSize;
    TableSize=n;
    Table=new(Table,TableSize) unsigned int;
    for(unsigned int i=nold; i < TableSize; i++) Table[i]=0;
  }
	
  if(Table[n-1] == 0) Table[n-1]=init_rcplan(n);
  int j=Table[n-1]-1;
	
  int inc2real=((inc2 == 1) ? inc2 : 2*inc2);
	
  rfftwnd_real_to_complex(rcplan[j],nk,(fftw_real *) data,inc1,inc2real,
			  NULL,1,1);
	
  signscalefft(data,n/2+1,-isign,nk,inc1,inc2,scale);
}

