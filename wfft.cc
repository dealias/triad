#include "options.h"
#include "fft.h"

static ifstream ifwisdom;
static ofstream ofwisdom;

static char *wisdom_name="wisdom.txt";

static fftw_plan *plan=NULL, *planinv=NULL;
static unsigned int nplan=0;

// Create plans for size n
static int init_plan(unsigned int n)
{
	nplan++;
	plan=new(plan,nplan) fftw_plan;
	planinv=new(planinv,nplan) fftw_plan;
	
	ifwisdom.open(wisdom_name);
	fftw_import_wisdom(get_wisdom,ifwisdom);
	ifwisdom.close();
	
	int options=FFTW_MEASURE | FFTW_USE_WISDOM | FFTW_IN_PLACE;
	plan[nplan-1]=fftw_create_plan(n, 1, options);
	planinv[nplan-1]=fftw_create_plan(n, -1, options);
	
	ofwisdom.open(wisdom_name);
	fftw_export_wisdom(put_wisdom,ofwisdom);
	ofwisdom.close();
	
	return nplan;
}

void mfft(Complex *data, unsigned int log2n, int isign, unsigned int nk,
		  unsigned int inc1, unsigned int inc2, Real scale, int)
{
	static unsigned int TableSize=0;
	static unsigned int *Table=NULL;
	static Complex *Work=NULL;
	
	unsigned int n=1 << log2n;
	
	if(inc1 == 0) inc1=nk;
	
	if(n > TableSize) {
		unsigned int nold=TableSize;
		TableSize=n;
		Table=new(Table,TableSize) unsigned int;
		for(unsigned int i=nold; i < TableSize; i++) Table[i]=0;
		Work=new(Work,n) Complex;
	}
	
    if(Table[n-1] == 0) Table[n-1]=init_plan(n);
	
	int j=Table[n-1]-1;
	fftw((isign == 1) ? plan[j] : planinv[j],nk,data,inc1,inc2,Work,0,0);
	
	if(scale != 1.0) {
		unsigned int kstop=nk*inc2;
		if(inc1 == 1) {
			for(unsigned int k=0; k < kstop; k += inc2) {
				Complex *p0=data+k, *pstop=p0+n;
#pragma ivdep			
				for(Complex *p=p0; p < pstop; p++) {
					p->re *= scale;
					p->im *= scale;
				}
			}
		} else {
			Complex *pstop=data+n*inc1;
			for(Complex *p=data; p < pstop; p += inc1) {
#pragma ivdep			
				for(unsigned int k=0; k < kstop; k += inc2) {
					p[k].re *= scale;
					p[k].im *= scale;
				}
			}
		}
	}
}

void fft(Complex *data, unsigned int log2n, int isign, Real scale, int)
{
	mfft(data,log2n,isign,1,1,1,scale);
}
