#include "options.h"
#include "fft.h"

extern "C" void dcft(const int& init,
					 Complex *x, const int& inc1x, const int& inc2x,
					 Complex *y, const int& inc1y, const int& inc2y,
					 const int& n, const int& m, const int& isign,
					 const double& scale, double *aux1, const int& naux1,
					 double *aux2, const int& naux2); 

static int zero=0, one=1;
static double scale=1.0;

void mfft(Complex *data, unsigned int log2n, int isign, unsigned int nk,
		  unsigned int inc1, unsigned int inc2, int)
{
	int j,naux;
	static int TableSize=0;
	static int *nTable,nkTable,*nauxTable;
	static double **aux1[2],**aux2[2];
	unsigned int n=1 << log2n;
	isign = -isign;
	
	for(j=0; j < TableSize; j++) 
		if(n == nTable[j] && nk == nkTable[j]) break;
	
    if(j == TableSize) {
		TableSize++;
		nTable=new(nTable,TableSize);
		nkTable=new(nkTable,TableSize);
		aux1=new(aux1,TableSize);
		aux2=new(aux2,TableSize);
		
		if(n <= 2048) naux=20000;
		else naux=20000+2.28*n;
		naux += (2*n+256)*((nk > 64) ? nk : 64);
		nTable[j]=n;
		nkTable[j]=nk;
		nauxTable[j]=naux;
	
		for(int isign=-1; isign <= 1; isign += 2) {
			int i=(isign == -1);
			aux1[j][i]=new(aux1[j][i],naux) Real;
			aux2[j][i]=new(aux2[j][i],naux) Real;
			dcft(one,data,inc1,inc2,data,inc1,inc2,n,nk,isign,scale,
				 aux1[j][i],naux,aux2[j][i],naux);
		}
	}
	
	int i=(isign == -1);
	naux=nauxTable[j];
	dcft(zero,data,inc1,inc2,data,inc1,inc2,n,nk,isign,scale,
		 aux1[j][i],naux,aux2[j][i],naux);
}

void fft(Complex *data, unsigned int log2n, int isign, int)
{
	mfft(data,log2n,isign,1,1,1,1);
}
