#include "utils.h"

extern "C" void drcft(const int& init, Complex *x, const int& inc2x,
				  Complex *y, const int& inc2y, const int& n,
				  const int& m, const int& isign, const double& scale,
				  double *aux1, const int& naux1,
				  double *aux2, const int& naux2); 

extern "C" void dcrft(const int& init, Complex *x, const int& inc2x,
				  Complex *y, const int& inc2y, const int& n,
				  const int& m, const int& isign, const double& scale,
				  double *aux1, const int& naux1,
				  double *aux2, const int& naux2);  

static int mone=-1, zero=0, one=1;

void init_aux(int n, double *& aux1, int& naux1, double *& aux2, int& naux2)
{
	if(n <= 4096) {naux1=22000; naux2=20000;}
	else {naux1=20000+1.64*n; naux2=20000+1.14*n;}
	
	aux1=new(aux1,naux1) Real;
	aux2=new(aux2,naux2) Real;
}
	 

void rfft_br(Complex *data, unsigned int log2n) {
	static int naux1,naux2,nlast=0;	
	static double *aux1,*aux2;
	static double scale=1.0;
	unsigned int n=1 << (log2n+1);
	
	if(n != nlast) {
		nlast=n;
		init_aux(n,aux1,naux1,aux2,naux2);
		drcft(one,data,zero,data,zero,n,one,one,scale,aux1,naux1,aux2,naux2);
	}
	drcft(zero,data,zero,data,zero,n,one,one,scale,aux1,naux1,aux2,naux2);
}

void rfft_brinv(Complex *data, unsigned int log2n) {
	static int naux1,naux2,nlast=0;	
	static double *aux1,*aux2;
	static double scale=0.5;
	unsigned int n=1 << (log2n+1);
	
	if(n != nlast) {
		nlast=n;
		init_aux(n,aux1,naux1,aux2,naux2);
		dcrft(one,data,zero,data,zero,n,one,mone,scale,aux1,naux1,aux2,naux2);
	}
	dcrft(zero,data,zero,data,zero,n,one,mone,scale,aux1,naux1,aux2,naux2);
}

