#include "utils.h"

const double twopi=2*PI;

void convolve(Complex *H, Complex *F, Complex *G, unsigned int m, unsigned
			  int log2n);
void convolve_direct(Complex *H, Complex *F, Complex *G, unsigned int m);
	
int main()
{	
	unsigned int i;
	unsigned int log2n=5;
	Complex f[17],g[17];
	Complex h0[]={2,Complex(3,1),Complex(4,-2)};
	
	Complex h[]={-5,Complex(3,1),Complex(4,-2),Complex(-3,1),Complex(0,-2),Complex(0,1),Complex(4,0),Complex(-3,-1),Complex(1,2),Complex(2,1),Complex(3,1)};
	
	unsigned int m=sizeof(h)/sizeof(Complex);

	for(i=0; i < m; i++) f[i]=h[i];
	for(i=0; i < m; i++) g[i]=h[i];
	convolve(f,f,g,m,log2n);
	
	for(i=0; i < m; i++) cout << f[i] << endl;
	cout << endl;
	
	for(i=0; i < m; i++) f[i]=h[i];
	for(i=0; i < m; i++) g[i]=h[i];
	convolve_direct(h,f,g,m);
	for(i=0; i < m; i++) cout << h[i] << endl;
}



