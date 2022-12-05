#include "convolution.h"

// Compile with: g++ -I ~/tri convtest.cc fftw++.cc -lfftw3

using std::cout;
using std::endl;

int main()
{	
  size_t i;
  size_t n=32;
  size_t np=n/2+1;
  
  Complex f[np],g[np],h[np];
  Complex d[]={-5,Complex(3,1),Complex(4,-2),Complex(-3,1),Complex(0,-2),Complex(0,1),Complex(4,0),Complex(-3,-1),Complex(1,2),Complex(2,1),Complex(3,1)};
	
  size_t m=sizeof(d)/sizeof(Complex);
  convolution convolve(n,m);

  for(i=0; i < m; i++) f[i]=d[i];
  for(i=0; i < m; i++) g[i]=d[i];
  convolve.fft(h,f,g);
	
  for(i=0; i < m; i++) cout << h[i] << endl;
  cout << endl;
	
  for(i=0; i < m; i++) f[i]=d[i];
  for(i=0; i < m; i++) g[i]=d[i];
  convolve.direct(h,f,g);
  for(i=0; i < m; i++) cout << h[i] << endl;
}
