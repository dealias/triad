#include "Array.h"
#include "fftw++.h"

// Compile with g++ fftw++example3.cc fftw++.cc -lfftw3

using std::cout;
using Array::array2;

int main() {
  unsigned int n=4, m=5;
  size_t align=sizeof(Complex);
  
  array2<Complex> f(n,m,align);
  
  fft2d Forward2(-1,f);
  fft2d Backward2(1,f);
  
  for(unsigned int i=0; i < n; i++) 
    for(unsigned int j=0; j < m; j++) 
      f(i,j)=i+j;
	
  cout << f << endl;
  
  cout << endl;
  
  Forward2.fft(f);
  Backward2.fftNormalized(f);
  
  cout << f << endl;
}
