#include "Array.h"
#include "fftw++.h"

// Compile with g++ fftw++example3.cc fftw++.cc -lfftw3

using std::cout;
using Array::array3;

int main() {
  unsigned int n=4, m=5, p=6;
  size_t align=sizeof(Complex);
  
  array3<Complex> f(n,m,p,align);
  
  fft3d Forward3(-1,f);
  fft3d Backward3(1,f);
  
  for(unsigned int i=0; i < n; i++) 
    for(unsigned int j=0; j < m; j++) 
      for(unsigned int k=0; k < p; k++) 
      f(i,j,k)=i+j+k;
	
  cout << f << endl;
  
  cout << endl;
  
  Forward3.fft(f);
  Backward3.fftNormalized(f);
  
  cout << f << endl;
}
