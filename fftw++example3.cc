#include "Array.h"
#include "fftw++.h"

// Compile with g++ fftw++example3.cc fftw++.cc -lfftw3

using namespace std;
using namespace Array;

int main() {
  unsigned int n=4, m=5, p=6;
  size_t align=sizeof(Complex);
  
  array3<Complex> g(n,m,p,align);
  
  fft3d Forward3(-1,g);
  fft3d Backward3(1,g);
  
  for(unsigned int i=0; i < n; i++) 
    for(unsigned int j=0; j < m; j++) 
      for(unsigned int k=0; k < p; k++) 
      g(i,j,k)=i+j+k;
	
  cout << g << endl;
  
  cout << endl;
  
  Forward3.fft(g);
  Backward3.fftNormalized(g);
  
  cout << g << endl;
}
