#include "fftw++.h"

// link with -lfftw3

bool fftw::Wise=false;
const char *fftw::WisdomName="wisdom.txt";
ifstream fftw::ifWisdom;
ofstream fftw::ofWisdom;

using namespace std;

int main() {
  unsigned int n=4;
  Complex *f=fftnew(n);
  
  for(unsigned int i=0; i < n; i++) f[i]=i;
	
  fft1d Forward(n,-1);
  Forward.fft(f);
	
  for(unsigned int i=0; i < n; i++) cout << f[i] << endl;
}
