#include "options.h"
#include "fft.h"

const double twopi=2.0*PI;

Complex expfactor[1000];

Complex gfft(Complex *data, unsigned int n, int isign, unsigned level)
{
  if(n == 1) return data[0];
  unsigned n2=n/2;
  return gfft(data,n2,isign,level+1)+
    expfactor[level]*gfft(data+(1 << level),n2,isign,level+1);
}

void fft(Complex *data, unsigned int log2n, int isign, Real scale,
	 int bitreverse) 
{
  unsigned n=1 << log2n;
  Complex out[n];
  for(unsigned k=0; k < n; k++) {
    unsigned l=0;
    for(unsigned m=n; m >= 1; m >>= 1) {
      Real phase=isign*twopi/m;
      expfactor[l++]=expi(k*phase);
    }
    out[k]=gfft(data,n,isign,0);
  }
  for(unsigned i=0; i < n; i++) data[i] = out[i] * scale;
}

int main()
{	
  unsigned int i,m,log2n;
	
  cout.precision(12);
  
  cout << "Input the base-2 logarithm of n: " << endl;
  cin >> log2n;
  cout << "Input the number of loops (1 for Demo): " << endl;
  cin >> m;
	
  unsigned int n=1 << log2n;
  Complex *f=new Complex[n];
	
  for(i=0; i < n; i++) f[i]=Complex(i,i*i);
	
  if(m == 1) { // Demo
    cout << "Data:" << endl;
	
    for(i=0; i < n; i++) cout << f[i] << endl;
    cout << endl;
	
    fft(f,log2n,1,1.0);
	
    cout << "Transform:" << endl;
	
    for(i=0; i < n; i++) cout << f[i] << endl;
    cout << endl;
	
    fft(f,log2n,-1,1.0/n);
		
    cout << "Inverse transform:" << endl;
    for(i=0; i < n; i++) cout << f[i] << endl;
		
  } else { // Timing Test
		
    Real ninv=1.0/n;
    for(i=0; i < m; i++) {
      fft(f,log2n,1,1.0);
      fft(f,log2n,-1,ninv);
    }
  }
}



