#include "options.h"
#include "kernel.h"
#include "Matrix.h"
#include "utils.h"

#define log2 __log2
#define pow2 __pow2

// Compute the matrix exponential of a square matrix
template<class T>
const Array2<T> exp(const Array2<T>& a)
{
  unsigned int n=a.Nx();
  assert(n == a.Ny());
  const Array2<Var> A(a);
	
  unsigned int n2=n*n;
  //	static DynVector<T> temp(n2);
  //	if(n2 > temp.Alloc()) temp.Resize(n2);
  Array2<T> B(n,n);
	
  Real Amax;
  for(unsigned int j=0; j < n2; j++) Amax=max(Amax,abs2(A(j)));
  Amax=sqrt(Amax);
	
  Real delta=1.0e-11;
  int j=max(0,1+((int) (floor(log2(Amax))+0.5)));
	
  Real scale=pow2(-j);
  B=A*scale;
  int q=1;
  Real epsilon=8.0;
  while (epsilon > delta) {
    epsilon /= 8*(4*q*q-1);
    q++;
  }
  q--;
	
  Array2<T> D(n,n);
  Array2<T> N(n,n);
  Array2<T> X(n,n);
	
  D.Identity();
  N.Identity();
  X.Identity();
	
  Real c=1.0;
  Real sign=-1.0;
	
  for(unsigned int k=1; k <= q; k++) {
    c=c*(q-k+1)/((Real) (2*q-k+1)*k);
    X=B*X;
    N += c*X;
    D += sign*c*X;
    sign *= -1.0;
  }
	
  Divide(B,D,N); // B=D^{-1} N
	
  for(unsigned int k=0; k < j; k++) B *= B;
	
  A.Purge();
  B.Hold();
  return B;
}	

int main()
{
  int n=3;
  Array2<Complex> A(n,n);
  //	cin >> A;
  A=Identity(n,n,Complex(0.0,0.0));
  cout.precision(15);
  cout << exp(A) << endl;
}

