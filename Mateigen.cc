#include "precision.h"
#include "Complex.h"
#include "Matrix.h"

extern "C" void zgees_(const char& jobvs, const char& sort, void (*select)(),
		       const int& n, Complex *a, const int& lda, int *sdim,Complex *w,
		       Complex *vs, const int& ldvs, Complex *work, const int& lwork, 
		       Real *rwork, char *bwork, int *info);

namespace Array {
  
// Compute the eigenvalues of a square complex matrix
const Array1<Complex> eigen(const Array2<Complex>& A)
{
  unsigned int n=A.Nx();
  assert(n == A.Ny());
	
  Array2<Complex> B(n,n);
  Array1<Complex> E(n);
	
  char jobvs='N';
  char sort='N';
  void (*select)()=NULL;
  B=A;
  int lda=n;
  int sdim;
  Complex *w=E;
  Complex *vs=NULL;
  int ldvs=1;
  int lwork=(n ? 4*n : 1);
  Real *rwork=new Real[n];
  Complex work0[1];
  char *bwork=NULL;
  int info;
  Complex *work=new Complex[lwork];
	
  zgees_(jobvs,sort,select,n,B,lda,&sdim,w,vs,ldvs,work,lwork,rwork,bwork,
	 &info);
	
  A.Purge();
  E.Hold();
  return E;
}	

}

using namespace Array;

#if 1
int main()
{
  int n=3;
  Array2<Complex> A(n,n);
  cin >> A;
  cout << endl;
  cout << A << endl;
  cout << endl;
  cout << eigen(A) << endl;
}
#endif
