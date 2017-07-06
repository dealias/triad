#include "utils.h"
#include "Matrix.h"

using namespace Array;

int main()
{
  int n=3;
  Array2<Complex> A(n,n);
//  A=Identity<Complex>(n,n);
  	cin >> A;
  cout.precision(15);
  cout << exp(A) << endl;
}

