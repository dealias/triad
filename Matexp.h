#include "utils.h"
//#include "options.h"
#include "Matrix.h"

using namespace Array;

int main()
{
  int n=3;
  Array2<Complex> A(n,n);
  	cin >> A;
//  A=Identity<Complex>(n,n);
  cout.precision(15);
  cout << exp(A) << endl;
}

