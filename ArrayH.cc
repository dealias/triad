#include <iostream>
#include "Complex.h"
#include "ArrayH.h"

using namespace Array;
using namespace std;

int main()
{
  int n=5;
  int xorigin=(n-1)/2;
  array1H<Complex> A(n);
  for(unsigned i=xorigin; i < n; ++ i)
    A[i]=Complex(i,10+i);

  for(unsigned i=0; i < n; ++ i)
    cout << A.get(i) << endl;

  int m=3;
  array2H<Complex> B(n,m);
  for(int i=0; i < n; ++ i) {
    for(int j=0; j < m; ++ j) {
      //B[i][j]=Complex(10+i,1+j);
      B.set(i,j,Complex(10+i,1+j));
    }
  }
  B.set(0,0,Complex(1,1));
  B.set(0,-2,Complex(0,1));

  cout << "\nB contains" << endl;
  cout << B << endl;
  cout << "but, expanded, produces" << endl;
  for(int i=0; i < n; ++ i) {
    for(int j=-m+1; j < m; ++ j) {
      cout << B.get(i,j) << " ";
    }
    cout << endl;
  }

  return 0;
}   
