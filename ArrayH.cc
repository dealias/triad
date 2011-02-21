#include <iostream>
#include "Complex.h"
#include "Arrayh.h"
#include "ArrayL.h"

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
  for(unsigned i=0; i < n; ++ i) {
    for(unsigned j=0; j < m; ++ j) {
      B[i][j]=Complex(10+i,1+j);
    }
  }
  cout << "\nB contains only" << endl;
  cout << B << endl;
  cout << "\nbut, expanded, produces" << endl;
  for(int i=0; i < n; ++ i) {
    for(int j=-m+1; j < m; ++ j) {
      cout << B.get(i,j) << " ";
    }
    cout << endl;
  }

  return 0;
}   
