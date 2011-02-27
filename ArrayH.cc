#include <iostream>
#include "Complex.h"
#include "ArrayH.h"

using namespace Array;
using namespace std;

int main()
{
  unsigned int m=3;
  int mi=m;
  unsigned int n=2*m-1;
  array1H<Complex> A(m);
  for(unsigned int i=0; i < m; ++i)
    A(i)=Complex(i,10+i);

  for(int i=1-mi; i < mi; ++i)
    cout << A[i] << endl;

  /*
  int m=3;
  array2H<Complex> B(n,m);
  for(int i=0; i < n; ++i) {
    for(int j=0; j < m; ++j) {
      //B[i][j]=Complex(10+i,1+j);
      B.set(i,j,Complex(10+i,1+j));
    }
  }
  B.set(0,0,Complex(1,1));
  B.set(0,-2,Complex(0,1));

  cout << "\nB contains" << endl;
  cout << B << endl;
  cout << "but, expanded, produces" << endl;
  for(int i=0; i < n; ++i) {
    for(int j=-m+1; j < m; ++j) {
      cout << B.get(i,j) << " ";
    }
    cout << endl;
  }
  */

  return 0;
}
