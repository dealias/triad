#include <iostream>
#include "Complex.h"
using namespace std;
#include "ArrayH.h"

using namespace Array;

int main()
{
  unsigned int m=3;
  int mi=m;
  unsigned int n=2*m-1;
  Array1H<Complex> A(m);
  
  for(unsigned int i=0; i < m; ++i)
    A.set(i,Complex(10+i,i));
  
//  for(int i=1-mi; i <= 0; ++i)
//    A.set(i,Complex(10-i,i));

  for(int i=1-mi; i < mi; ++i)
    cout << A(i) << endl;

  Array2H<Complex> B(mi,mi);
  for(int i=-mi+1; i < mi; ++i) {
    for(int j=0; j < mi; ++j) {
      B.set(i,j,Complex(10+i,j));
    }
  }
  
  cout << endl;
  for(int i=-mi+1; i < mi; ++i) {
    for(int j=-mi+1; j < mi; ++j) {
      cout << B(i,j) << " ";
    }
    cout << endl;
  }

  cout << endl;
  for(int i=-mi+1; i < mi; ++i) {
    Array1HH<Complex> Bi=B[i];
    for(int j=-mi+1; j < mi; ++j) {
      cout << Bi(j) << " ";
    }
    cout << endl;
  }

  return 0;
}
