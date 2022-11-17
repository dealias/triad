#include <iostream>
#include "Complex.h"
using namespace std;
#include "array2h.h"

using namespace Array;

int main()
{
  int mx=4;
  int my=3;

  array2h<Complex> A(mx,my);
  for(int i=-mx+1; i < mx; ++i) {
    array1<Complex> Ai=A[i];
    for(int j=i <= 0 ? 1 : 0; j < my; ++j) {
      Ai[j]=Complex(i,j);
    }
  }

  cout << endl;

  for(int i=-mx+1; i < mx; ++i) {
    for(int j=i <= 0 ? 1 : 0; j < my; ++j) {
      cout << A(i,j) << " ";
    }
    cout << endl;
  }

  cout << endl;
  cout << A << endl;

  return 0;
}
