#include "Matrix.h"

using namespace Array;

inline double abs2(double x)
{
  return x*x;
}

int main()
{
  array2<double> E(3,3);
  array2<double> I(3,3);
  array2<double> A(3,3);
  array2<double> One(3,3);
  array2<double> B(6,6);
	
  I=Identity<double>(3,3,0.0);
  One=1.0;
	
  cout << I << endl;
  cout << One << endl;
	
  for(int i=0; i < 3; i++)
    for(int j=0; j < 3; j++)
      A(i,j)=i+j;
	
  for(int i=0; i < 6; i++)
    for(int j=0; j < 6; j++)
      B(i,j)=i+j;
	
//  cin >> A;
	
  cout << A << endl;
  cout << A*A << endl;
  cout << B*B << endl;
	
  A=A*I;
  cout << A << endl;
  // E=A^{-1} I
  Divide(E,A,I);
  cout << E << endl;
  cout << E*A << endl;

  A=A*I+1.0;
  //  Or, for better optimization, use an equivalent specialized routine:
  //	MultAdd(A,A,I,One);
	
  cout << A << endl;

}
