#include "Matrix.h"

using namespace Array;

int main()
{
  array2<double> A(3,3);
  array2<double> Ainv(3,3);
  array2<double> I(3,3);
  array2<double> One(3,3);
	
  I=Identity<double>(3,3);
  One=1.0;
	
  cout << I << endl;
  cout << One << endl;
	
  for(int i=0; i < 3; i++)
    for(int j=0; j < 3; j++)
      A(i,j)=i+j;
	
  cin >> A;
  cout << A << endl;
  
  cout << A*A << endl;
	
  A=A*I;
  cout << A << endl;
  // E=A^{-1} I
  Divide(Ainv,A,I);
  cout << Ainv << endl;
  cout << Ainv*A << endl;

  A=A*I+1.0;
  //  Or, for better optimization, use an equivalent specialized routine:
  //	MultAdd(A,A,I,One);
	
  cout << A << endl;

}
