#include <iostream>
#include "ArrayL.h"

using namespace Array;

int main()
{
  uint n=5;
	
  array2L<double> A(n);

  for(uint i=0; i < n; i++) {
    for(uint j=0; j <= i; j++) A(i,j)=10*i+j;
  }

  std::cout << A << std::endl;
	
  array2L<double> B(n);
	
  // Optimized version:
	
  for(uint i=0; i < n; i++) {
    array1<double> Bi=B[i];
    for(uint j=0; j <= i; j++) Bi[j]=10*i+j;
  }

  std::cout << B << std::endl;
   
}   
