#include <iostream>
#include "DynVector.h"

//using namespace Array;

int main()
{
  uint n=5;
	
  DynVector<int> A(n);

  for(uint i=0; i < n; i++) {
    A[i] = i+1;
  }

  std::cout << A << std::endl;

  A.Push(-1);

  std::cout << A << std::endl;
  
  A.sort();

  std::cout << A << std::endl;
}   
