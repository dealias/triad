// g++ -g -Wall -c newton.cc 
// g++ -g -Wall -c newtonexample.cc
// g++ -g -Wall -o newton newton.o newtonexample.o


#include <iostream>
#include <cmath>
#include "newton.h"
#include "precision.h"

using namespace std;

Real f(Real x) 
{
  return x*x*x-1.0;
}

Real dfdx(Real x)
{
  return 3.0*x*x;
}

int main()
{
  Real x=1.0;
  cout << "enter a guess:" << endl;
  cin >> x;
  Real xguess=x;
  
  // bisect:
  //Real xr=2.0;
  //newtonbisect(x,xr,f,dfdx,true,100);

  // newton-rhapson:
  newton(x,f,dfdx,true,100);

  if(x != xguess)
    cout << "you were wrong." << endl;
  else
    cout << "you totally cheated." << endl;
  return 0;
}
