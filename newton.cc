#include "precision.h"
#include <iostream>

using std::cin;
using std::cout;
using std::cerr;
using std::endl;

bool newt(Real &x, Real (*f)(Real x), Real (*dfdx)(Real x), bool verbose=false,
	  unsigned int MaxIterations=100) {
  unsigned int i=0;
  if(verbose) {
    cerr.precision(16);
    cerr << x << endl;
  }

  Real lastdiff;
  Real diff=DBL_MAX;
  do {
    Real x0=x;
    
    x -= (*f)(x)/(*dfdx)(x);

    lastdiff=diff;
    
    if(verbose) cerr << x << endl;
    
    diff=fabs(x-x0);
    if(++i == MaxIterations) {
      cerr << "WARNING: Newton-Raphson iteration did not converge." << endl;
      return false;
    }
  } while(diff != 0.0 && diff < lastdiff);

  if(verbose) cerr << endl;
  return true;
}

#if 0
Real f(Real x) 
{
  return x*x-2.0;
}

Real dfdx(Real x)
{
  return 2.0*x;
}

main()
{
  Real x=1.0;
  cin >> x;
  newt(x,f,dfdx,true);
}
#endif


