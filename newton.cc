#include <iostream>
#include <cmath>

#include "newton.h"

using namespace std;

// Root solve by Newton-Raphson
bool newton(Real &x, Real (*f)(Real x), Real (*dfdx)(Real x),
	    bool verbose=false, unsigned int MaxIterations=100)
{
  static const Real epsilon=1000.0*DBL_EPSILON;

  unsigned int i=0;
  if(verbose) {
    cerr.precision(16);
    cerr << endl << "newton0: " << x << endl;
  }

  Real lastdiff;
  Real diff=DBL_MAX;
  do {
    Real x0=x;
    
    x -= (*f)(x)/(*dfdx)(x);

    lastdiff=diff;
    
    if(verbose) cerr << endl << "newton: " << x << endl;
    
    diff=fabs(x-x0);
    if(++i == MaxIterations) {
      cerr << "WARNING: Newton-Raphson iteration did not converge." << endl;
      return false;
    }
  } while (diff != 0.0 && (diff < lastdiff || diff > epsilon*fabs(x)));

  if(verbose) cerr << endl;
  return true;
}

// Root solve by Newton-Raphson bisection

bool newton(Real &x1, Real x2, Real (*f)(Real x), Real (*dfdx)(Real x),
	    bool verbose=false, unsigned int MaxIterations=100)
{
  static const Real epsilon=1000.0*DBL_EPSILON;
  cerr.precision(16);
  
  Real f1=(*f)(x1);
  if(f1 == 0.0) return x1;
  Real f2=(*f)(x2);
  if(f2 == 0.0) return x2;
	
  if((f1 > 0.0 && f2 > 0.0) || (f1 < 0.0 && f2 < 0.0)) {
    cerr << "ERROR: root not bracketed, f(x1)=" << f1 << ", f(x2)=" << f2 
	 << endl;
    return false;
  }

  Real x=0.5*(x1+x2);
  Real dxold=fabs(x2-x1);
  if(f1 > 0.0) {
    Real temp=x1;
    x1=x2;
    x2=temp;
  }
	
  if(verbose) {
    cerr << endl << "Midpoint: " << x << endl;
  }

  Real dx=dxold;
  Real y=(*f)(x);
  Real dy=(*dfdx)(x);
  for(unsigned int j=1; j <= MaxIterations; j++) {
    if(((x-x2)*dy-y)*((x-x1)*dy-y) >= 0.0 || fabs(2.0*y) > fabs(dxold*dy)) {
      dxold=dx;
      dx=0.5*(x2-x1);
      x=x1+dx;
      if(verbose) cerr << endl << "Bisection: " << x << endl;
      if(x1 == x) return true;
    } else {
      dxold=dx;
      dx=y/dy;
      Real temp=x;
      x -= dx;
      if(verbose) cerr << endl << "Newton-Raphson: " << x << endl;
      if(temp == x) {x1=x; return true;}
    }
    if(fabs(dx) < epsilon*fabs(x)) {x1=x; return true;}
    y=(*f)(x);
    dy=(*dfdx)(x);
    if(y < 0.0) x1=x;
    else x2=x;
  }
  cerr << "WARNING: Newton-Raphson bisection did not converge." << endl;
  return false;
}


#ifdef TEST
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
  cin >> x;
  newton(x,2.0,f,dfdx,true);
  return 0;
}
#endif


