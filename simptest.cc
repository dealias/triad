#include <iostream>
#include <cmath>

#include "precision.h"

using std::cout;
using std::endl;

int simpfast(Real (*)(Real), Real a, Real b, Real acc, Real& sum, Real dxmax, int& iflag);

main()
{
double sum,error;
double f(double x);
int iflag;

cout.precision(12);
simpfast(f,0.0,1.0,1e-10,sum,1,iflag);
cout << "sum=" << sum << ", iflag=" << iflag << endl;
}

double n=0.3;
double K=0.4;

double f(double x)
{
//  double sx=sin(x);
//  return 1.0/(sqrt(1.0-K*K*sx*sx));
//  return 1.0/((1+n*sx*sx)*sqrt(1.0-K*K*sx*sx));
//  return 4.0*sqrt(1.0-x*x);
  if(x == 0) return 0.0;
  return 2.0/sqrt(1.0-x*x);
}
