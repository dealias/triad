#ifndef __extensions_h__
#define __extensions_h__ 1

#include <math.h>

#if MATH_EXTENSION
inline double hypot(double x, double y)
{
  return sqrt(x*x+y*y);
}

inline double log2(double x)
{
  static const double Log2=log(2);
  return log(x)/Log2;
}

inline double pow2(double x)
{
  return pow(2,x);
}

inline void sincos(const double x, double *sinx, double *cosx)
{
  *sinx=sin(x); *cosx=cos(x);
}
#endif

inline double sgn(const double x)
{
  return (x == 0.0 ? 0.0 : (x > 0.0 ? 1.0 : -1.0));
}

#ifdef HAVE_LONGEXP
inline long double exp(long double x)
{
  return expl(x);
}

inline long double expratio1(long double x)
{
  return (x != 0.0) ? expm1l(x)/x : 1.0;
} 

inline long double expratio1l(long double x)
{
  if(fabsl(x) > 1.0) return (expl(x)-1.0)/x;
  x *= 0.0625;
  register long double x2=x*x;
  register long double x3=x2*x;
  register long double x4=x2*x2;
  register long double x8=x4*x4;
  register long double
    sum=1+x*Coeff[1]+x2*Coeff[2]+x3*Coeff[3]+x4*Coeff[4]+x4*x*Coeff[5]
    +x4*x2*Coeff[6]+x4*x3*Coeff[7]+x8*Coeff[8]+x8*x*Coeff[9];
  register long double y=sum+1.0;
  register long double y2=y*y;
  register long double y4=y2*y2;
  return sum*(y+1.0)*(y2+1.0)*(y4+1.0)*(y4*y4+1.0);
}
#endif

#endif

