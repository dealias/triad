#ifndef __extensions_h__
#define __extensions_h__ 1

#include <math.h>

inline double hypot(double x, double y)
{
  return sqrt(x*x+y*y);
}

inline void sincos(const double x, double *sinx, double *cosx)
{
  *sinx=sin(x); *cosx=cos(x);
}

inline double sgn(const double x)
{
  return (x == 0.0 ? 0.0 : (x > 0.0 ? 1.0 : -1.0));
}

#endif

