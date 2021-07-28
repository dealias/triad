/* Robust efficient routines to compute phi_n(x) for n=1 to 4.
   Copyright (C) 2005-2021
   John C. Bowman, University of Alberta

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. */

#ifndef __phi_h__
#define __phi_h__ 1

#define __PHI_H_VERSION__ 1.02

#include <cmath>

static const long double Coeff[]={1.0L,1.0L/2.0L,1.0L/6.0L,1.0L/24.0L,
                                  1.0L/120.0L,1.0L/720.0L,1.0L/5040.0L,
                                  1.0L/40320.0L,1.0L/362880.0L,
                                  1.0L/3628800.0L,1.0L/39916800.0L,
				  1.0L/479001600.0L,1.0L/6227020800.0L,
				  1.0L/87178291200.0L,1.0L/1307674368000.0L,
				  1.0L/20922789888000.0L,
                                  1.0L/355687428096000.0L,
				  1.0L/6402373705728000.0L,
				  1.0L/121645100408832000.0L,
				  1.0L/2432902008176640000.0L,
				  1.0L/51090942171709440000.0L,
				  1.0L/1124000727777607680000.0L};

// phi1(x)=(exp(x)-1)/x
#ifndef NEED_EXPM1
inline double phi1(double x)
{
  return (x != 0.0) ? expm1(x)/x : 1.0;
}
#else
// Use the identity phi1(2x)=(0.5*x*phi1(x)+1)*phi1(x)
inline double phi1(double x)
{
  if(fabs(x) > 0.78) return (exp(x)-1.0)/x;
  x *= 0.125;
  long double x2=x*x;
  long double x3=x2*x;
  long double x4=x2*x2;
  long double y=1+x*Coeff[1]+x2*Coeff[2]+x3*Coeff[3]+x4*Coeff[4]+
    x4*x*Coeff[5]+x4*x2*Coeff[6]+x4*x3*Coeff[7]+x4*x4*Coeff[8];
  y *= 0.5*x*y+1.0;
  y *= x*y+1.0;
  return y*(2.0*x*y+1.0);
}
#endif

// phi2(x)=(exp(x)-1-x)/(x^2);
// Use the identity phi2(2x)=0.25*(x*phi2(x)+1)^2+0.5*phi2(x);
inline double phi2(double x)
{
  if(fabs(x) > 1.0) return (exp(x)-x-1.0)/(x*x);
  x *= 0.125;
  long double x2=x*x;
  long double x3=x2*x;
  long double x5=x2*x3;
  long double y=Coeff[1]+x*Coeff[2]+x2*Coeff[3]+x3*Coeff[4]+x2*x2*Coeff[5]
      +x5*Coeff[6]+x3*x3*Coeff[7]+x5*x2*Coeff[8]+x5*x3*Coeff[9];
  long double y1=x*y+1.0;
  y=0.25*y1*y1+0.5*y;
  y1=x*y+0.5;
  y=y1*y1+0.5*y;
  y1=2.0*x*y+0.5;
  return y1*y1+0.5*y;
}

// phi3(x)=(exp(x)-1-x-x^2/2)/(x^3)
// Use the identity phi3(2x)=0.125*phi2(x)*(x*phi2(x)+2)+0.25*phi3(x)
// where phi2(x)=x*phi3(x)+0.5
inline double phi3(double x)
{
  if(fabs(x) > 1.6) return (exp(x)-0.5*x*x-x-1.0)/(x*x*x);
  x *= 0.125;
  long double x2=x*x;
  long double x3=x2*x;
  long double x5=x2*x3;
  long double y=Coeff[2]+x*Coeff[3]+x2*Coeff[4]+x3*Coeff[5]
    +x2*x2*Coeff[6]+x5*Coeff[7]+x3*x3*Coeff[8]+x5*x2*Coeff[9]
    +x5*x3*Coeff[10];
  long double y2=x*y+0.5;
  y=0.125*y2*(x*y2+2)+0.25*y;
  y2=2.0*x*y+0.5;
  y=0.25*y2*(x*y2+1.0)+0.25*y;
  y2=4.0*x*y+0.5;
  return 0.25*y2*(2.0*x*y2+1.0)+0.25*y;
}

// phi4(x)=(exp(x)-1-x-x^2/2-x^3/6)/(x^4)
// Use the identity phi4(2x)=0.0625*(x*phi3(x)+0.5)^2+0.125*(phi3(x)+phi4(x));
// where phi3(x)=x*phi4(x)+1/6
inline double phi4(double x)
{
  if(fabs(x) > 1.6) return (exp(x)-Coeff[2]*x*x*x-0.5*x*x-x-1.0)/(x*x*x*x);
  x *= 0.125;
  long double x2=x*x;
  long double x3=x2*x;
  long double x4=x2*x2;
  long double x5=x2*x3;
  long double y=Coeff[3]+x*Coeff[4]+x2*Coeff[5]+x3*Coeff[6]
    +x4*Coeff[7]+x5*Coeff[8]+x3*x3*Coeff[9]+x5*x2*Coeff[10]
    +x4*x4*Coeff[11];
  long double y3=x*y+Coeff[2];
  long double y2=x*y3+0.5;
  y=0.0625*y2*y2+0.125*(y3+y);
  y3=2.0*x*y+Coeff[2];
  y2=0.5*x*y3+0.125;
  y=y2*y2+0.125*(y3+y);
  y3=4.0*x*y+Coeff[2];
  y2=x*y3+0.125;
  return y2*y2+0.125*(y3+y);
}

#endif
