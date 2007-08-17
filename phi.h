static const long double Coeff[]={1.0,1.0/2.0,1.0/6.0,1.0/24.0,1.0/120.0,
				  1.0/720.0,1.0/5040.0,1.0/40320.0,
				  1.0/362880.0,1.0/3628800.0,1.0/39916800.0,
				  1.0/479001600.0,1.0/6227020800.0,
				  1.0/87178291200.0,1.0/1307674368000.0,
				  1.0/20922789888000.0,1.0/355687428096000.0,
				  1.0/6402373705728000.0,
				  1.0/121645100408832000.0,
				  1.0/2432902008176640000.0,
				  1.0/51090942171709440000.0,
				  1.0/1124000727777607680000.0};


// phi1(x)=(exp(x)-1)/x

#ifndef NEED_EXPM
inline double phi1(double x)
{
  return (x != 0.0) ? expm1(x)/x : 1.0;
} 
#else
inline double phi1(double x)
{
  if(fabs(x) > 1.0) return (exp(x)-1.0)/x;
  x *= 0.0625;
  register long double x2=x*x;
  register long double x3=x2*x;
  register long double x4=x2*x2;
  register long double
    sum=1+x*Coeff[1]+x2*Coeff[2]+x3*Coeff[3]+x4*Coeff[4]+x4*x*Coeff[5]
    +x4*x2*Coeff[6]+x4*x3*Coeff[7]+x4*x4*Coeff[8];
  register long double y=sum+1.0;
  register long double y2=y*y;
  register long double y4=y2*y2;
  return sum*(y+1.0)*(y2+1.0)*(y4+1.0)*(y4*y4+1.0);
}
#endif

// phi2(x)=(exp(x)-1-x)/(x*x);

inline double phi2(double x)
{
  register long double x2=x*x;
  if(fabs(x) > 1.0) return (exp(x)-x-1.0)/x2;
  register long double x3=x2*x;
  register long double x5=x2*x3;
  if(fabs(x) < 0.1) 
    return Coeff[1]+x*Coeff[2]+x2*Coeff[3]+x3*Coeff[4]+x2*x2*Coeff[5]
      +x5*Coeff[6]+x3*x3*Coeff[7]+x5*x2*Coeff[8]+x5*x3*Coeff[9];
  else {
    register long double x7=x5*x2;
    register long double x8=x7*x;
    return Coeff[1]+x*Coeff[2]+x2*Coeff[3]+x3*Coeff[4]+x2*x2*Coeff[5]
      +x5*Coeff[6]+x3*x3*Coeff[7]+x7*Coeff[8]+x8*Coeff[9]
      +x8*x*Coeff[10]+x5*x5*Coeff[11]+x8*x3*Coeff[12]+x7*x5*Coeff[13]+
      x8*x5*Coeff[14]+x7*x7*Coeff[15]+x8*x7*Coeff[16]+x8*x8*Coeff[17];
  }
}

// phi3(x)=(exp(x)-1-x-0.5*x*x)/(x*x*x)

inline double phi3(double x)
{
  register long double x2=x*x;
  register long double x3=x2*x;
  if(fabs(x) > 1.6) return (exp(x)-0.5*x2-x-1.0)/x3;
  register long double x5=x2*x3;
  if(fabs(x) < 0.1) 
    return Coeff[2]+x*Coeff[3]+x2*Coeff[4]+x3*Coeff[5]
      +x2*x2*Coeff[6]+x5*Coeff[7]+x3*x3*Coeff[8]+x5*x2*Coeff[9]
      +x5*x3*Coeff[10];
  else {
    register long double x7=x5*x2;
    register long double x8=x7*x;
    register long double x16=x8*x8;
    return Coeff[2]+x*Coeff[3]+x2*Coeff[4]+x3*Coeff[5]
      +x2*x2*Coeff[6]+x5*Coeff[7]+x3*x3*Coeff[8]+x5*x2*Coeff[9]
      +x5*x3*Coeff[10]+x8*x*Coeff[11]
      +x5*x5*Coeff[12]+x8*x3*Coeff[13]+x7*x5*Coeff[14]
      +x8*x5*Coeff[15]+x7*x7*Coeff[16]+x8*x7*Coeff[17]+x16*Coeff[18]
      +x16*x*Coeff[19]+x16*x2*Coeff[20];
  }
}
