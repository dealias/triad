#define PI M_PI

inline double hypot(const double x, const double y)
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

inline double sgn1(const double x)
{
	return (x >= 0.0 ? 1.0 : -1.0);
}
