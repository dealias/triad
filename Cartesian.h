#ifndef __Cartesian_h__
#define __Cartesian_h__ 1

#include "Basis.h"

class Cartesian {
public:	
	int x,y;	// wavenumber components
	
	Cartesian(int x0=0, int y0=0) : x(x0), y(y0) {}
	Real mag2() {return x*x+y*y;}
	Real K() {return sqrt(this->mag2());}
	Real Th() {return atan2(y,x);}
	Real Kx() {return x;}
	Real Ky() {return y;}
};

inline Cartesian operator - (const Cartesian& a)
{
  return Cartesian(-a.x, -a.y);
}

inline int operator > (const Cartesian& a, const Cartesian& b)
{
return (a.x > b.x || (a.x == b.x && a.y > b.y));	
}

inline int operator == (const Cartesian& a, const Cartesian& b)
{
return (a.x == b.x && a.y == b.y);	
}

inline int operator != (const Cartesian& a, const Cartesian& b)
{
return (a.x != b.x || a.y != b.y);	
}

inline Cartesian operator + (const Cartesian& a, const Cartesian& b)
{
return Cartesian(a.x+b.x, a.y+b.y);
}

inline Cartesian operator - (const Cartesian& a, const Cartesian& b)
{
return Cartesian(a.x-b.x, a.y-b.y);
}

inline Cartesian operator * (const Cartesian& a, const Cartesian& b)
{
return Cartesian(a.x*b.x, a.y*b.y);
}

inline Cartesian operator / (const Cartesian& a, const Cartesian& b)
{
return Cartesian(a.x/b.x, a.y/b.y);
}

inline ostream& operator << (ostream& os, const Cartesian& b) {
	os << "(" << b.x << ", " << b.y << ")";
	return os;
}
	
// For Navier-Stokes turbulence (velocity normalization):

template<class T>
inline Mc Basis<T>::Ckpq(T& k, T& p, T& q)
{
	return (p.Kx()*q.Ky()-p.Ky()*q.Kx())*(q.mag2()-p.mag2())/
		sqrt(k.mag2()*p.mag2()*q.mag2());
}

template<class T>
inline int Basis<T>::InGrid(Cartesian& m)
{
	return (m.x != 0 || m.y != 0) && 
		(-period.x <= m.x && m.x <= period.x) &&
		(-period.y <= m.y && m.y <= period.y);
}

#endif

