#ifndef __Cartesian_h__
#define __Cartesian_h__ 1

#include "Basis.h"

extern int Nx;
extern int Ny;

class Cartesian {
public:	
	int x,y;	// wavenumber components
	
	Cartesian(int column=0, int row=0) : x(column), y(row) {}
	Real K2() {return x*x+y*y;}
	Real K() {return sqrt(this->K2());}
	Real Th() {return atan2(y,x);}
	Real Kx() {return x;}
	Real Ky() {return y;}
	
	int Column() {return x;}
	int Row() {return y;}
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
	
extern Cartesian *CartesianMode;

inline int Basis<Cartesian>::InGrid(Cartesian& m)
{
	return (m != Cartesian(0,0)) && 
		(low.Column() <= m.Column() && m.Column() <= high.Column()) &&
		(low.Row() <= m.Row() && m.Row() <= high.Row());
}

// For Navier-Stokes turbulence (stream function normalization):

template<class T>
inline Mc Basis<T>::Mkpq(T& k, T& p, T&)
{
	return (k.Ky()*p.Kx()-k.Kx()*p.Ky())*p.K2()/k.K2();
}

// Factor which converts |y|^2 to energy in this normalization:
inline Real Basis<Cartesian>::Normalization(int i) {return K2(i);}

#endif

