#ifndef __Cartesian_h__
#define __Cartesian_h__ 1

#include "Basis.h"

extern int Nx;
extern int Ny;

extern int xoffset;
extern int Nx0,NRows,NPad,NPadTop;
extern unsigned int log2Nxb,log2Nyb; // Order of FFT in each direction
extern double scale;

class Cartesian {
public:	
	int x,y;	// wavenumber components
	
	Cartesian() {}
	Cartesian(int column, int row) : x(column), y(row) {}
	Real K2() const {return x*x+y*y;}
	Real K() const {return sqrt(K2());}
	Real Th() const {
#if _CRAY		
		if(x == 0.0 && y == 0.0) return 0.0;
#endif		
		return atan2(y,x);
	}
	Real X() const {return x;}
	Real Y() const {return y;}
	
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

#if _CRAY || 1
void CartesianPad(Var *to, Var *from);
void CartesianUnPad(Var *to, Var *from);
#else
inline void CartesianPad(Var *to, const Var *from)
{
	to += xoffset;
	*(to++)=0.0;
	set(to,from,Nx0);
	to += Nx0; from += Nx0;
	for(int j=0; j < NRows-1; j++) {
		Var *tostop=to+NPad;
		for(; to < tostop; to++) *to=0.0;
		set(to,from,Nx);
		to += Nx; from += Nx;
	}
	Var *tostop=to+NPadTop;
	for(; to < tostop; to++) *to=0.0;
}

inline void CartesianUnPad(Var *to, const Var *from)
{
	from += xoffset+1;
	set(to,from,Nx0);
	to += Nx0; from += Nx0+NPad;
	for(int j=1; j < NRows; j++) {
		set(to,from,Nx);
		to += Nx; from += Nx+NPad;
	}
}
#endif

// For Navier-Stokes turbulence (stream function normalization):
template<class T>
inline Mc Mkpq(T& k, T& p, T&)
{
	return (k.Y()*p.X()-k.X()*p.Y())*p.K2()/k.K2();
}

// Factor which converts |y|^2 to 2.0*energy in this normalization:
inline Real Basis<Cartesian>::Normalization(int i) {return K2(i);}

#endif

