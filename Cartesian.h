#ifndef __Cartesian_h__
#define __Cartesian_h__ 1

#include "Basis.h"

extern int Nx;
extern int Ny;

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
		Real th=atan2(y,x);
		if(!reality && th < 0) return th+twopi;
		if(x < 0 && y == 0.0) th=-PI;
		return th;
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

extern int NRows,NPad;
extern int *RowBoundary;
extern Var *ZeroBuffer; 

#if _CRAY
int CartesianPad(Var *to, Var *from);
void CartesianUnPad(Var *to, Var *from);
#else
inline int CartesianPad(Var *to, const Var *from)
{
	for(int i=0; i < NRows; i++) {
		int ncol=RowBoundary[i+1]-RowBoundary[i];
		set(to,from,ncol);
		to += ncol; from += ncol;
		set(to,ZeroBuffer,NPad);
		to += NPad;
	}
	return to-psibuffer;
}

inline void CartesianUnPad(Var *to, const Var *from)
{
	for(int i=0; i < NRows; i++) {
		int ncol=RowBoundary[i+1]-RowBoundary[i];
		set(to,from,ncol);
		to += ncol; from += ncol+NPad;
	}
}
#endif

// For Navier-Stokes turbulence (stream function normalization):
template<class T>
inline Mc Mkpq(T& k, T& p, T&)
{
	return (k.Y()*p.X()-k.X()*p.Y())*p.K2()/k.K2();
}

// Factor which converts |y|^2 to energy in this normalization:
inline Real Basis<Cartesian>::Normalization(int i) {return K2(i);}

#endif

