#ifndef __Cartesian_h__
#define __Cartesian_h__ 1

#include "Mode.h"
#include "Basis.h"

extern int Nx;
extern int Ny;
extern int *ModeBin;
extern int reality;

extern int Nx0,NRows,NPad,NPadTop,xoffset;
extern int Nevolved;
extern unsigned int log2Nxb,log2Nyb; // Order of FFT in each direction
extern double scale;

extern Real krmin;
extern Real krmin2;

class Cartesian : public Mode {
public:	
	int x,y;	// wavenumber components
	
	Cartesian() {}
	Cartesian(int column, int row) : x(column), y(row) {}
	int Column() {return x;}
	int Row() {return y;}
	
	Real X() const {return krmin*x;}
	Real Y() const {return krmin*y;}
	Real K2() const {return krmin2*(x*x+y*y);}
	Real K() const {return sqrt(K2());}
	Real Th() const {
#if _CRAY		
		if(x == 0.0 && y == 0.0) return 0.0;
#endif		
		return atan2(y,x);
	}
	
	int ModeIndex() {
		if(y > 0 || (y == 0 && x > 0)) return y*Nx+x-1;
		else return Nevolved-y*Nx-x-1; // Reflected Modes;
	}
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

inline ixstream& operator >> (ixstream& s, Cartesian& b)
{
	s >> b.x >> b.y;
	return s;
}
	
inline oxstream& operator << (oxstream& s, const Cartesian& b)
{
	s << b.x << b.y;
	return s;
}

inline ostream& operator << (ostream& s, const Cartesian& b)
{
	s << "(" << b.x << ", " << b.y << ")";
	return s;
}
	
extern Cartesian *CartesianMode;

inline int Basis<Cartesian>::InGrid(Cartesian& m)
{
	return (m != Cartesian(0,0)) && 
		(low.Column() <= m.Column() && m.Column() <= high.Column()) &&
		(low.Row() <= m.Row() && m.Row() <= high.Row());
}

void set_fft_parameters();
void DiscretePad(Var *to, Var *from, Real *norm);

#if _CRAYMVP
void CartesianPad(Var *to, Var *from);
void CartesianUnPad(Var *to, Var *from);
#else
inline void CartesianPad(Var *to, const Var *from)
{
	*to=0.0;
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

// Factor which converts |y|^2 to 2.0*energy in stream function normalization:
inline Real Basis<Cartesian>::Normalization(int i) {
	return Linearity->Denominator(ModeOf(i));
}

#endif
