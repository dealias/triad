#ifndef __Cartesian1_h__
#define __Cartesian1_h__ 1

#include "Mode.h"
#include "Basis.h"

extern int Nx;
extern int *ModeBin;

extern int Nx0,Nxb,Nxp,NPad;
extern int Nevolved;
extern unsigned int log2Nxb; // Order of FFT in each direction
extern double scale;

extern Real krmin;
extern Real krmin2;

class Cartesian1 : public Mode {
public:	
	int x;
	
	Cartesian1() {}
	Cartesian1(int column) : x(column) {}
	int Column() {return x;}
	
	Real X() const {return krmin*x;}
	Real K() const {return X();}
	Real K2() const {return krmin2*x*x;}
	
	int ModeIndex() {
		if(x > 0) return x-1; // Evolved Modes
		else return Nevolved-x-1; // Reflected Modes;
	}
};

inline Cartesian1 operator - (const Cartesian1& a)
{
	return Cartesian1(-a.x);
}

inline int operator > (const Cartesian1& a, const Cartesian1& b)
{
	return (a.x > b.x);
}

inline int operator == (const Cartesian1& a, const Cartesian1& b)
{
	return (a.x == b.x);
}

inline int operator != (const Cartesian1& a, const Cartesian1& b)
{
	return (a.x != b.x);
}

inline Cartesian1 operator + (const Cartesian1& a, const Cartesian1& b)
{
	return Cartesian1(a.x+b.x);
}

inline Cartesian1 operator - (const Cartesian1& a, const Cartesian1& b)
{
	return Cartesian1(a.x-b.x);
}

inline Cartesian1 operator * (const Cartesian1& a, const Cartesian1& b)
{
	return Cartesian1(a.x*b.x);
}

inline Cartesian1 operator / (const Cartesian1& a, const Cartesian1& b)
{
	return Cartesian1(a.x/b.x);
}

inline ixstream& operator >> (ixstream& s, Cartesian1& b)
{
	s >> b.x;
	return s;
}
	
inline oxstream& operator << (oxstream& s, const Cartesian1& b)
{
	s << b.x;
	return s;
}

inline ostream& operator << (ostream& s, const Cartesian1& b)
{
	s << b.x;
	return s;
}
	
extern Cartesian1 *Cartesian1Mode;

inline int Basis<Cartesian1>::InGrid(Cartesian1& m)
{
	return (m != Cartesian1(0)) && 
		(low.Column() <= m.Column() && m.Column() <= high.Column());
}

void set_fft_parameters();
void DiscretePad(Var *to, Var *from, Real *norm);

#if _CRAYMVP
void Cartesian1Pad(Var *to, Var *from);
void Cartesian1UnPad(Var *to, Var *from);
#else
inline void Cartesian1Pad(Var *to, const Var *from)
{
	*(to++)=0.0;
	set(to,from,Nx0);
	to += Nx0; from += Nx0;
	Var *tostop=to+NPad;
	for(; to < tostop; to++) *to=0.0;
}

inline void Cartesian1UnPad(Var *to, const Var *from)
{
	from++;
	set(to,from,Nx0);
}
#endif

// Factor which converts |y|^2 to 2.0*energy in stream function normalization:
inline Real Basis<Cartesian1>::Normalization(int) {
	return 1.0;
}

#endif
