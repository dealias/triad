#ifndef __Polar_h__
#define __Polar_h__ 1

#include "types.h"
#include "Bin.h"
#include "Geometry.h"

// Polar vocabulary
extern int Nr;
extern int Nth;
extern Real krmin;
extern Real krmax;
extern Real kthmin;

// Global polar variables
extern Real kthmax;
extern Real kthneg;

class Polar {
public:	
	Real r,th;	// central wavenumber components
//	int a;		// field index
	
	Polar(Real r0=0.0, Real th0=0.0) : r(r0), th(th0) {}

	Polar& operator = (const Polar& y) {
		r=y.r; th=y.th; return *this; 
	} 
};

inline int operator > (const Polar& x, const Polar& y)
{
return (x.r > y.r || (x.r == y.r && x.th > y.th));	
}

inline int operator == (const Polar& x, const Polar& y)
{
return (x.r == y.r && x.th == y.th);	
}

inline Polar operator + (const Polar& x, const Polar& y)
{
return Polar(x.r+y.r, x.th+y.th);
}

inline Polar operator - (const Polar& x, const Polar& y)
{
return Polar(x.r-y.r, x.th-y.th);
}

inline Polar operator * (const Polar& x, const Polar& y)
{
return Polar(x.r*y.r, x.th*y.th);
}

inline Polar operator / (const Polar& x, const Polar& y)
{
return Polar(x.r/y.r, x.th/y.th);
}

inline void extract(Bin<Polar> *k, double& kl, double& kg,
					double& al, double& ag)
{
	kl=k->min.r; kg=k->max.r; al=k->min.th; ag=k->max.th;
}

inline void build(Bin<Polar> *k, double kl, double kg, double al, double ag)
{
	k->min.r=kl; k->max.r=kg; k->min.th=al; k->max.th=ag;
}

inline ostream& operator << (ostream& os, const Polar& y) {
	os << "(" << y.r << ", " << y.th << ")";
	return os;
}
	
inline Real Bin<Polar>::Area()
{
	return 0.5*(max.r*max.r-min.r*min.r)*(max.th-min.th);
}

inline Real Bin<Polar>::K()
{
	return cen.r;
}

inline Real Bin<Polar>::Th()
{
	return cen.th;
}

inline Real Bin<Polar>::Kx()
{
	return cen.r*cos(cen.th);
}

inline Real Bin<Polar>::Ky()
{
	return cen.r*sin(cen.th);
}

inline int coangular(Bin<Polar> *k, Bin<Polar> *p)
{
	Real tol=10.0*REAL_EPSILON*pi;
	
	return (fabs(k->min.th - p->min.th) < tol 
			&& fabs(k->max.th - p->max.th) < tol);
}

int simpfast(Real (*)(Real), Real a, Real b, Real acc, Real& ans, Real dxmax,
			 int& iflag);
		 
typedef double POLAR_FCN(double, double, double,
						 double, double, double, double, double);

POLAR_FCN Jkpq;
		
Real BinAverage(Bin<Polar> *k, Bin<Polar> *p, Bin<Polar> *q,
				POLAR_FCN *f0, Real acc0);

// For Navier-Stokes turbulence (velocity normalization):

inline Mc Partition<Polar>::Ckpq(Polar, Polar P, Polar Q)
{
	return (Q.r-P.r)*(Q.r+P.r);
}

#endif
