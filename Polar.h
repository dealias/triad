#ifndef __Polar_h__
#define __Polar_h__ 1

#include "Partition.h"
#include "Cartesian.h"

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
	Real r,th;	// wavenumber components
//	int a;		// field index
	
	Polar(Real r0=0.0, Real th0=0.0) : r(r0), th(th0) {}

	Polar& operator = (const Polar& y) {
		r=y.r; th=y.th; return *this; 
	} 
	
	Real K() const {return r;}
	Real K2() const {return r*r;}
	Real Th() const {return th;}
	Real X() const {return r*cos(th);}
	Real Y() const {return r*sin(th);}
};

inline Real InInterval(const Cartesian& m, const Polar& a, const Polar& b)
{
	const Real xK2=m.K2();
	if(xK2 < a.K2() || xK2 >= b.K2()) return 0.0;
	
	const Real A=a.Th();
	const Real B=b.Th();
	Real opposite, xTh=m.Th();
	
	if(reality) opposite=0.0;
	else {
		if(xTh < 0) xTh += twopi;
		opposite=pi;
	}
	
	if(A == kthneg && xTh >= opposite) xTh -= twopi;
	if(B == kthmax && xTh < opposite) xTh += twopi;
	
	if(xTh < A || xTh >= B) return 0.0;
	return 1.0;
}

inline int operator == (const Polar& x, const Polar& y)
{
	return x.r == y.r && x.th == y.th;	
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

inline void extract(Bin<Polar,Cartesian> *k, double& kl, double& kg,
					double& al, double& ag)
{
	kl=k->min.r; kg=k->max.r; al=k->min.th; ag=k->max.th;
}

inline void build(Bin<Polar,Cartesian> *k, double kl, double kg, double al,
				  double ag)
{
	k->min.r=kl; k->max.r=kg; k->min.th=al; k->max.th=ag;
}

inline ostream& operator << (ostream& os, const Polar& y) {
	os << "(" << y.r << ", " << y.th << ")";
	return os;
}
	

inline Real Bin<Polar,Cartesian>::Area()
{
	if(discrete) return area;
	else return 0.5*(max.r*max.r-min.r*min.r)*(max.th-min.th);
}

inline int coangular(Bin<Polar,Cartesian> *k, Bin<Polar,Cartesian> *p)
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
		
Real BinAverage(Bin<Polar,Cartesian> *k, Bin<Polar,Cartesian> *p,
				Bin<Polar,Cartesian> *q, POLAR_FCN *f0, Real acc0);

// For Navier-Stokes turbulence (velocity normalization):

inline Mc Partition<Polar,Cartesian>::Ckpq(Polar&, Polar& P, Polar& Q)
{
	return (Q.K()-P.K())*(Q.K()+P.K());
}

// For Navier-Stokes turbulence (velocity normalization):
inline Mc Jkpq(Cartesian& k, Cartesian& p, Cartesian& q)
{
	return (k.Y()*p.X()-k.X()*p.Y())/sqrt(k.K2()*p.K2()*q.K2());
}

// Factor which converts |y|^2 to 2.0*energy in this normalization:
inline Real Partition<Polar,Cartesian>::Normalization(int) {return 1.0;}

#endif
