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

inline Real InInterval(const Cartesian& x, const Polar& a, const Polar& b)
{
	const Real width=10000.0;
	const Real fuzz=width*REAL_EPSILON;
	
	const Real K2fuzz=(b.K2()-a.K2())*fuzz;
	const Real a1=a.K2()-K2fuzz, a2=a.K2()+K2fuzz;
	const Real b1=b.K2()-K2fuzz, b2=b.K2()+K2fuzz;
	
	const Real xK2=x.K2();
	if(xK2 < a1 || xK2 >= b2) return 0.0;
	
	const Real Thfuzz=(b.Th()-a.Th())*fuzz;
	const Real A1=a.Th()-Thfuzz, A2=a.Th()+Thfuzz;
	const Real B1=b.Th()-Thfuzz, B2=b.Th()+Thfuzz;
	Real opposite, xTh=x.Th();
	
	if(reality) {
		opposite=0.0;
	} else {
		if(xTh < 0) xTh += twopi;
		opposite=pi;
	}
	
	if(a.Th() == kthneg && xTh >= opposite) xTh -= twopi;
	if(b.Th() == kthmax && xTh < opposite) xTh += twopi;
	
	if(xTh < A1 || xTh >= B2) return 0.0;
	if(xK2 >= a2 && xK2 < b1 && xTh >= A2 && xTh < B1) return 1.0;
	
	Real factor=1.0;
	if(xK2 < a2 && a.K() > krmin) factor *= (xK2-a1)/(a2-a1);
	else if(xK2 >= b1 && b.K() < krmax) factor *= (xK2-b1)/(b2-b1);
	if(xTh < A2) factor *= (xTh-A1)/(A2-A1);
	if(xTh >= B1) factor *= (xTh-B1)/(B2-B1);
	
	return factor;
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
