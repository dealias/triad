#ifndef __PolarBin_h__
#define __PolarBin_h__ 1

#include "Partition.h"
#include "Polar.h"
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

inline int InInterval(const Cartesian& m, const Polar& a, const Polar& b)
{
	const Real xK2=m.K2();
	if(xK2 < a.K2() || xK2 >= b.K2()) return 0;
	
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
	
	if(xTh < A || xTh >= B) return 0;
	return 1;
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

inline Real Bin<Polar,Cartesian>::Area()
{
	if(discrete) return (Real) nmode;
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