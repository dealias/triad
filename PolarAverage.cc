/* These routines calculate the geometric weight factors for a polar
   coordinate partition of the Fourier wavenumber space. The input
   arguments to the function BinAverage() are the six coordinates
   specifying the magnitude and angular limits of the three bins which are
   integrated over.  */

#include "Polar.h"

static int cnt; // Index to keep track of which variables have been swapped.

inline void SortBins(Bin<Polar,Cartesian> *& k, Bin<Polar,Cartesian> *& q,
					 int level) {
	// Don't optimize these away.
	volatile Polar Dk=k->Delta(), Dq=q->Delta();

	if(Dk.th > Dq.th || (Dk.th == Dq.th && Dk.r > Dq.r)) { 
		swap(k,q); cnt+=level;
	}
}

static Real acc; // desired relative accuracy
static Real eps,epsplus1;
static int perm; // index to our permutation of the original input parameters.
static double ComputeAverage(Bin<Polar,Cartesian> *k, Bin<Polar,Cartesian> *p,
							 Bin<Polar,Cartesian> *q);

static POLAR_FCN *f; // pointer to function $\fkpq$
static int permcnt[]= {0,3,1,2,5,4,2,1};

Real BinAverage(Bin<Polar,Cartesian> *k, Bin<Polar,Cartesian> *p,
				Bin<Polar,Cartesian> *q, POLAR_FCN *f0, Real acc0)
{
	f=f0; acc=acc0;
	
	eps=100.0*DBL_EPSILON;
	epsplus1=1.0+eps;

	cnt=0;
	SortBins(k,q,4);
	SortBins(p,q,2);

	// w.l.o.g. $\Db\ge\Da$
	if(k->Delta().th > p->Delta().th) {swap(k,p); cnt++;}

	perm=permcnt[cnt]; // index to permutation obtained from above sorts.
	return ComputeAverage(k,p,q);
}

inline void checklim(double& al, double& ag, double offset)
{
	al-=offset; ag-=offset;
	if(al > ag) msg(ERROR,"Invalid limits (%.15g > %.15g)",al,ag);
	if(ag-al > twopi) ag=al+twopi;
}

// This routine returns the branch of a in the interval [br,br+2\pi]. 
inline double branch(double a, double br)
{
	while(a < br) a+=twopi;
	while(a > br+twopi) a-=twopi;
	return a;
}

// This routine returns the principal angle of a (from the interval [0,2\pi]). 
inline double pbranch(double a)
{
	return branch(a,0.0);
}

// input parameters: $\kl,\kg,\pl,\pg,\ql,\qg,\al,\ag,\bl,\bg,\gl,\gg$
static double kl,kg,pl,pg,ql,qg,al,ag,bl,bg,gl,gg; 
static double dalpha,dbeta; // $\Da,\Db$
static double alpb,roff,ql2,qg2,dxmax3,a0,qlo,qhi;
static double glminusal[3],ggminusal[3];
int PolarAverageCount;
static int iflag;

static double pint(double);
static double kint(double);
static double rint0(double);

static double ComputeAverage(Bin<Polar,Cartesian> *k, Bin<Polar,Cartesian> *p,
							 Bin<Polar,Cartesian> *q)
{
	double a,b,ans,dxmax,offset;
	int iflag;
	const double qeps=1.0E-8;

	extract(k, kl, kg, al, ag);
	extract(p, pl, pg, bl, bg);
	extract(q, ql, qg, gl, gg);
	
	if(kg == kl || pg == pl) return(0.0);
	if(kg < kl || pg < pl) msg(ERROR,"Invalid wavenumber magnitude limits");

	offset=gl;
	checklim(al,ag,offset); checklim(bl,bg,offset); checklim(gl,gg,offset);
	
	dalpha=ag-al;
	dbeta=bg-bl;
	alpb=pbranch(al);

	glminusal[0]=gl-alpb;
	ggminusal[0]=gg-alpb;
	for(int n=1; n <= 2; n++) {
		glminusal[n]=glminusal[0]+n*twopi;
		ggminusal[n]=ggminusal[0]+n*twopi;
	}

	qlo=ql*(1.0-qeps);
	qhi=qg*(1.0+qeps);

	roff=pbranch(bl-al);
	a0=al; // Used in rint0

	dxmax3=pi; // dxmax for rint0

	ql2=ql*ql;
	qg2=qg*qg;
	PolarAverageCount=0;

	a=pl;
	b=pg;
	dxmax=0.4*(qg-ql);

	// Call adaptive Simpson integration routine.
	// See documentation in simp for a description of possible errors.
	if(!simpfast(pint,a,b,acc,ans,dxmax,iflag))
		msg(OVERRIDE,"First nested call to simp returned code %d",iflag);
	return ans;
}

static double p,p2,ql2mp2,qg2mp2;

static double pint(double p0)
{
	double a,b,ans,dxmax;

	p=p0;
	p2=p*p;
	ql2mp2=ql2-p2; // Used later in solving \Eq(rrestrict)
	qg2mp2=qg2-p2;

	a=kl;
	b=kg;
	dxmax=0.4*(qg-ql);

	if(!simpfast(kint,a,b,acc,ans,dxmax,iflag))
		msg(OVERRIDE,"Second nested call to simp returned code %d",iflag);
	return p*ans;
}

static double k,k2,k2pp2,twokp,term,sum;
static int coeffa,coeffb;

inline void rbarterm(double a, double b, int coeffa0, int coeffb0,
					 double term0)
{
	double ans;
	coeffa=coeffa0; coeffb=coeffb0; term=term0;
	
	if(b > a*(1.0+sgn1(a)*eps)) {
		if(!simpfast(rint0,a,b,acc,ans,dxmax3,iflag))
			msg(OVERRIDE,"Third nested call to simp returned code %d",iflag);
		sum+=ans;
	}
}

void rbarint(double r1, double r2)
{
	int j,loopcount=1;
	double r1true;
	
	if(r2 < 0.0) {r1+=twopi; r2+=twopi;}
	r1true=r1;

	if(r1 < 0.0) {r1=0.0; loopcount=2;}
	// split the interval into two subintervals and carry out the following
	// procedure for each subinterval

	for(j=1; j <= loopcount; ++j) { // evaluate \Eq(rbarint)
		// $\abgr=\Db-\rb$, $\ablr=0$
		rbarterm(max(dbeta-dalpha,r1),min(dbeta,r2),0,-1,dbeta);

		// $\abgr=\Da$, $\ablr=0$
		rbarterm(max(0.0,r1),min(dbeta-dalpha,r2),0,0,dalpha);

		// $\abgr=\Da, $\ablr=-r$
		rbarterm(max(-dalpha,r1-twopi),min(0.0,r2-twopi),-1,0,dalpha);

		r1=twopi+r1true; // Setup for second sub-interval (in case of
		r2=twopi; // a split $[\rb_1,\rb_2]$ interval)
	}
}

static double kint(double k0)
{
	double cosr1,cosr2,r1=0.0,r2=0.0,r3,r4;

	k=k0;
	k2=k*k;
	k2pp2=k2+p2;
	twokp=2.0*k*p;
	if(twokp != 0.0) {
		// solve \Eq(rrestrict) in the interval $[0,\pi]$
		cosr1=(qg2mp2-k2)/twokp; 
		cosr2=(ql2mp2-k2)/twokp;
		if(fabs(cosr1) <= 1.0) r1=acos(cosr1);
		else r1=(cosr1 > 0.0) ? 0.0 : pi; 
		// if cosr1 is out of range then $r_1$ is set to $0$ or $\pi$
		if(fabs(cosr2) <= 1.0) r2=acos(cosr2);
		else r2=(cosr2 > 0.0) ? 0.0 : pi;
	} // Note $0\le r_1 \le r_2 \le \pi$
	
	if(p == 0.0) {
		r1=0.0;
		if(ql <= k && k < qg) r2=pi; // Special case of magnitude restriction
		else r2=r1; // magnitude restriction is violated
	}
	
	if(k==0.0) {
		r1=0.0;
		if(ql <= p && p < qg) r2=pi;
		else r2=r1;
	}

	r3=twopi-r2; r4=twopi-r1; // $r_{3,4}$ are from interval $[\pi,2\pi]$
	r1-=roff; r2-=roff; // Subtract $r_{off}$ to obtain $\rb$ values
	r3-=roff; r4-=roff;
	sum=0.0;
	
	rbarint(r1,r2); // integrate over $[\rb_1,\rb_2]$ interval
	rbarint(r3,r4); // integrate over $[\rb_3,\rb_4]$ interval

	return k*sum;
}

static double rint0(double rbar)
{
	double q,theta,x,y,sum;
	double alphabl,alphabg,r,cosr,sinr,b0,g0,lo,hi;
	int n;
	//	PolarAverageCount++; // total number of function evaluations
	r=rbar+roff; // $r=\rb+r_{off}$
	alphabg=term+coeffb*rbar; // $\abgr=term+coeffb*\rb$
	alphabl=coeffa*rbar; // $\ablr=term+coeffb*\rb$

	sincos(r,&sinr,&cosr);

	q=sqrt(k2pp2+twokp*cosr); // $q=\abs{\vk+\vp}$
	if(q < qlo || q >= qhi) msg(ERROR,"Magnitude of q is out of bounds");

	y=-p*sinr; x=-k-p*cosr;
	if(x == 0.0) theta=0.0; // atan2(y,0.0) is undefined.
	else theta=atan2(y,x); // returns angle in $[-\pi,\pi]$
	if(theta < 0.0) theta+=twopi; // we want $\th\in[0,2\pi]$

	b0=rbar+bl; g0=theta+al;
	// $\a=\ab+\a_0, \b=\ab+b_0, \g=\ab+g_0$

	sum=0.0;
	for(n=0; n <= 2; n++) { // n indexes offsets of $0$, $2\pi$, $4\pi$
		hi=min(alphabg,ggminusal[n]-theta);
		lo=max(alphabl,glminusal[n]-theta);

		// \Eq(angular)
		if(hi > lo*epsplus1) switch(perm) {
		case 0:
			sum+=(*f)(k,p,q,a0,b0,g0,lo,hi); break;
		case 1:
			sum+=(*f)(k,q,p,a0,g0,b0,lo,hi); break;
		case 2:
			sum+=(*f)(p,q,k,b0,g0,a0,lo,hi); break;
		case 3:
			sum+=(*f)(p,k,q,b0,a0,g0,lo,hi); break;
		case 4:
			sum+=(*f)(q,k,p,g0,a0,b0,lo,hi); break;
		case 5:
			sum+=(*f)(q,p,k,g0,b0,a0,lo,hi); break;
		}
	}
	return sum;
}
