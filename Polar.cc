#include "kernel.h"
#include "Polar.h"

// Polar vocabulary
int Nr=3;
int Nth=6;
Real krmin=1.0;
Real krmax=2.0;
Real kthmin=0.0;

// Global polar variables
Real kthmax;
Real kthneg;

char *Partition<Polar,Cartesian>::Name() {return "Polar";}

char *Partition<Polar,Cartesian>::WeightFileName(char *suffix) {
	char *filename=new char[strlen(Problem->Abbrev())+strlen(Name())+
	strlen(suffix)+100];
	sprintf(filename,"%s/%s/%dx%d_%g:%g_%g%s%s%s",Problem->Abbrev(),
			downcase(Name()),Nr,Nth,krmin,krmax,kthmin,reality ? "R" : "",
			discrete ? "D" : "",suffix);
	return filename;
}

void angular_grid(Bin<Polar,Cartesian> *grid, Real kthmin, Real kthmax,
				  Real delta, int start, int stop)
{
	Real kthbnd=kthmin;
	Real kth=kthmin+0.5*delta;
	int j;
	
	for(j=start; j < stop; j++) {
		grid[j].min.th=kthbnd;
		grid[j].cen.th=kth;
		kthbnd+=delta; kth+=delta;
	}
	for(j=start; j < stop-1; j++) grid[j].max.th=grid[j+1].min.th;
	grid[stop-1].max.th=kthmax;
}	

void radial_grid(Bin<Polar,Cartesian> *grid, Real krmin, Real krmax, Real
				 delta, int start, int stop)
{
	Real krbnd=krmin;
	Real kr=krmin*sqrt(delta);
	int j;
	
	for(j=start; j < stop; j++) {
		grid[j].min.r=krbnd;
		grid[j].cen.r=kr;
		krbnd*=delta; kr*=delta;
	}
	for(j=start; j < stop-1; j++) grid[j].max.r=grid[j+1].min.r;
	grid[stop-1].max.r=krmax;
}	

void Partition<Polar,Cartesian>::MakeBins()
	// Radial: logarithmic, Angular: uniform
{
	Bin<Polar,Cartesian> *grid,*p;
	Real delta;
	int i,j,Nthpos;
	
	grid=new Bin<Polar,Cartesian> [max(Nr,Nth)];
	
	if(reality) {
		Nthpos=Nth/2;
		if(2*Nthpos != Nth)
			msg(ERROR,"The reality condition requires that Nth be even");
		kthneg=kthmin-pi;
		kthmax=kthmin+pi;
	} else {
		Nthpos=Nth;
		kthneg=kthmin;
		kthmax=kthmin+twopi;
	}
	
	delta=twopi/Nth;

	// primary grid
	angular_grid(grid,kthmin,kthmax,delta,0,Nthpos);
	
	// reflected grid
	if(reality) angular_grid(grid,kthneg,kthmin,delta,Nthpos,Nth);
	
	delta=pow(krmax/krmin,1.0/Nr);	
	radial_grid(grid,krmin,krmax,delta,0,Nr);
	
	n=Nr*Nth;
	p=bin=new Bin<Polar,Cartesian>[n];
	Nmode=Nr*Nthpos;
	nindependent=(reality || Nth % 2) ? Nmode : n/2;
	
	for(j=0; j < Nth; j++) {
		for(i=0; i < Nr; i++) {
			p->min.r=grid[i].min.r;
			p->cen.r=grid[i].cen.r;
			p->max.r=grid[i].max.r;
			p->min.th=grid[j].min.th;
			p->cen.th=grid[j].cen.th;
			p->max.th=grid[j].max.th;
			p++;
		}
	}
	
	if(p-bin != n) 
		msg(ERROR,"Calculated number and actual number of bins disagree."); 

	if(discrete) for (i=0; i < n; i++) bin[i].MakeModes();
	return;
}

void Bin<Polar,Cartesian>::MakeModes()
{
	int cx=(int) cen.X(), cy=(int) cen.Y();
	int w=1, last_nmode=-1;
	Count(Cartesian(cx,cy));
	while(nmode != last_nmode) {
		last_nmode=nmode;
		for(int i=cx-w; i <= cx+w; i++) {
			Count(Cartesian(i,cy+w));
			Count(Cartesian(i,cy-w));
		}
		for(int j=cy-w+1; j < cy+w; j++) {
			Count(Cartesian(cx+w,j));
			Count(Cartesian(cx-w,j));
		}
		w++;
	}
	for(int i=0; i < nmode; i++) cout << mode[i] << endl;
	cout << endl;
}

static const Real linacc=0.01;
static const Real dxmax=REAL_MAX;

Real frequency(const Polar& v);
Real growth(const Polar& v);

static Real k0;
static Bin<Polar,Cartesian> b;

static Real growthk(Real th)
{
	return growth(Polar(k0,th));
}

static Real frequencyk(Real th)
{
	return frequency(Polar(k0,th));
}

static Real ThAveragedGrowth(Real k) {
	Real ans;
	int iflag;
	k0=k;
	if(!simpfast(growthk,b.min.th,b.max.th,linacc,ans,pi,iflag))
		msg(ERROR,"Simp returned code %d",iflag);
	return k*ans;
}

static Real ThAveragedFrequency(Real k) {
	Real ans;
	int iflag;
	k0=k;
	if(!simpfast(frequencyk,b.min.th,b.max.th,linacc,ans,pi,iflag))
		msg(ERROR,"Simp returned code %d",iflag);
	return k*ans;
}

void BinAveragedLinearity(Real& nu)
{
	Real ans;
	int iflag;
	if(!simpfast(ThAveragedGrowth,b.min.r,b.max.r,linacc,
				 ans,dxmax,iflag)) msg(ERROR,"Simp returned code %d",iflag);
	nu=-ans;
}

void BinAveragedLinearity(Complex& nu)
{
	Real ans;
	int iflag;
	BinAveragedLinearity(ans);
	nu=ans;
	if(!simpfast(ThAveragedFrequency,b.min.r,b.max.r,linacc,
				 ans,dxmax,iflag)) msg(ERROR,"Simp returned code %d",iflag);
	nu+=I*ans;
}

Nu Partition<Polar,Cartesian>::Linearity(int i)
{
	Nu nu;
	b=bin[i];
	BinAveragedLinearity(nu);
	return nu/Area(i);
}

Mc Partition<Polar,Cartesian>::
ComputeBinAverage(Bin<Polar,Cartesian> *k, Bin<Polar,Cartesian> *p,
				  Bin<Polar,Cartesian> *q)
{
	Real sum;
	if(discrete) {
		sum=0.0;
		Cartesian *mk=k->mode.Base(), *mkstop=mk+k->nmode; 
		Cartesian *mp=p->mode.Base(), *mpstop=mp+p->nmode; 
		for(; mk < mkstop; mk++) {
			for(; mp < mpstop; mp++) {
				Cartesian mq=-*mk-*mp;
				if(q->InBin(mq)) sum += Mkpq(*mk,*mp,mq);
			}
		}
	} else {
		Real acc=1.0E-2;
		sum=BinAverage(k,p,q,Jkpq,acc);
	}
	return sum;
}

// For Navier-Stokes turbulence (velocity normalization):

Mc Jkpq(double k, double, double, double, double b, double g,
		double lo, double hi)
{
	return (hi-lo)*sin(g-b)/k;
}
