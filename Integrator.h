#ifndef __Integrator_h__
#define __Integrator_h__ 1

#define INTEGRATOR(key) {new Entry<key,IntegratorBase>(#key,IntegratorTable);}

inline void IntegratorBase::Source(Var *src, Var *y, double t)
{
	if(NonlinearSrc) (*NonlinearSrc)(src,y,t);
	if(LinearSrc) (*LinearSrc)(src,y,t);
	if(ConstantSrc) (*ConstantSrc) (src,y,t);
}

inline void calc_error(const Var initial, const Var norm, const Var pred,
					   const Var corr, double& errmax)
{
	Real denom,error;
	denom=max(abs2(norm),abs2(initial));
	if(denom && (pred != initial)) {
        error=abs2(corr-pred)/denom;
	    if(error > errmax) errmax=error;
	}
}

inline Solve_RC IntegratorBase::CheckError(double errmax)
{
	if(errmax > tolmax2) return UNSUCCESSFUL;
	if(errmax < tolmin2) return ADJUST;
	return SUCCESSFUL;
}

template<class T> 
inline void set(T *to, T *from, int n) {
	memcpy(to,from,sizeof(*from)*n);
}

class Euler : public IntegratorBase {
protected:	
public:
	void Allocate(int n) {ny=n; source=new Var[n];}
	char *Name() {return "Euler";}
	Solve_RC Solve(Var *, double, double);
};

class PC : public IntegratorBase {
protected:
	Var *y,*y1,*source0;
public:
	void Allocate(int n) {ny=n; source=new Var[n];
						  y1=y=new Var [n]; source0=new Var [n];}
	char *Name() {return "Predictor-Corrector";}
	Solve_RC Solve(Var *, double, double);
	virtual void Predictor(Var *, double, double);
	virtual int Corrector(Var *, double, double&, int, int);
	virtual int StandardCorrector(Var *y0, double dt, double& errmax,
								  int start, int stop) {
		return PC::Corrector(y0,dt,errmax,start,stop);
	}
};

class RK2 : public PC {
	double halfdt;
public:
	char *Name() {return "Second-Order Runge-Kutta";}
	void TimestepDependence(double);
	void Predictor(Var *, double, double);
	int Corrector(Var *, double, double&, int, int);
	int StandardCorrector(Var *y0, double dt, double& errmax,
								  int start, int stop) {
		return RK2::Corrector(y0,dt,errmax,start,stop);
	}
};

class RK4 : public PC {
protected:
	Var *source1,*source2;
	double halfdt,sixthdt;
public:
	void Allocate(int n) {PC::Allocate(n);
						  source1=new Var [n]; source2=new Var [n];}
	char *Name() {return "Fourth-Order Runge-Kutta";}
	void TimestepDependence(double);
	void Predictor(Var *, double, double);
	int Corrector(Var *, double, double&, int, int);
	int StandardCorrector(Var *y0, double dt, double& errmax,
								  int start, int stop) {
		return RK4::Corrector(y0,dt,errmax,start,stop);
	}
};

class RK5 : public RK4 {
protected:
	Var *y2,*y3,*y4,*y5;
	Var *source3,*source4;
	double a1,a2,a3,a4,a5;
	double b10;
	double b20,b21;
	double b30,b31,b32;
	double b40,b41,b42,b43;
	double b50,b51,b52,b53,b54;
	double c0,c2,c3,c5;
	double d0,d2,d3,d4,d5;
	double pgrow, pshrink;
public:
	void Allocate(int n) {RK4::Allocate(n);
						  y2=new Var [n]; y3=new Var [n]; y4=y; y5=y;
						  source3=source1; source1=source; source4=new Var [n];
						  pgrow=0.5*0.2; pshrink=0.5*0.25;}
	char *Name() {return "Fifth-Order Runge-Kutta";}
	void TimestepDependence(double);
	void ExtrapolateTimestep (double errmax) {
		if(errmax < tolmin2) {if(errmax) stepfactor=pow(tolmin2/errmax,pgrow);}
		else if(tolmin2) stepinverse=stepfactor=pow(tolmin2/errmax,pshrink);
		if(errmax <= tolmax2) errmax=0.0; // Force a time step adjustment.
	}
	void Predictor(Var *, double, double);
	int Corrector(Var *, double, double&, int, int);
	int StandardCorrector(Var *y0, double dt, double& errmax,
								  int start, int stop) {
		return RK5::Corrector(y0,dt,errmax,start,stop);
	}
};

class Exact : public RK5 {
public:
	char *Name() {return "Exact";}
	int Microfactor(){return 100;}
};

#endif
