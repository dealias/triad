#ifndef __Integrator_h__
#define __Integrator_h__ 1

#define INTEGRATOR(key) \
{(void) new Entry<key,IntegratorBase>(#key,IntegratorTable);}

inline void IntegratorBase::CalcError(const Var& initial, const Var& norm0, 
									  const Var& pred, const Var& corr)
{
	if(pred != initial) {
		Real error=max(divide0(abs2(corr-pred),
						  max(abs2(norm0),abs2(initial))));
		if(error > errmax) errmax=error;
	}
}

inline Solve_RC IntegratorBase::CheckError()
{
	if(dynamic >= 0) {
		if(errmax > tolmax2) return UNSUCCESSFUL;
		if(errmax < tolmin2) return ADJUST;
		return SUCCESSFUL;
	}
	if(++dynamic == 0) dynamic=1;
	return SUCCESSFUL;
}

class Euler : public IntegratorBase {
public:
	const char *Name() {return "Euler";}
	Solve_RC Solve(double, double);
};

class PC : public IntegratorBase {
protected:
	int new_y0;
	Var *y,*y1,*source0;
	double halfdt;
	int order;
	double pgrow, pshrink;
public:
	PC() {order=2;}
	void Allocate(int n) {
		IntegratorBase::Allocate(n);
		y=y1=new Var [n]; source0=new(n) (Var); new_y0=1;
		pgrow=0.5/order; pshrink=0.5/(order-1);
	}
	const char *Name() {return "Predictor-Corrector";}
	Solve_RC Solve(double, double);
	
	void TimestepDependence(double dt) {
		halfdt=0.5*dt;
	}
	virtual void ExtrapolateTimestep () {
		if(errmax < tolmin2) {
			if(errmax) growfactor=pow(tolmin2/errmax,pgrow);
		} else if(tolmin2) shrinkfactor=growfactor=pow(tolmin2/errmax,pshrink);
		growfactor=min(growfactor,stepfactor);
		shrinkfactor=max(shrinkfactor,stepinverse);
		if(errmax <= tolmax2) errmax=0.0; // Force a time step adjustment.
	}
	virtual void Predictor(double, double, int, int);
	virtual int Corrector(double, int, int, int);
	virtual void StandardPredictor(double t, double dt, int start, int stop) {
		PC::Predictor(t,dt,start,stop);
	}
	virtual int StandardCorrector(double dt, int dynamic, int start,
								  int stop) {
		return PC::Corrector(dt,dynamic,start,stop);
	}
};

class Midpoint : public PC {
	const char *Name() {return "Midpoint Rule";}
	int Corrector(double, int, int, int);
};

class LeapFrog : public PC {
	Var *yp,*yp0;
	double oldhalfdt,lasthalfdt;
public:
	void Allocate(int n) {
		PC::Allocate(n); 
		yp=new Var[n]; yp0=new Var[n]; lasthalfdt=0.0;
	}
	void TimestepDependence(double dt) {
		if(lasthalfdt == 0.0) for(int j=0; j < ny; j++) yp[j] = y0[j];
		halfdt=0.5*dt;
	}
	const char *Name() {return "LeapFrog";}
	void Predictor(double, double, int, int);
	int Corrector(double, int, int, int);
	void StandardPredictor(double t, double dt, int start, int stop) { 
		LeapFrog::Predictor(t,dt,start,stop);
	}
	int StandardCorrector(double dt, int dynamic, int start, int stop) {
		return LeapFrog::Corrector(dt,dynamic,start,stop);
	}
};

class RK2 : public PC {
protected:	
public:
	const char *Name() {return "Second-Order Runge-Kutta";}
	void TimestepDependence(double dt) {
		halfdt=0.5*dt;
	}
	void Predictor(double, double, int, int);
	int Corrector(double, int, int, int);
	void StandardPredictor(double t, double dt, int start, int stop) {
		RK2::Predictor(t,dt,start,stop);
	}
	int StandardCorrector(double dt, int dynamic, int start, int stop) {
		return RK2::Corrector(dt,dynamic,start,stop);
	}
};

class RK4 : public PC {
protected:
	Var *source1,*source2;
	double sixthdt;
public:
	RK4() {order=4;}
	void Allocate(int n) {
		PC::Allocate(n); source1=new(n) (Var); source2=new(n) (Var) ;
	}
	const char *Name() {return "Fourth-Order Runge-Kutta";}
	void TimestepDependence(double);
	void Predictor(double, double, int, int);
	int Corrector(double, int, int, int);
	void StandardPredictor(double t, double dt, int start, int stop) {
		RK4::Predictor(t,dt,start,stop);
	}
	int StandardCorrector(double dt, int dynamic, int start, int stop) {
		return RK4::Corrector(dt,dynamic,start,stop);
	}
};

class RK5 : public RK4 {
protected:
	Var *y2,*y3,*y4;
	Var *source3,*source4;
	double a1,a2,a3,a4,a5;
	double b10;
	double b20,b21;
	double b30,b31,b32;
	double b40,b41,b42,b43;
	double b50,b51,b52,b53,b54;
	double c0,c2,c3,c5;
	double d0,d2,d3,d4,d5;
public:
	RK5() {order=5;}
	void Allocate(int n) {
		RK4::Allocate(n); y2=y4=y; y3=new Var [n];
		source3=source1; source1=NULL; source4=new(n) Var;
	}
	const char *Name() {return "Fifth-Order Runge-Kutta";}
	void TimestepDependence(double);
	void Predictor(double, double, int, int);
	int Corrector(double, int, int, int);
	void StandardPredictor(double t, double dt, int start,int stop) {
		RK5::Predictor(t,dt,start,stop);
	}
	int StandardCorrector(double dt, int dynamic, int start, int stop) {
		return RK5::Corrector(dt,dynamic,start,stop);
	}
};

class Exact : public RK5 {
public:
	const char *Name() {return "Exact";}
	int Microfactor(){return 100;}
};

#endif
