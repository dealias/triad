#ifndef __Integrator_h__
#define __Integrator_h__ 1

#define INTEGRATOR(key) {new Entry<key,IntegratorBase>(#key,IntegratorTable);}

inline void IntegratorBase::Source(Var *src, Var *y, double t)
{
	Problem->NonLinearSrc(src,y,t);
	Problem->LinearSrc(src,y,t);
}

inline void calc_error(const Var& initial, const Var& norm, const Var& pred,
					   const Var& corr, double& errmax)
{
	Real denom,error;
	denom=max(norm2(norm),norm2(initial));
	if(denom && (pred != initial)) {
        error=norm2(corr-pred)/denom;
	    if(error > errmax) errmax=error;
	}
}

inline Solve_RC IntegratorBase::CheckError(double errmax)
{
	if(errmax > tolmax2) return UNSUCCESSFUL;
	if(errmax < tolmin2) return ADJUST;
	return SUCCESSFUL;
}

class Euler : public IntegratorBase {
public:
	char *Name() {return "Euler";}
	Solve_RC Solve(double, double);
};

class PC : public IntegratorBase {
protected:
	int new_y0;
	Var *y,*y1,*source0;
public:
	void Allocate(int n) {
		IntegratorBase::Allocate(n);
		y=y1=new Var [n]; source0=new Var [n](0); new_y0=1;
	}
	char *Name() {return "Predictor-Corrector";}
	Solve_RC Solve(double, double);
	virtual void Predictor(double, double, int, int);
	virtual int Corrector(double, double&, int, int);
	virtual void StandardPredictor(double t, double dt, int start,
								   int stop) {
		PC::Predictor(t,dt,start,stop);
	}
	virtual int StandardCorrector(double dt, double& errmax,
								  int start, int stop) {
		return PC::Corrector(dt,errmax,start,stop);
	}
};

class LeapFrog : public PC {
	Var *yp,*yp0;
	double halfdt,oldhalfdt,lasthalfdt;
public:
	void Allocate(int n) {
		PC::Allocate(n); 
		yp=new Var[n]; yp0=new Var[n]; lasthalfdt=0.0;
	}
	void TimestepDependence(double) {
		if(lasthalfdt == 0.0) for(int j=0; j < ny; j++) yp[j] = y0[j];
		halfdt=0.5*dt;
	}
	char *Name() {return "LeapFrog";}
	void Predictor(double, double, int, int);
	int Corrector(double, double&, int, int);
	void StandardPredictor(double t, double dt, int start, int stop) {
		LeapFrog::Predictor(t,dt,start,stop);
	}
	int StandardCorrector(double dt, double& errmax, int start, int stop) {
		return LeapFrog::Corrector(dt,errmax,start,stop);
	}
};

class RK2 : public PC {
protected:	
	double halfdt;
public:
	char *Name() {return "Second-Order Runge-Kutta";}
	void TimestepDependence(double);
	void Predictor(double, double, int, int);
	int Corrector(double, double&, int, int);
	void StandardPredictor(double t, double dt, int start,
						   int stop) {
		RK2::Predictor(t,dt,start,stop);
	}
	int StandardCorrector(double dt, double& errmax, int start,
						  int stop) {
		return RK2::Corrector(dt,errmax,start,stop);
	}
};

class RK4 : public PC {
protected:
	Var *source1,*source2;
	double halfdt,sixthdt;
public:
	void Allocate(int n) {
		PC::Allocate(n); source1=new Var [n](0); source2=new Var [n](0);
	}
	char *Name() {return "Fourth-Order Runge-Kutta";}
	void TimestepDependence(double);
	void Predictor(double, double, int, int);
	int Corrector(double, double&, int, int);
	void StandardPredictor(double t, double dt, int start,
						   int stop) {
		RK4::Predictor(t,dt,start,stop);
	}
	int StandardCorrector(double dt, double& errmax, int start,
						  int stop) {
		return RK4::Corrector(dt,errmax,start,stop);
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
	double pgrow, pshrink;
public:
	void Allocate(int n) {
		RK4::Allocate(n); y2=y4=y; y3=new Var [n];
		source3=source1; source1=NULL; source4=new Var [n](0);
		pgrow=0.5*0.2; pshrink=0.5*0.25;
	}
	char *Name() {return "Fifth-Order Runge-Kutta";}
	void TimestepDependence(double);
	void ExtrapolateTimestep (double errmax) {
		if(errmax < tolmin2) {if(errmax) stepfactor=pow(tolmin2/errmax,pgrow);}
		else if(tolmin2) stepinverse=stepfactor=pow(tolmin2/errmax,pshrink);
		if(errmax <= tolmax2) errmax=0.0; // Force a time step adjustment.
	}
	void Predictor(double, double, int, int);
	int Corrector(double, double&, int, int);
	void StandardPredictor(double t, double dt, int start,int stop) {
		RK5::Predictor(t,dt,start,stop);
	}
	int StandardCorrector(double dt, double& errmax, int start,int stop) {
		return RK5::Corrector(dt,errmax,start,stop);
	}
};

class Exact : public RK5 {
public:
	char *Name() {return "Exact";}
	int Microfactor(){return 100;}
};

#endif
