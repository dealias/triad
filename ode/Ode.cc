#include "options.h"
#include "kernel.h"

char *method="Ode";
char *integrator="PC";

// Global variables
Nu *nu,*nu_inv;

// Local vocabulary declarations and default values
static Nu nu0=1.0;

class OdeVocabulary : public VocabularyBase {
public:
	char *Name() {return "Ode";}
	char *Abbrev() {return "Ode";}
	char *Directory() {return "";}
	OdeVocabulary();
};

class Ode : public ProblemBase {
public:
	void InitialConditions();
	void Ode::Initialize();
	void Output(int it);
	void NonLinearSrc(Var *, Var *, double);
	void LinearSrc(Var *, Var *, double);
	static void ExponentialLinearity(Var *, Var *, double);
};

class E_Euler : public Euler {
protected:
	Nu *expinv,*onemexpinv;
public:
	void Allocate(int);
	char *Name() {return "Exponential Euler";}
	void TimestepDependence(double);
	Solve_RC Solve(double, double);
	void Source(Var *src, Var *Y, double t) {
		Problem->NonLinearSrc(src,Y,t);
		Ode::ExponentialLinearity(src,Y,t);
	}
};

class I_Euler : public Euler {
protected:
	Nu *expinv,*onemexpinv;
public:
	void Allocate(int);
	char *Name() {return "Euler w/Integrating Factor";}
	void TimestepDependence(double);
	Solve_RC Solve(double, double);
	void Source(Var *src, Var *Y, double t) {
		Problem->NonLinearSrc(src,Y,t);
	}
};

class RB1 : public Euler {
public:
	char *Name() {return "First-Order Rosenbrock";}
	Solve_RC Solve(double, double);
};

class I_PC : public PC {
protected:
	Nu *expinv;
public:
	void Allocate(int);
	char *Name() {return "Predictor-Corrector w/Integrating Factor";}
	void TimestepDependence(double);
	void Predictor(double, double, int, int);
	int Corrector(double, int, int, int);
	void Source(Var *src, Var *Y, double t) {
		Problem->NonLinearSrc(src,Y,t);
	}
};

class E_PC : public PC {
protected:
	Nu *expinv,*onemexpinv;
	double dtinv;
public:
	void Allocate(int);
	char *Name() {return "Exponential Predictor-Corrector";}
	void TimestepDependence(double);
	void Predictor(double, double, int, int);
	int Corrector(double, int, int, int);
	void Source(Var *src, Var *Y, double t) {
		Problem->NonLinearSrc(src,Y,t);
		Ode::ExponentialLinearity(src,Y,t);
	}
};

OdeVocabulary::OdeVocabulary()
{
	Vocabulary=this;
	
	VOCAB(nu0,0.0,0.0);
	
	METHOD(Ode);
	
	INTEGRATOR(I_Euler);
	INTEGRATOR(E_Euler);
	INTEGRATOR(RB1);
	INTEGRATOR(I_PC);
	INTEGRATOR(E_PC);
}

OdeVocabulary Ode_Vocabulary;

ofstream fout;

void Ode::InitialConditions()
{
	ny=1;
	Nmoment=0;
	y=new Var[ny];
	nu=new Var[ny];
	
	y[0]=1.0;
	nu[0]=nu0;
	
	dynamic=0;
	open_output(fout,dirsep,downcase(undashify(Integrator->Abbrev())));
}

void Ode::Initialize()
{
	fout << "# " << Integrator->Name() << endl;
}

void Ode::Output(int)
{
	fout << t << "\t" << y[0] << "\t" << endl;
}

void Ode::NonLinearSrc(Var *source, Var *, double)
{
	source[0]=cos(y[0]);
}

void Ode::LinearSrc(Var *source, Var *y, double)
{
	source[0] -= nu[0]*y[0];
}

void Ode::ExponentialLinearity(Var *source, Var *, double)
{
	source[0] *= nu_inv[0];
}

void E_Euler::Allocate(int n)
{
	Euler::Allocate(n);
	expinv=new Nu[nyprimary];
	onemexpinv=new Nu[nyprimary];
	
	nu_inv=new Nu[nyprimary];

	for(int j=0; j < nyprimary; j++) {
		if(nu[j] != 0.0) nu_inv[j]=1.0/nu[j];
		else nu_inv[j]=1.0;
	}
}

void E_Euler::TimestepDependence(double dt)
{
	for(int j=0; j < nyprimary; j++) {
		onemexpinv[j]=-expm1(-nu[j]*dt);
		expinv[j]=1.0-onemexpinv[j];
		if(nu[j] == 0.0) onemexpinv[j]=dt;
	}
}

Solve_RC E_Euler::Solve(double t, double dt)
{
	Source(source,y0,t);
	Problem->Transform(y0,t,dt,yi);
	for(int j=0; j < nyprimary; j++)
		y0[j]=expinv[j]*y0[j]+onemexpinv[j]*source[j];
	Problem->BackTransform(y0,t+dt,dt,yi);
	return SUCCESSFUL;
}

void I_Euler::Allocate(int n)
{
	Euler::Allocate(n);
	expinv=new Nu[nyprimary];
}

void I_Euler::TimestepDependence(double dt)
{
	for(int j=0; j < nyprimary; j++) expinv[j]=exp(-nu[j]*dt);
}

Solve_RC I_Euler::Solve(double t, double dt)
{
	Source(source,y0,t);
	Problem->Transform(y0,t,dt,yi);
	for(int j=0; j < nyprimary; j++)
		y0[j]=(y0[j]+dt*source[j])*expinv[j];
	Problem->BackTransform(y0,t+dt,dt,yi);
	return SUCCESSFUL;
}

Solve_RC RB1::Solve(double t, double dt)
{
	Source(source,y0,t);
	Problem->Transform(y0,t,dt,yi);
	for(int j=0; j < nyprimary; j++)
		y0[j]=y0[j]+dt*source[0]/(1.0-dt*(-sin(y0[j])-nu[0]));
	Problem->BackTransform(y0,t+dt,dt,yi);
	return SUCCESSFUL;
}

void I_PC::Allocate(int n)
{
	PC::Allocate(n);
	expinv=new Nu[nyprimary];
}

void I_PC::TimestepDependence(double dt)
{
	for(int j=0; j < nyprimary; j++) expinv[j]=exp(-nu[j]*dt);
}

void I_PC::Predictor(double t, double dt, int start, int stop)
{
	for(int j=start; j < stop; j++) y1[j]=(y0[j]+dt*source0[j])*expinv[j];
	Source(source,y1,t+dt);
}

int I_PC::Corrector(double dt, int dynamic, int start, int stop)
{
	const double halfdt=0.5*dt;
	if(dynamic) for(int j=start; j < stop; j++) {
		Var pred=y[j];
		y[j]=y0[j]*expinv[j]+halfdt*(source0[j]*expinv[j]+source[j]);
		if(!errmask || errmask[j]) 
			CalcError(y0[j]*expinv[j],y[j],pred,y[j]);
	} else for(int j=start; j < stop; j++) {
		y[j]=y0[j]*expinv[j]+halfdt*(source0[j]*expinv[j]+source[j]);
	}
	return 1;
}

void E_PC::Allocate(int n)
{
	PC::Allocate(n);
	expinv=new Nu[nyprimary];
	onemexpinv=new Nu[nyprimary];
	
	nu_inv=new Nu[nyprimary];

	for(int j=0; j < nyprimary; j++) {
		if(nu[j] != 0.0) nu_inv[j]=1.0/nu[j];
		else nu_inv[j]=1.0;
	}
}

void E_PC::TimestepDependence(double dt)
{
	for(int j=0; j < nyprimary; j++) {
		onemexpinv[j]=-expm1(-nu[j]*dt);
		expinv[j]=1.0-onemexpinv[j];
		if(nu[j] == 0.0) onemexpinv[j]=dt;
	}
	dtinv=1.0/dt;
}

void E_PC::Predictor(double t, double, int start, int stop)
{
	for(int j=start; j < stop; j++)	
		y1[j]=expinv[j]*y0[j]+onemexpinv[j]*source0[j];
	Source(source,y1,t+dt);
}

int E_PC::Corrector(double, int dynamic, int start, int stop)
{
	int j;
	for(j=start; j < stop; j++) {
		source[j]=0.5*(source0[j]+source[j]);
		y[j]=expinv[j]*y0[j]+onemexpinv[j]*source[j];
	}
	if(dynamic)
		for(j=start; j < stop; j++)
			if(!errmask || errmask[j]) 
				CalcError(y0[j]*dtinv,source[j],source0[j],source[j]);
	return 1;
}
