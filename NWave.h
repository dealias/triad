#ifndef __NWave_h__
#define __NWave_h__ 1

#include "kernel.h"
#include "Geometry.h"

extern Nu *nu,*nu_inv;
extern Real *nuR_inv,*nuI;
extern Real *forcing;
extern Var *randomforce;

const int Nmoment=7;	

Source_t PrimitiveNonlinearity;
Source_t StandardLinearity;
Source_t ExponentialLinearity;
void ConservativeExponentialLinearity(Real *source, Real *psi, double t);
void ConservativeExponentialLinearity(Complex *source, Complex *psi, double t);
Source_t ConstantForcing;

extern Real *K,*Area,*y2;
extern Real continuum_factor;
void compute_invariants(Real *y2, int Npsi, Real& E, Real& Z, Real& P);
void display_invariants(Real E, Real Z, Real P);

class NWave : public ProblemBase {
protected:
public:
	Table<GeometryBase> *GeometryTable;
	
	NWave();
	char *Name();
	char *Abbrev();
	void InitialConditions();
	void Initialize();
	void OpenOutput();
	void Output(int it);
	void FinalOutput();
	
	GeometryBase *NewGeometry(char *key) {return GeometryTable->Locate(key);}
};

extern NWave *GeometryProblem;

class C_Euler : public Euler {
	Var *y;
	Real *lastdiff;
public:
	void Allocate(int n) {ny=n; source=new Var[n]; y=new Var[n];
	lastdiff=new Real[n];}
	char *Name() {return "Conservative Euler";}
	Solve_RC Solve(Real *, double, double);
	Solve_RC Solve(Complex *, double, double);
};

class CorrectC_PC {
public:
	int Correct(const Real y0, const Real y1, Real& y,
				const Real source0, const Real source,
				const double dt);
	int Correct(const Complex y0, const Complex y1, Complex& y,
				const Complex source0, const Complex source,
				const double dt);
};

class C_PC : public PC, public CorrectC_PC {
public:
	void Allocate(int n) {PC::Allocate(n); if(hybrid) y1=new Var[n];}
	char *Name() {return "Conservative Predictor-Corrector";}
	int Corrector(Var *, double, double&, int, int);
};

class E_PC : public PC {
protected:
	Nu *expinv,*onemexpinv;
public:
	void Allocate(int);
	char *Name() {return "Exponential Predictor-Corrector";}
	void TimestepDependence(double);
	void Predictor(Var *, double, double);
	int Corrector(Var *, double, double&, int start, int stop);
	void Source(Var *src, Var *var, double t) {
		if(NonlinearSrc) (*NonlinearSrc)(src,var,t);
		ExponentialLinearity(src,var,t);
	}
};

class E_RK2 : public E_PC {
protected:
	Nu *expinv1,*onemexpinv1;
	double halfdt;
public:
	void Allocate(int);
	char *Name() {return "Exponential Second-Order Runge-Kutta";}
	void TimestepDependence(double);
	void Predictor(Var *, double, double);
	int Corrector(Var *, double, double&, int start, int stop);
};

class E_RK4 : public E_RK2 {
protected:
	Var *source1,*source2;
public:
	void Allocate(int);
	char *Name() {return "Exponential Fourth-Order Runge-Kutta";}
	void Predictor(Var *, double, double);
	int Corrector(Var *, double, double&, int start, int stop);
};

class CE_PC : public E_PC, public CorrectC_PC {
protected:
	Real *onemexpinv;
public:
	void Allocate(int);
	char *Name() {return "Conservative Exponential Predictor-Corrector";}
	void TimestepDependence(double);
	void Predictor(Var *, double, double);
	int Corrector(Var *, double, double&, int, int);
	void Source(Var *src, Var *var, double t) {
		if(NonlinearSrc) (*NonlinearSrc)(src,var,t);
		ConservativeExponentialLinearity(src,var,t);
	}
};

class C_RK2 : public RK2, public CorrectC_PC {
public:
	void Allocate(int n) {RK2::Allocate(n); if(hybrid) y1=new Var[n];}
	char *Name() {return "Conservative Second-Order Runge-Kutta";}
	int Corrector(Var *, double, double&, int, int);
	int Correct(const Real y0, const Real y1, Real& y,
				const Real source0, const Real source, const double dt);
	int Correct(const Complex y0, const Complex y1, Complex& y,
				const Complex source0, const Complex source, const double dt);
};

class C_RK4 : public RK4 {
protected:
	double sixthdt2;
public:
	char *Name() {return "Conservative Fourth-Order Runge-Kutta";}
	void TimestepDependence(double);
	int Corrector(Var *, double, double&, int, int);
	int Correct(const Real y0, Real& y, const Real source0,
				const Real source1, const Real source2, 
				const Real source, const double dt);
	int Correct(const Complex y0, Complex& y, const Complex source0,
				const Complex source1, const Complex source2, 
				const Complex source, const double dt);
};

class C_RK5 : public RK5 {
public:
	void Allocate(int n) {RK5::Allocate(n); y4=new Var[n];
	if(hybrid) y5=new Var[n];}
	char *Name() {return "Conservative Fifth-Order Runge-Kutta";}
	int Corrector(Var *, double, double&, int, int);
	int Correct(const Real y0, const Real y2, const Real y3,
				const Real y4, const Real y5, Real& y,
				const Real source0, const Real source2, 
				const Real source3, const Real source4,
				const Real source, const double dt, double& errmax);
	int Correct(const Complex y0, const Complex y2, const Complex y3,
				const Complex y4, const Complex y5, Complex& y,
				const Complex source0, const Complex source2, 
				const Complex source3, const Complex source4,
				const Complex source, const double dt, double& errmax);
};

#endif
