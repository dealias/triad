#ifndef __NWave_h__
#define __NWave_h__ 1

#include "kernel.h"
#include "Geometry.h"
#include "Partition.h"
#include "Basis.h"

extern int Npsi;   // number of explictly evolved modes
extern int NpsiR;  // number of reflected modes
extern int Ntotal; // total number of (evolved+reflected) modes

extern Nu *nu,*nu_inv;
extern Real *nuR_inv,*nuI;
extern Real *forcing;

Source_t PrimitiveNonlinearitySR;
Source_t PrimitiveNonlinearity;
Source_t PrimitiveNonlinearityFFT;
Source_t StandardLinearity;
Source_t ExponentialLinearity;
void ConservativeExponentialLinearity(Real *source, Real *psi, double t);
void ConservativeExponentialLinearity(Complex *source, Complex *psi, double t);
Source_t ConstantForcing;

extern Real continuum_factor;
void compute_invariants(Var *y, int Npsi, Real& E, Real& Z, Real& P);
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
	void Allocate(int n) {
		ny=n; source=new Var[n]; y=new Var[n]; lastdiff=new Real[n];
	}
	char *Name() {return "Conservative Euler";}
	Solve_RC Solve(double, double);
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
	int Corrector(double, double&, int, int);
};

class I_PC : public PC {
protected:
	Nu *expinv;
public:
	void Allocate(int);
	char *Name() {return "Predictor-Corrector w/Integrating Factor";}
	void TimestepDependence(double);
	void Predictor(double, double, int, int);
	int Corrector(double, double&, int, int);
	void Source(Var *src, Var *var, double t) {
		if(NonlinearSrc) (*NonlinearSrc)(src,var,t);
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
	int Corrector(double, double&, int start, int stop);
	void Source(Var *src, Var *var, double t) {
		if(NonlinearSrc) (*NonlinearSrc)(src,var,t);
		ExponentialLinearity(src,var,t);
	}
};


class CE_PC : public E_PC, public CorrectC_PC {
protected:
	Real *onemexpinv;
public:
	void Allocate(int);
	char *Name() {return "Conservative Exponential Predictor-Corrector";}
	void TimestepDependence(double);
	void Predictor(double, double, int, int);
	int Corrector(double, double&, int, int);
	void Source(Var *src, Var *var, double t) {
		if(NonlinearSrc) (*NonlinearSrc)(src,var,t);
		ConservativeExponentialLinearity(src,var,t);
	}
};

class I_RK2 : public RK2 {
protected:
	Nu *expinv;
public:
	void Allocate(int);
	char *Name() {return "Second-Order Runge-Kutta w/Integrating Factor";}
	void TimestepDependence(double);
	void Predictor(double, double, int, int);
	int Corrector(double, double&, int, int);
	void Source(Var *src, Var *var, double t) {
		if(NonlinearSrc) (*NonlinearSrc)(src,var,t);
	}
};

class C_RK2 : public RK2, public CorrectC_PC {
public:
	void Allocate(int n) {RK2::Allocate(n); if(hybrid) y1=new Var[n];}
	char *Name() {return "Conservative Second-Order Runge-Kutta";}
	int Corrector(double, double&, int, int);
	int Correct(const Real y0, const Real y1, Real& y,
				const Real source0, const Real source, const double dt);
	int Correct(const Complex y0, const Complex y1, Complex& y,
				const Complex source0, const Complex source, const double dt);
};

class I_RK4 : public RK4 {
protected:
	Nu *expinv;
	double dtinv,thirddt;
public:
	void Allocate(int);
	char *Name() {return "Fourth-Order Runge-Kutta w/Integrating Factor";}
	void TimestepDependence(double);
	void Predictor(double, double, int, int);
	int Corrector(double, double&, int, int);
	void Source(Var *src, Var *var, double t) {
		if(NonlinearSrc) (*NonlinearSrc)(src,var,t);
	}
};

class C_RK4 : public RK4 {
protected:
	double sixthdt2;
public:
	char *Name() {return "Conservative Fourth-Order Runge-Kutta";}
	void TimestepDependence(double);
	int Corrector(double, double&, int, int);
	int Correct(const Real y0, Real& y, const Real source0,
				const Real source1, const Real source2, 
				const Real source, const double dt);
	int Correct(const Complex y0, Complex& y, const Complex source0,
				const Complex source1, const Complex source2, 
				const Complex source, const double dt);
};

class I_RK5 : public RK5 {
protected:
	Nu *expinv1,*expinv2,*expinv3,*expinv4,*expinv5;
public:
	void Allocate(int);
	char *Name() {return "Fifth-Order Runge-Kutta w/Integrating Factor";}
	void TimestepDependence(double);
	void Predictor(double, double, int, int);
	int Corrector(double, double&, int, int);
	void Source(Var *src, Var *var, double t) {
		if(NonlinearSrc) (*NonlinearSrc)(src,var,t);
	}
};

class C_RK5 : public RK5 {
public:
	void Allocate(int n) {RK5::Allocate(n); y2=new Var[n]; y4=new Var[n];}
	char *Name() {return "Conservative Fifth-Order Runge-Kutta";}
	int Corrector(double, double&, int, int);
	inline void Correct(const Real y0, Real& y2, Real &y3,
						const Real y4, Real& y,
						const Real source0, const Real source2, 
						const Real source3, const Real source4,
						const Real source, const double dt, int& invertible);
	inline void Correct(const Complex y0, Complex& y2, Complex &y3,
						const Complex y4, Complex& y,
						const Complex source0, const Complex source2, 
						const Complex source3, const Complex source4,
						const Complex source, const double dt,int& invertible);
};


#endif
