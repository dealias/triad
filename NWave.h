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

extern Real continuum_factor;
void compute_invariants(Var *y, int Npsi, Real& E, Real& Z, Real& P);
void display_invariants(Real E, Real Z, Real P);

class NWaveVocabulary : public VocabularyBase {
public:
	char *Name();
	char *Abbrev();
	NWaveVocabulary();
	Table<GeometryBase> *GeometryTable;
	GeometryBase *NewGeometry(char *key) {return GeometryTable->Locate(key);}
};

enum Linearity {STANDARD,EXPONENTIAL,CONSERVATIVE_EXPONENTIAL};

class NWave : public ProblemBase {
	Source_t (NWave::*Linearity);
public:
	void InitialConditions();
	void Initialize();
	void OpenOutput();
	void Output(int it);
	void FinalOutput();
	
	void LinearSrc(Var *src, Var *y, double t) {Linearity(src,y,t);}
	
	Source_t StandardLinearity;
	Source_t ExponentialLinearity;
	void ConservativeExponentialLinearity(Real *source, Real *psi, double t);
	void ConservativeExponentialLinearity(Complex *source, Complex *psi,
										  double t);
	void SetLinearity(int type) {
		switch(type) {
		case STANDARD:
			Linearity=StandardLinearity;
			break;
		case EXPONENTIAL:
			Linearity=ExponentialLinearity;
			break;
		case CONSERVATIVE_EXPONENTIAL:
			Linearity=ConservativeExponentialLinearity;
			break;
		}
	}
		
	NWave() {SetLinearity(STANDARD);}
};

class SR : public NWave {
public:	
	SR() {}
	char *Name() {return "Spectral Reduction";}
	Source_t NonLinearSrc;
};

class Convolution : public NWave {
public:
	Convolution() {}
	char *Name() {return "Convolution";}
	Source_t NonLinearSrc;
};

class PS : public NWave {
public:
	PS() {
		if(!reality) msg(ERROR,"Pseudospectral approximation needs reality=1");
	}
	char *Name() {return "Pseudospectral";}
	Source_t NonLinearSrc;
};

class C_Euler : public Euler {
	Var *y,*lastdiff;
public:
	void Allocate(int n) {
		IntegratorBase::Allocate(n);
		ny=n; source=new Var[n]; y=new Var[n]; lastdiff=new Var[n];
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
	void Source(Var *src, Var *y, double t) {Problem->NonLinearSrc(src,y,t);}
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
	void Source(Var *src, Var *y, double t) {Problem->NonLinearSrc(src,y,t);}
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
	void Source(Var *src, Var *y, double t) {Problem->NonLinearSrc(src,y,t);}
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
	Nu *expinv,*expinv5;
public:
	void Allocate(int);
	char *Name() {return "Fifth-Order Runge-Kutta w/Integrating Factor";}
	void Predictor(double, double, int, int);
	int Corrector(double, double&, int, int);
	void Source(Var *src, Var *y, double t) {Problem->NonLinearSrc(src,y,t);}
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
