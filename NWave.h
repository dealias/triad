#ifndef __NWave_h__
#define __NWave_h__ 1

#include "kernel.h"
#include "Geometry.h"
#include "Partition.h"
#include "Basis.h"
#include "Linearity.h"

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
	Table<LinearityBase> *LinearityTable;
	GeometryBase *NewGeometry(char *key) {
		return GeometryTable->Locate(key);
	}
	LinearityBase *NewLinearity(char *key) {
		return LinearityTable->Locate(key);
	}
};

class NWave : public ProblemBase {
public:
	NWave() {}
	void InitialConditions();
	void Initialize();
	void OpenOutput();
	void Output(int it);
	void FinalOutput();
	void LinearSrc(Var *, Var *, double);
	
	static void ExponentialLinearity(Var *, Var *, double);
	static void ConservativeExponentialLinearity(Real *, Real *, double );
	static void ConservativeExponentialLinearity(Complex *, Complex *, double);
};

class SR : public NWave {
public:	
	SR() {}
	char *Name() {return "Spectral Reduction";}
	void NonLinearSrc(Var *source, Var *psi, double);
};

class Convolution : public NWave {
public:
	Convolution() {}
	char *Name() {return "Convolution";}
	void NonLinearSrc(Var *source, Var *psi, double);
};

class PS : public NWave {
public:
	PS() {
		if(!reality) msg(ERROR,"Pseudospectral approximation needs reality=1");
	}
	char *Name() {return "Pseudospectral";}
	void NonLinearSrc(Var *source, Var *psi, double);
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
	int Corrector(double, int, int, int);
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
	void Source(Var *src, Var *Y, double t) {Problem->NonLinearSrc(src,Y,t);}
};

#if _CRAY
typedef NWave _dummy_;
#endif

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
		NWave::ExponentialLinearity(src,Y,t);
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
	int Corrector(double, int, int, int);
	void Source(Var *src, Var *Y, double t) {
		Problem->NonLinearSrc(src,Y,t);
		NWave::ConservativeExponentialLinearity(src,Y,t);
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
	int Corrector(double, int, int, int);
	void Source(Var *src, Var *Y, double t) {Problem->NonLinearSrc(src,Y,t);}
};

class C_RK2 : public RK2, public CorrectC_PC {
public:
	void Allocate(int n) {RK2::Allocate(n); if(hybrid) y1=new Var[n];}
	char *Name() {return "Conservative Second-Order Runge-Kutta";}
	int Corrector(double, int, int, int);
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
	int Corrector(double, int, int, int);
	void Source(Var *src, Var *Y, double t) {Problem->NonLinearSrc(src,Y,t);}
};

class C_RK4 : public RK4 {
protected:
	double sixthdt2;
public:
	char *Name() {return "Conservative Fourth-Order Runge-Kutta";}
	void TimestepDependence(double);
	int Corrector(double, int, int, int);
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
	int Corrector(double, int, int, int);
	void Source(Var *src, Var *Y, double t) {Problem->NonLinearSrc(src,Y,t);}
};

class C_RK5 : public RK5 {
public:
	void Allocate(int n) {RK5::Allocate(n); y2=new Var[n]; y4=new Var[n];}
	char *Name() {return "Conservative Fifth-Order Runge-Kutta";}
	int Corrector(double, int, int, int);
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
