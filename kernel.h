#ifndef __kernel_h__
#define __kernel_h__ 1

#include <iostream.h>
#include <fstream.h>
#include <stdio.h>
#include <limits.h>
#include <errno.h>
#include <time.h>
#include <string.h>

#include "utils.h"
#include "types.h"
#include "DynVector.h"
#include "Table.h"

// Global variables
extern double t;
extern double last_dump;
extern int iteration;
extern int invert_cnt;

extern char *run;
extern char *integrator;
extern char *approximation;

// Global vocabulary
extern int itmax; 
extern double tmax;
extern double dt;
extern int average;
extern int verbose;
extern int dynamic;
extern int digits;
extern int restart;
extern double polltime;
extern int output;
extern int hybrid;
extern int override;
extern int discrete;

class ParamBase {
	int nvar;
public:
	ParamBase() {}
	virtual void Display(ostream& os)=0;
	virtual void GraphicsOutput(ostream& os)=0;
	virtual void Output(ostream& os)=0;
	virtual void SetStr(char *)=0;		// Set from string	
	virtual char *Name()=0;
};

enum Solve_RC {NONINVERTIBLE=-1,UNSUCCESSFUL,SUCCESSFUL,ADJUST};

typedef void Source_t(Var *, Var *, double);

class IntegratorBase {
protected:
	char *abbrev;
	int ny,nyconserve;
	Var *source;
	double tolmax2,tolmin2;
	double stepfactor,stepinverse,stepnoninverse;
	double dtmax;
	int itmax,microsteps;
	int microprocess;
	int verbose;
	typedef void Source_t(Var *, Var *, double);
	Source_t *LinearSrc,*NonlinearSrc,*ConstantSrc;

public:	
	void SetAbbrev(char *abbrev0) {abbrev=abbrev0;}
	char *Abbrev() {return abbrev;}
	void SetParam(double tolmax, double tolmin, double stepfactor0,
				  double stepnoninvert, double dtmax0, int itmax0,
				  int microsteps0, int nyconserve0, int verbose0) {
		if(tolmax < tolmin) msg(ABORT,"tolmax < tolmin"); 
		tolmax2=tolmax*tolmax;
		tolmin2=tolmin*tolmin;
		stepfactor=stepfactor0;
		stepinverse=1.0/stepfactor;
		stepnoninverse=1.0/stepnoninvert;
		dtmax=dtmax0;
		itmax=itmax0;
		microsteps=microsteps0*Microfactor();
		nyconserve=nyconserve0;
		verbose=verbose0;
	}
	void Integrate(Var *const y, double& t, double tmax,
				   Source_t *const LinearSrc0, Source_t *const NonlinearSrc0,
				   Source_t *const ConstantSrc, double& dt,
				   const double sample);
	void ChangeTimestep(double& dt, const double dtnew, const double t);
	
	Solve_RC CheckError(double errmax);
	virtual void Allocate(int)=0;
	virtual char *Name()=0;
	virtual Solve_RC Solve(Var *, double, double)=0;
	virtual void Source(Var *, Var *, double);
	virtual int Microfactor() {return 1;}
	virtual void TimestepDependence(double) {}
};

Compare_t IntegratorCompare;
KeyCompare_t IntegratorKeyCompare;
extern IntegratorBase *Integrator;

class ApproximationBase {
protected:	
	char *abbrev;
public:	
	void SetAbbrev(char *abbrev0) {abbrev=abbrev0;}
	char *Abbrev() {return abbrev;}
	virtual char *Name()=0;
	virtual void SetSrcRoutines(Source_t **LinearSrc, Source_t **NonlinearSrc,
								Source_t **ConstantSrc)=0;
};

Compare_t ApproximationCompare;
KeyCompare_t ApproximationKeyCompare;
extern ApproximationBase *Approximation;

class ProblemBase {
protected:
	Var *y;
	DynVector<ParamBase *> ParamList;
	int ny,nyconserve;
public:	
	Table<IntegratorBase> *IntegratorTable;
	Table<ApproximationBase> *ApproximationTable;
	
	ProblemBase();
	
	ParamBase *Locate(char *key, int *match_type);
	void ParamAdd(ParamBase *p);
	void Parse(char *s);
	void Assign(const char *s);
	void Sort();
	void List(ostream& os);
	void Dump(ostream& os);
	void GraphicsDump(ostream& os);

	IntegratorBase *NewIntegrator(char *key) {
		char *abbrev=undashify(key);
		IntegratorBase *p=IntegratorTable->Locate(abbrev);
		p->SetAbbrev(abbrev);
		return p;
	}
	
	ApproximationBase *NewApproximation(char *key) {
		ApproximationBase *p=ApproximationTable->Locate(key);
		p->SetAbbrev(upcase(key));
		return p;
	}
	
	Var *Vector() {return y;}
	int Size() {return ny;}
	int Nconserve() {return nyconserve;}
	virtual char *Name()=0;
	virtual char *Abbrev()=0;
	virtual void InitialConditions()=0;
	virtual void Initialize() {};
	virtual void Output(int it)=0;
	virtual void FinalOutput() {}
	virtual char *FileName(const char* delimiter="", const char *suffix="");
	virtual int Microprocess() {return 0;}
};

extern ProblemBase *Problem;

#include "Param.h"
#include "Integrator.h"
#include "Approx.h"

#define PLURAL(x) ((x)==1 ? "" : "s")

int poll();
void read_init();
void set_timer();
void statistics();
void lock();
void unlock();
void testlock();
void dump(int it, int final, double tmax);
void open_output(ofstream& fout, const char *delimiter, char *suffix,
				 int append=restart);

#endif
