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
#include "DynVector.h"
#include "Table.h"

// Global variables
extern double t;
extern double last_dump;
extern int iteration;
extern int invert_cnt;

extern char *run;
extern char *method;
extern char *integrator;

// Global vocabulary
extern int itmax; 
extern double tmax;
extern double dt;
extern int Nmoment;
extern int dynamic;
extern int digits;
extern int restart;
extern double polltime;
extern int output;
extern int hybrid;
extern int override;
extern int verbose;
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

class ProblemBase {
protected:
	Var *y;
	int ny;
	char *abbrev;
public:	
	void SetAbbrev(char *abbrev0) {abbrev=abbrev0;}
	char *Abbrev() {return abbrev;}
	Var *Vector() {return y;}
	int Size() {return ny;}
	
	virtual char *Name() {return "";}
	virtual void SetLinearity(Source_t *);
	virtual void LinearSrc(Var *, Var *, double) {}
	virtual void NonLinearSrc(Var *, Var *, double) {}
	virtual void Transform(Var *, double, double) {}
	virtual void BackTransform(Var *, double, double) {}
	virtual void Initialize() {}
	virtual void FinalOutput() {}
	virtual int Microprocess() {return 0;}
	
	virtual void InitialConditions()=0;
	virtual void Output(int it)=0;
};
	
Compare_t ProblemCompare;
KeyCompare_t ProblemKeyCompare;
extern ProblemBase *Problem;

typedef void Source_t(Var *, Var *, double);

class IntegratorBase {
protected:
	char *abbrev;
	int ny,nyprimary;
	Var *y0, *source;
	double tolmax2,tolmin2;
	double stepfactor,stepinverse,stepnoninverse;
	double dtmax;
	int itmax,microsteps;
	int microprocess;
	int verbose;
public:	
	typedef void Source_t(Var *, Var *, double);
	void SetAbbrev(char *abbrev0) {abbrev=abbrev0;}
	char *Abbrev() {return abbrev;}
	void SetParam(double tolmax, double tolmin, double stepfactor0,
				  double stepnoninvert, double dtmax0, int itmax0,
				  int microsteps0, int verbose0) {
		if(tolmax < tolmin) msg(ABORT,"tolmax < tolmin"); 
		tolmax2=tolmax*tolmax;
		tolmin2=tolmin*tolmin;
		stepfactor=stepfactor0;
		stepinverse=1.0/stepfactor;
		stepnoninverse=1.0/stepnoninvert;
		dtmax=dtmax0;
		itmax=itmax0;
		microsteps=microsteps0*Microfactor();
		verbose=verbose0;
	}
	void Integrate(Var *const y, double& t, double tmax, double& dt, 
				   const double sample);
	void ChangeTimestep(double& dt, const double dtnew, const double t);
	
	Solve_RC CheckError(double errmax);
	virtual void Allocate(int n);
	virtual char *Name()=0;
	virtual Solve_RC Solve(double, double)=0;
	virtual void Source(Var *, Var *, double);
	virtual int Microfactor() {return 1;}
	virtual void TimestepDependence(double) {}
};

Compare_t IntegratorCompare;
KeyCompare_t IntegratorKeyCompare;
extern IntegratorBase *Integrator;

class VocabularyBase {
protected:
	DynVector<ParamBase *> ParamList;
public:	
	VocabularyBase();
	Table<ProblemBase> *ProblemTable;
	Table<IntegratorBase> *IntegratorTable;
	
	ParamBase *Locate(char *key, int *match_type);
	void ParamAdd(ParamBase *p);
	void Parse(char *s);
	void Assign(const char *s);
	void Sort();
	void List(ostream& os);
	void Dump(ostream& os);
	void GraphicsDump(ostream& os);
	virtual char *Name()=0;
	virtual char *Abbrev()=0;

	ProblemBase *NewProblem(char *key) {
		ProblemBase *p=ProblemTable->Locate(key);
		p->SetAbbrev(upcase(key));
		return p;
	}
	
	IntegratorBase *NewIntegrator(char *key) {
		char *abbrev=undashify(key);
		IntegratorBase *p=IntegratorTable->Locate(abbrev);
		p->SetAbbrev(abbrev);
		return p;
	}
	
	virtual char *FileName(const char* delimiter="", const char *suffix="");
};

extern VocabularyBase *Vocabulary;

#include "Param.h"
#include "Integrator.h"

#define METHOD(key) {new Entry<key,ProblemBase> (#key,ProblemTable);}

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
