#ifndef __kernel_h__
#define __kernel_h__ 1

#include <iostream.h>
#include <fstream.h>
#include <stdio.h>
#include <limits.h>
#include <errno.h>
#include <time.h>
#include <string.h>
#include <strstream.h>

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
	int *errmask;
public:	
	ProblemBase() {errmask=NULL;}
	void SetAbbrev(char *abbrev0) {abbrev=abbrev0;}
	char *Abbrev() {return abbrev;}
	Var *Vector() {return y;}
	int Size() {return ny;}
	int *ErrorMask() {return errmask;}
	
	virtual char *Name() {return "";}
	virtual void LinearSrc(Var *, Var *, double) {}
	virtual void NonLinearSrc(Var *, Var *, double) {}
	virtual void Transform(Var *, double, double, Var *&) {}
	virtual void BackTransform(Var *, double, double, Var *) {}
	virtual void Initialize() {}
	virtual void FinalOutput() {}
	virtual int Microprocess() {return 0;}
	
	virtual void InitialConditions()=0;
	virtual void Output(int it)=0;
};
	
Compare_t ProblemCompare;
KeyCompare_t ProblemKeyCompare;
extern ProblemBase *Problem;

class IntegratorBase {
protected:
	char *abbrev;
	int ny,nyprimary;
	Var *y0, *yi, *source;
	double errmax;
	int *errmask;
	double tolmax2,tolmin2;
	double stepfactor,stepinverse,stepnoninverse;
	double dtmax;
	int itmax,microsteps;
	int microprocess;
	int verbose;
public:	
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
	void ChangeTimestep(double& dt, double dtnew, const double t,
						const double sample);
	
	void CalcError(const Var& initial, const Var& norm, const Var& pred,
				   const Var& corr);
	Solve_RC CheckError();
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
	virtual char *Directory() {
		strstream buf;
		buf << Abbrev() << dirsep << ends;
		buf.freeze();
		return buf.str();
	}

	ProblemBase *NewProblem(char *& key) {
		ProblemBase *p=ProblemTable->Locate(key);
		p->SetAbbrev(key);
		return p;
	}
	
	IntegratorBase *NewIntegrator(char *& key) {
		undashify(key,key);
		IntegratorBase *p=IntegratorTable->Locate(key);
		p->SetAbbrev(key);
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

// Cfront can't seem to handle a template here.
inline void open_output(ofstream& fout, const char *delimiter, char *suffix,
						int append=restart)
{
	char *filename=Vocabulary->FileName(delimiter,suffix);
	if(append) fout.open(filename,fout.app); // Append to end of output file.
	else fout.open(filename);
	if(!fout) msg(ERROR,"Output file %s could not be opened",filename);
	fout.precision(digits);
}

inline void open_output(oxstream& fout, const char *delimiter, char *suffix,
						int append=restart)
{
	char *filename=Vocabulary->FileName(delimiter,suffix);
	if(append) fout.open(filename,fout.app); // Append to end of output file.
	else fout.open(filename);
	if(!fout) msg(ERROR,"Output file %s could not be opened",filename);
}

#endif
