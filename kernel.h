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

#define __ArrayExtensions
#define __ExternalArrayExit
void __ArrayExit(char *x);
	
// Global variables
extern double t;
extern double last_dump;
extern int iteration;
extern int invert_cnt;

extern const char *run;
extern const char *method;
extern const char *integrator;

// Global vocabulary
extern int itmax; 
extern double tmax;
extern double dt;
extern double tprecision;
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
  virtual void SetStr(const char *)=0;		// Set from string	
  virtual const char *Name()=0;
};

enum Solve_RC {NONINVERTIBLE=-1,UNSUCCESSFUL,SUCCESSFUL,ADJUST};

class ProblemBase {
 protected:
  Var *y;
  unsigned int ny;
  const char *abbrev;
  int *errmask;
 public:	
  ProblemBase() {errmask=NULL;}
  void SetAbbrev(const char *abbrev0) {abbrev=abbrev0;}
  const char *Abbrev() {return abbrev;}
  Var *Vector() {return y;}
  unsigned int Size() {return ny;}
  int *ErrorMask() {return errmask;}
	
  virtual const char *Name() {return "";}
	
  virtual void Source(Var *src, Var *y, double t) {};
  virtual void Transform(Var *, double, double, Var *&) {}
  virtual void BackTransform(Var *, double, double, Var *) {}
  virtual void Stochastic(Var *, double, double) {}
  virtual void Initialize() {}
  virtual void FinalOutput() {}
  virtual int Microprocess() {return 0;}
	
  virtual void InitialConditions() {};
  virtual void Output(int it) {};
};
	
Compare_t ProblemCompare;
KeyCompare_t ProblemKeyCompare;
extern ProblemBase *Problem;

class IntegratorBase {
 protected:
  const char *abbrev;
  ProblemBase *Problem;
  unsigned int ny;
  Var *y0, *yi, *source;
  double errmax;
  int *errmask;
  double tolmax2,tolmin2;
  double stepfactor,stepinverse,stepnoninverse;
  double growfactor, shrinkfactor;
  double dtmin,dtmax;
  int itmax,microsteps;
  int microprocess;
  int verbose;
  int dynamic;
 public:	
  void SetAbbrev(const char *abbrev0) {abbrev=abbrev0;}
  const char *Abbrev() {return abbrev;}
  void SetParam(double tolmax, double tolmin, double stepfactor0,
		double stepnoninvert, double dtmin0, double dtmax0,
		int itmax0, int microsteps0, int verbose0, int dynamic0) {
    if(tolmax < tolmin) msg(ERROR_GLOBAL,"tolmax < tolmin"); 
    tolmax2=tolmax*tolmax;
    tolmin2=tolmin*tolmin;
    growfactor=stepfactor=stepfactor0;
    shrinkfactor=stepinverse=1.0/stepfactor;
    stepnoninverse=1.0/stepnoninvert;
    dtmin=dtmin0;
    dtmax=dtmax0;
    itmax=itmax0;
    microsteps=microsteps0*Microfactor();
    verbose=verbose0;
    dynamic=dynamic0;
  }
  void Integrate(ProblemBase& problem, Var *const y, double& t, double tmax,
		 double& dt, const double sample, int& iteration);
  void ChangeTimestep(double& dt, double dtnew, const double t,
		      const double sample);
	
  virtual void Source(Var *src, Var *y, double t) {Problem->Source(src,y,t);}
	
  void CalcError(const Var& initial, const Var& norm, const Var& pred,
		 const Var& corr);
  Solve_RC CheckError();
  virtual void Allocate(int n);
  virtual const char *Name()=0;
  virtual Solve_RC Solve(double, double&)=0;
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
  virtual ~VocabularyBase() {}
  Table<ProblemBase> *ProblemTable;
  Table<IntegratorBase> *IntegratorTable;
	
  ParamBase *Locate(const char *key, int *match_type);
  void ParamAdd(ParamBase *p);
  void Parse(char *s);
  void Assign(const char *s, int warn=1);
  void Sort();
  void List(ostream& os);
  void Dump(ostream& os);
  void GraphicsDump(ostream& os);
  virtual const char *Name()=0;
  virtual const char *Abbrev()=0;
  virtual const char *Directory() {return "";}

  ProblemBase *NewProblem(const char *& key) {
    ProblemBase *p=ProblemTable->Locate(key);
    p->SetAbbrev(key);
    return p;
  }
	
  IntegratorBase *NewIntegrator(const char *& key) {
    char *key2=strdup(key);
    undashify(key,key2);
    const char *key0=key2;
    IntegratorBase *p=IntegratorTable->Locate(key0);
    p->SetAbbrev(key0);
    return p;
  }
	
  virtual const char *FileName(const char* delimiter="", 
			       const char *suffix="");
};

extern VocabularyBase *Vocabulary;

#include "Param.h"
#include "Integrator.h"

#define METHOD(key) (void) new Entry<key,ProblemBase> (#key,ProblemTable);

#define PLURAL(x) ((x)==1 ? "" : "s")

int poll();
void read_init();
void set_timer();
void statistics(int it);
void lock();
void unlock();
void testlock();
void dump(int it, int final, double tmax);

template<class T>
inline void open_output(T& fout, const char *delimiter, const char *suffix,
			int append)
{
  const char *filename=Vocabulary->FileName(delimiter,suffix);
  if(append) fout.open(filename,fout.app); // Append to end of output file.
  else fout.open(filename);
  if(!fout) msg(ERROR,"Output file %s could not be opened",filename);
  fout.precision(digits);
  errno=0;
}

template<class T>
inline void open_output(T& fout, const char *delimiter, const char *suffix)
{
  open_output(fout,delimiter,suffix,restart);
}	

#endif
