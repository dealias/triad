#ifndef __kernel_h__
#define __kernel_h__ 1

#include <iostream>
#include <fstream>
#include <cstdio>
#include <climits>
#include <cerrno>
#include <ctime>
#include <string>
#include <sstream>

#include "utils.h"
#include "Table.h"

#ifdef NDEBUG
const bool DEBUG=false;
#else
const bool DEBUG=true;
#endif
void check_compatibility(const bool debug);

// Global variables
extern const double ProblemVersion;
extern double t;
extern double last_dump;
extern long long iteration;
extern int invert_cnt;

extern const char *run;
extern const char *method;
extern const char *integrator;

// Global vocabulary
extern long long itmax; 
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
extern int threads;

using Array::array1;
using Array::Allocate;
using Array::Dimension;
using Array::Set;

typedef array1<Var>::opt vector;
typedef array1<Real>::opt rvector;
typedef array1<int>::opt ivector;
typedef array1<unsigned int>::opt uvector;

typedef array1<vector> vector2;
typedef array1<vector2> vector3;
typedef array1<rvector> rvector2;

enum Solve_RC {NONINVERTIBLE=-1,UNSUCCESSFUL,SUCCESSFUL,ADJUST};

class ProblemBase {
 protected:
  vector y; // Array of all field data
  DynVector<unsigned int> NY; // number of variables in each field
  array1<vector> Y; // array of dependent fields
  uvector index; // array of offsets to start of each field
  unsigned int ny;
  const char *abbrev;
  ivector errmask;
  unsigned int nfields;
  bool stochastic;
  bool prepareOutput;
  
 public:	
  size_t align;
  
  ProblemBase() : stochastic(false), prepareOutput(false), align(0) {
    Array::Null(errmask);
  }
  virtual ~ProblemBase() {}
  void SetAbbrev(const char *abbrev0) {abbrev=abbrev0;}
  const char *Abbrev() {return abbrev;}
  vector yVector() {return y;}
  vector2 YVector() {return Y;}
  unsigned int Size() {return ny;}
  unsigned int Size(int field) {return NY[field];}
  unsigned int Start(unsigned int field) {return field < nfields ? index[field] : 0;}
  unsigned int Stop(unsigned int field) {return field < nfields ? index[field]+NY[field] : 0;}
  DynVector<unsigned int>* Sizes() {return &NY;}
  
  unsigned int NFields() {return nfields;}

  const ivector& ErrorMask() {return errmask;}
  
  void Allocator(size_t Align=0) {
    align=Align;
    nfields=NY.Size();
    Y.Allocate(nfields);
    Allocate(index,nfields);
    ny=0;
    for(unsigned int i=0; i < nfields; i++) ny += NY[i];
    Allocate(y,ny,align);
    Var *p=y;
    unsigned int count=0;
    for(unsigned int i=0; i < nfields; i++) {
      unsigned int n=NY[i];
      Dimension(Y[i],n,p);
      index[i]=count;
      count += n;
      p += n;
    }
  }
	
  virtual const char *Name() {return "";}
	
  virtual void Source(const vector2& Src, const vector2& Y, double t)=0;
  virtual void Transform(const vector2&, double, double, vector2&) {}
  virtual void BackTransform(const vector2&, double, double, const vector2&) {}
  virtual void Stochastic(const vector2&, double, double) {}
  virtual void Stochastic(bool b) {stochastic=b;}
  virtual bool Stochastic() {return stochastic;}
  virtual void Initialize() {}
  virtual void Setup() {}
  virtual void FinalOutput() {}
  virtual int Microprocess() {return 0;}
	
  virtual void InitialConditions() {};
  virtual void Output(int) {};
  virtual void PrepareOutput(bool b) {prepareOutput=b;}
  virtual bool PrepareOutput() {return prepareOutput;}
};
	
Compare_t ProblemCompare;
KeyCompare_t ProblemKeyCompare;
extern ProblemBase *Problem;

void poll();
void read_init();
void set_timer();
void statistics(double t, double dt, int it);
void lock();
void unlock();
void testlock();
void dump(double t, int it, int final, double tmax);
void SaveParameters();
  
#include "Param.h"

#endif
