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
extern const double problem_version;
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

using Array::array1;
using Array::Allocate1;
using Array::Dimension1;
using Array::Set1;

typedef array1<Var>::opt vector;
typedef array1<Real>::opt rvector;
typedef array1<int>::opt ivector;
typedef array1<unsigned int>::opt uvector;

typedef array1<vector> vector2;

class ProblemBase {
 protected:
  vector y; // Array of all field data
  DynVector<unsigned int> NY; // number of variables in each field
  array1<vector > Y; // array of dependent fields
  unsigned int ny;
  const char *abbrev;
  ivector errmask;
 public:	
  ProblemBase() {errmask=__NULLARRAY}
  virtual ~ProblemBase() {}
  void SetAbbrev(const char *abbrev0) {abbrev=abbrev0;}
  const char *Abbrev() {return abbrev;}
  vector Vector() {return y;}
  vector2 Vector2() {return Y;}
  unsigned int Size() {return ny;}
  DynVector<unsigned int> *Index() {return &NY;}

  const ivector& ErrorMask() {return errmask;}
  
  void Allocate() {
    unsigned int nfields=NY.Size();
    Y.Allocate(nfields);
    ny=0;
    for(unsigned int i=0; i < nfields; i++) ny += NY[i];
    Allocate1(y,ny);
    Var *p=y;
    for(unsigned int i=0; i < nfields; i++) {
      Dimension1(Y[i],NY[i],p);
      p += NY[i];
    }
  }
	
  virtual const char *Name() {return "";}
	
  virtual void Source(const vector2& Src, const vector2& Y, double t)=0;
  virtual void Transform(const vector2&, double, double, vector2&) {}
  virtual void BackTransform(const vector2&, double, double, const vector2&) {}
  virtual void Stochastic(const vector2&, double, double) {}
  virtual void Initialize() {}
  virtual void Setup() {}
  virtual void FinalOutput() {}
  virtual int Microprocess() {return 0;}
	
  virtual void InitialConditions() {};
  virtual void Output(int) {};
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
void save_parameters();

#include "Param.h"

#endif
