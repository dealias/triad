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
extern size_t iteration;
extern size_t invert_cnt;

extern const char *run;
extern const char *method;
extern const char *integrator;

// Global vocabulary
extern size_t itmax;
extern double tmax;
extern double dt;
extern double tprecision;
extern double polltime;
extern size_t dynamic;
extern size_t digits;
extern size_t restart;
extern size_t output;
extern size_t hybrid;
extern size_t override;
extern size_t verbose;
extern size_t threads;

using Array::array1;
using Array::Allocate;
using Array::Dimension;
using Array::Set;

typedef array1<Var>::opt vector;
typedef array1<Real>::opt rvector;
typedef array1<size_t>::opt uvector;

typedef array1<vector> vector2;
typedef array1<vector2> vector3;
typedef array1<rvector> rvector2;

enum Solve_RC {NONINVERTIBLE=-1,UNSUCCESSFUL,SUCCESSFUL,ADJUST};

class ProblemBase {
 protected:
  vector y; // Array of all field data
  DynVector<size_t> NY; // number of variables in each field
  array1<vector> Y; // array of dependent fields
  uvector index; // array of offsets to start of each field
  size_t ny;
  const char *abbrev;
  uvector errmask;
  size_t nfields;
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
  size_t Size() {return ny;}
  size_t Size(size_t field) {return NY[field];}
  size_t Start(size_t field) {return field < nfields ? index[field] : 0;}
  size_t Stop(size_t field) {return field < nfields ? index[field]+NY[field] : 0;}
  DynVector<size_t>* Sizes() {return &NY;}

  size_t NFields() {return nfields;}

  const uvector& ErrorMask() {return errmask;}

  void Allocator(size_t Align=0) {
    align=Align;
    nfields=NY.Size();
    Y.Allocate(nfields);
    Allocate(index,nfields);
    ny=0;
    for(size_t i=0; i < nfields; i++) ny += NY[i];
    Allocate(y,ny,align);
    Var *p=y;
    size_t count=0;
    for(size_t i=0; i < nfields; i++) {
      size_t n=NY[i];
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
  virtual size_t Microprocess() {return 0;}

  virtual void InitialConditions() {};
  virtual void Output(size_t) {};
  virtual void PrepareOutput(bool b) {prepareOutput=b;}
  virtual bool PrepareOutput() {return prepareOutput;}
};

Compare_t ProblemCompare;
KeyCompare_t ProblemKeyCompare;
extern ProblemBase *Problem;

void poll();
void read_init();
void set_timer();
void statistics(double t, double dt, size_t it);
void lock();
void unlock();
void testlock();
void dump(double t, size_t it, size_t final, double tmax);
void SaveParameters();

#include "Param.h"

#endif
