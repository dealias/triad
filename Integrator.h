#ifndef __Integrator_h__
#define __Integrator_h__ 1

enum Solve_RC {NONINVERTIBLE=-1,UNSUCCESSFUL,SUCCESSFUL,ADJUST};

class IntegratorBase {
 protected:
  const char *abbrev;
  ProblemBase *Problem;
  unsigned int ny;
  vector y,source;
  vector2 Y,Src,YI,Yout; // arrays of dependent fields
  double t;
  double dt;
  double sample;
  double errmax;
  ivector errmask;
  double tolmax2,tolmin2;
  double stepfactor,stepinverse,stepnoninverse;
  double growfactor,shrinkfactor;
  double dtmin,dtmax;
  int itmax,microsteps;
  int microprocess;
  int verbose;
  int dynamic;
 public:	
  virtual ~IntegratorBase() {}
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
  void Integrate(double& t0, double tmax, double& dt0, double sample0,
		 int& iteration, unsigned long& nout);
  void ChangeTimestep(double dtnew);
	
  virtual void Source(const vector2& Src, 
		      const vector2& Y, double t) {
    Problem->Source(Src,Y,t);
  }
	
  void CalcError(const Var& initial, const Var& norm, const Var& pred,
		 const Var& corr);
  Solve_RC CheckError();
  
  void Alloc0(vector2& Y, const vector& y);
  void Alloc(vector2& Y, vector& y);
  
  void Alloc(vector2& A, vector& a, const vector2& B, const vector& b) {
    A.Dimension(B);
    Dimension1(a,b);
  }
  
  virtual void SetProblem(ProblemBase& problem);
  virtual void Allocate();
  virtual const char *Name()=0;
  virtual Solve_RC Solve()=0;
  virtual int Microfactor() {return 1;}
  virtual void TimestepDependence() {}
};

Compare_t IntegratorCompare;
KeyCompare_t IntegratorKeyCompare;
extern IntegratorBase *Integrator;

#define INTEGRATOR(key) \
{(void) new Entry<key,IntegratorBase>(#key,IntegratorTable);}

inline void IntegratorBase::CalcError(const Var& initial, const Var& norm0, 
				      const Var& pred, const Var& corr)
{
  if(initial != 0.0 && pred != initial) {
    Real error=max(divide0(abs2(corr-pred),
			   max(abs2(norm0),abs2(initial))));
    if(error > errmax) errmax=error;
  }
}

inline Solve_RC IntegratorBase::CheckError()
{
  if(dynamic >= 0) {
    if(errmax > tolmax2) return UNSUCCESSFUL;
    if(errmax < tolmin2) return ADJUST;
    return SUCCESSFUL;
  }
  if(++dynamic == 0) dynamic=1;
  return SUCCESSFUL;
}

class Euler : public IntegratorBase {
 public:
  const char *Name() {return "Euler";}
  Solve_RC Solve();
};

class SYM1 : public IntegratorBase {
 public:
  const char *Name() {return "First-Order Symplectic";}
  Solve_RC Solve();
};

class Midpoint : public IntegratorBase {
 protected:	
  vector y0;
  vector2 Y0;
  double halfdt;
 public:
  void Allocate() {
    IntegratorBase::Allocate();
    Alloc(Y0,y0);
    halfdt=0.5*dt;
  }
  const char *Name() {return "Midpoint Rule";}
  void TimestepDependence() {
    halfdt=0.5*dt;
  }
  Solve_RC Solve();
};

class AdamsBashforth : public IntegratorBase {
  vector y0,source0,source1;
  vector2 Y0,Src0,Src1;
  double a0,a1,a2;
  double b0,b1,b2;
  int init;
 public:
  void Allocate() {
    IntegratorBase::Allocate();
    Alloc(Y0,y0);
    Alloc(Src0,source0);
    Alloc(Src1,source1);
    init=2;
  }
  const char *Name() {return "Third-Order Adams-Bashforth";}
  Solve_RC Solve();
  void TimestepDependence() {
    a0=23.0/12.0*dt;
    a1=-4.0/3.0*dt;
    b0=a2=5.0/12.0*dt;
    b1=2.0/3.0*dt;
    b2=-1.0/12.0*dt;
  }
};

class PC : public IntegratorBase {
 protected:
  int new_y0;
  vector y0,source0;
  vector2 Y0,Src0;
  double halfdt;
  int order;
  double pgrow, pshrink;
 public:
  PC() {order=2;}
  void Allocate() {
    IntegratorBase::Allocate();
    Alloc(Y0,y0);
    Alloc(Src0,source0);
    new_y0=1;
    pgrow=0.5/order; pshrink=(order > 1) ? 0.5/(order-1) : 0;
  }
  const char *Name() {return "Predictor-Corrector";}
  Solve_RC Solve();
	
  void TimestepDependence() {
    halfdt=0.5*dt;
  }
  virtual void ExtrapolateTimestep () {
    if(errmax < tolmin2) {
      if(errmax) growfactor=pow(tolmin2/errmax,pgrow);
    } else if(tolmin2) shrinkfactor=growfactor=pow(tolmin2/errmax,pshrink);
    growfactor=min(growfactor,stepfactor);
    shrinkfactor=max(shrinkfactor,stepinverse);
    if(errmax <= tolmax2) errmax=0.0; // Force a time step adjustment.
  }
  virtual void Predictor();
  virtual int Corrector();
};

class LeapFrog : public PC {
  vector yp,yp0;
  vector2 YP,YP0; // array of dependent fields
  double oldhalfdt,lasthalfdt;
 public:
  void Allocate() {
    PC::Allocate(); 
    Alloc(YP,yp);
    Alloc(YP0,yp0);
    lasthalfdt=0.0;
  }
  void TimestepDependence() {
    if(lasthalfdt == 0.0) {
      Set1(yp,YP[0]);
      Set1(y,Y[0]);
      for(unsigned int j=0; j < ny; j++) yp[j] = y[j];
    }
    halfdt=0.5*dt;
  }
  const char *Name() {return "Leapfrog";}
  void Predictor();
  int Corrector();
};

class SYM2 : public PC {
 public:
  const char *Name() {return "Second-Order Symplectic";}
  void Predictor();
  int Corrector();
};

class RK2 : public PC {
 public:
  const char *Name() {return "Second-Order Runge-Kutta";}
  void Predictor();
  int Corrector();
};

class RK4 : public PC {
 protected:
  vector source1, source2;
  vector2 Src1,Src2;
  double sixthdt;
 public:
  RK4() {order=4;}
  void Allocate() {
    PC::Allocate();
    Alloc(Src1,source1);
    Alloc(Src2,source2);
  }
  const char *Name() {return "Fourth-Order Runge-Kutta";}
  void TimestepDependence();
  void Predictor();
  int Corrector();
};

class RK5 : public RK4 {
 protected:
  vector source3,source4;
  vector2 Src3,Src4;
  double a1,a2,a3,a4,a5;
  double b10;
  double b20,b21;
  double b30,b31,b32;
  double b40,b41,b42,b43;
  double b50,b51,b52,b53,b54;
  double c0,c2,c3,c5;
  double d0,d2,d3,d4,d5;
 public:
  RK5() {order=5;}
  void Allocate() {
    RK4::Allocate();
    Alloc(Src3,source3,Src1,source1);
    Alloc(Src4,source4);
  }
  const char *Name() {return "Fifth-Order Runge-Kutta";}
  void TimestepDependence();
  void Predictor();
  int Corrector();
};

class Exact : public RK5 {
 public:
  const char *Name() {return "Exact";}
  int Microfactor(){return 100;}
};

#endif
