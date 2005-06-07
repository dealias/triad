#ifndef __Integrator_h__
#define __Integrator_h__ 1


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
  int order;
  double pgrow, pshrink;
  
public:	
  
  IntegratorBase() {}
  
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
	
  virtual void Source(const vector2& Src, const vector2& Y, double t) {
    Problem->Source(Src,Y,t);
  }
	
  virtual void PSource(const vector2& Src, const vector2& Y, double t) {
    Source(Src,Y,t);
  }
	
  virtual void CSource(const vector2& Src, const vector2& Y, double t) {
    Source(Src,Y,t);
  }
  
  void CalcError(const Var& initial, const Var& norm, const Var& pred,
		 const Var& corr);
  Solve_RC CheckError();
  
  virtual void ExtrapolateTimestep () {
    if(errmax < tolmin2) {
      if(errmax) growfactor=pow(tolmin2/errmax,pgrow);
    } else if(tolmin2) shrinkfactor=growfactor=pow(tolmin2/errmax,pshrink);
    growfactor=min(growfactor,stepfactor);
    shrinkfactor=max(shrinkfactor,stepinverse);
    if(errmax <= tolmax2) errmax=0.0; // Force a time step adjustment.
  }
  
  void Alloc0(vector2& Y, vector& y);
  void Alloc(vector2& Y, vector& y);
  
  void Alloc(vector2& A, vector& a, const vector2& B, const vector& b) {
    A.Dimension(B);
    Dimension(a,b);
  }
  
  void SetProblem(ProblemBase& problem);
  
  unsigned int Start(int field) {return Problem->Start(field);}
  unsigned int Stop(int field) {return Problem->Stop(field);}
  
  void Allocator(const vector2& Y0, const ivector& mask);
  
  virtual void Allocator() {
    Allocator(Problem->YVector(),Problem->ErrorMask());
  }
  virtual void Allocator(ProblemBase& problem) {
    SetProblem(problem);
    Allocator();
  }
  
  virtual const char *Name()=0;
  virtual Solve_RC Solve()=0;
  virtual int Microfactor() {return 1;}
  virtual void TimestepDependence() {}
  
  virtual void Unswap() {if(Yout[0] != Y[0]) set(Yout[0],Y[0],ny);}
  
  void SetTime(double t0, double dt0) {
    t=t0;
    dt=dt0;
    TimestepDependence();
  }
};

Compare_t IntegratorCompare;
KeyCompare_t IntegratorKeyCompare;
extern IntegratorBase *Integrator;

#define INTEGRATOR(key)						\
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
    if(errmax > tolmax2) {
      return UNSUCCESSFUL;
    }
    if(errmax < tolmin2) {
      return ADJUST;
    }
  } else if(++dynamic == 0) dynamic=1;
  return SUCCESSFUL;
}

class Euler : public IntegratorBase {
public:
  const char *Name() {return "Euler";}
  Solve_RC Solve();
  void Predictor(unsigned int n0, unsigned int ny) {
    for(unsigned int j=n0; j < ny; j++) y[j] += dt*source[j];
  }
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
  void Allocator() {
  IntegratorBase::Allocator();
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
  AdamsBashforth() {order=3;}
  void Allocator() {
    IntegratorBase::Allocator();
    Alloc(Y0,y0);
    Alloc0(Src0,source0);
    Alloc0(Src1,source1);
  }
  const char *Name() {return "Third-Order Adams-Bashforth-Moulton";}
  Solve_RC Solve();
  void TimestepDependence() {
    init=2;
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
public:
  PC() {order=2;}
  void Allocator(bool base) {
    if(base) IntegratorBase::Allocator();
    Alloc(Y0,y0);
    Alloc0(Src0,source0);
    new_y0=1;
  }
  void Allocator() {
    Allocator(true);
  }
  const char *Name() {return "Predictor-Corrector";}
  Solve_RC Solve();
	
  void TimestepDependence() {
    halfdt=0.5*dt;
  }
  virtual void Predictor(unsigned int n0, unsigned int ny);
  virtual int Corrector(unsigned int n0, unsigned int ny);
};

class LeapFrog : public PC {
  vector yp,yp0;
  vector2 YP,YP0; // array of dependent fields
  double oldhalfdt,lasthalfdt;
public:
  void Allocator() {
  PC::Allocator(); 
  Alloc(YP,yp);
  Alloc(YP0,yp0);
  lasthalfdt=0.0;
  }
  const char *Name() {return "Leapfrog";}
  void TimestepDependence() {
    if(lasthalfdt == 0.0) {
      Set(yp,YP[0]);
      Set(y,Y[0]);
      for(unsigned int j=0; j < ny; j++) yp[j] = y[j];
    }
    halfdt=0.5*dt;
  }
  void Predictor(unsigned int n0, unsigned int ny);
  int Corrector(unsigned int n0, unsigned int ny);
};

class SYM2 : public PC {
public:
  const char *Name() {return "Second-Order Symplectic";}
  void Predictor(unsigned int n0, unsigned int ny);
  int Corrector(unsigned int n0, unsigned int ny);
};

class RK2 : public PC {
public:
  const char *Name() {return "Second-Order Runge-Kutta";}
  void Predictor(unsigned int n0, unsigned int ny);
  int Corrector(unsigned int n0, unsigned int ny);
};

class RK3 : public PC {
protected:  
  vector source1;
  vector2 Src1;
  double sixthdt;
public:
  RK3() {order=3;}
  void Allocator() {
    PC::Allocator();
    Alloc0(Src1,source1);
  }
  const char *Name() {return "Third-Order Runge-Kutta";}
  void TimestepDependence();
  void Predictor(unsigned int n0, unsigned int ny);
  int Corrector(unsigned int n0, unsigned int ny);
};

class RK4 : public RK3 {
protected:  
  vector source2;
  vector2 Src2;
public:
  RK4() {order=4;}
  void Allocator() {
    RK3::Allocator();
    Alloc0(Src2,source2);
  }
  const char *Name() {return "Fourth-Order Runge-Kutta";}
  void Predictor(unsigned int n0, unsigned int ny);
  int Corrector(unsigned int n0, unsigned int ny);
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
  void Allocator() {
    RK4::Allocator();
    Alloc(Src3,source3,Src1,source1);
    Alloc0(Src4,source4);
  }
  const char *Name() {return "Fifth-Order Runge-Kutta";}
  void TimestepDependence();
  void Predictor(unsigned int n0, unsigned int ny);
  int Corrector(unsigned int n0, unsigned int ny);
};

class Exact : public RK5 {
public:
  const char *Name() {return "Exact";}
  int Microfactor(){return 100;}
};

#endif
