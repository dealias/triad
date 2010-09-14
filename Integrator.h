#ifndef __Integrator_h__
#define __Integrator_h__ 1

class IntegratorBase {
protected:
  const char *abbrev;
  ProblemBase *Problem;
  unsigned int ny;
  vector y,source;
  vector2 Y,Src,YI,Yout; // arrays of dependent fields
  DynVector<unsigned int> NY; // number of variables in each field
  double t;
  double dt;
  double sample;
  double errmax;
  ivector errmask;
  double tolmax2,tolmin2;
  double stepfactor,stepinverse,stepnoninverse;
  double growfactor,shrinkfactor;
  double dtmin,dtmax;
  long long itmax;
  int microsteps;
  Real microfactor;
  int microprocess;
  int verbose;
  int dynamic;
  int order;
  double pgrow, pshrink;
  bool FSAL; // First Same As Last
  bool first;
  size_t align;
  
public:	
  
  IntegratorBase(int order=0, bool fsal=false) : order(order), FSAL(fsal) {}
  
  virtual ~IntegratorBase() {}
  void SetAbbrev(const char *abbrev0) {abbrev=abbrev0;}
  const char *Abbrev() {return abbrev;}
  void SetParam(double tolmax, double tolmin, double stepfactor0,
		double stepnoninverse0, double dtmin0, double dtmax0,
		long long itmax0, int microsteps0, Real microfactor0,
		int verbose0, int dynamic0) {
    if(tolmax < tolmin) msg(ERROR_GLOBAL,"tolmax < tolmin"); 
    tolmax2=tolmax*tolmax;
    tolmin2=tolmin*tolmin;
    growfactor=stepfactor=stepfactor0;
    shrinkfactor=stepinverse=1.0/stepfactor;
    stepnoninverse=1.0/stepnoninverse0;
    dtmin=dtmin0;
    dtmax=dtmax0;
    itmax=itmax0;
    microsteps=microsteps0*Microfactor();
    microfactor=microfactor0;
    verbose=verbose0;
    dynamic=dynamic0;
  }
  void SetParam(const IntegratorBase& I) {
    tolmax2=I.tolmax2;
    tolmin2=I.tolmin2;
    stepfactor=I.stepfactor;
    stepnoninverse=I.stepnoninverse;
    dtmin=I.dtmin;
    dtmax=I.dtmax;
    itmax=I.itmax;
    microsteps=I.microsteps;
    microfactor=I.microfactor;
    verbose=I.verbose;
    dynamic=I.dynamic;
  }
  void Integrate(double& t0, double tmax, double& dt0, double sample0,
		 long long& iteration, unsigned long& nout);
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
  virtual Solve_RC CheckError();
  
  double Errmax() {
    return errmax;
  }
  
  int Order() {
    return order;
  }

  virtual bool fsal() {return FSAL;}

  unsigned int Ny() {
    return ny;
  }
  
  // Add tolgood to inhibit time step adjustment.
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
  
  virtual void Allocator(const vector2& Y0, DynVector<unsigned int>* NY0,
			 const ivector& mask, size_t Align=0);
  virtual void Allocator() {}
  virtual void Allocator(ProblemBase& problem, size_t Align=0) {
    SetProblem(problem);
    align=Align;
    Allocator(problem.YVector(),problem.Sizes(),problem.ErrorMask(),align);
  }
  
  virtual const char *Name()=0;
  virtual Solve_RC Solve()=0;
  virtual int Microfactor() {return 1;}
  virtual void TimestepDependence() {}
  
  virtual void Unswap() {
    if(Yout.Size() && Yout != Y) set(Yout[0],Y[0],ny);
  }
  
  void SetTime(double t0, double dt0) {
    t=t0;
    dt=dt0;
  }
  
  void SetTime(double t0, double dt0, double errmax0) {
    SetTime(t0,dt0);
    errmax=errmax0;
  }
  
  const vector2& YVector() const {
    return Y;
  }
  
  void Sync(IntegratorBase *I) {
    Set(Y,I->YVector());
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
  Var diff=corr-pred;
  if(isfinite(norm0) && isfinite(diff)) {
    if(initial != 0.0 && pred != initial) {
      static const double epsilon=DBL_MIN/DBL_EPSILON;
      double error=max(abs2(diff)/(max(norm2(norm0),norm2(initial))+epsilon));
      if(error > errmax) errmax=error;
    }
  } else errmax=HUGE_VAL;
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
  virtual void Predictor(unsigned int n0, unsigned int ny) {
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
    Alloc(Y0,y0);
    halfdt=0.5*dt;
  }
  const char *Name() {return "Midpoint Rule";}
  void TimestepDependence() {
    halfdt=0.5*dt;
  }
  Solve_RC Solve();
};

class AB2 : public IntegratorBase {
  vector y0,source0,source1;
  vector2 Y0,Src0,Src1;
  double a0,a1;
  int init;
public:
  AB2() : IntegratorBase(2) {}
  void Allocator() {
    Alloc(Y0,y0);
    Alloc0(Src0,source0);
    Alloc0(Src1,source1);
  }
  const char *Name() {return "Second-Order Adams-Bashforth";}
  Solve_RC Solve();
  void TimestepDependence() {
    init=1;
    a0=1.5*dt;
    a1=-0.5*dt;
  }
};

class ABM3 : public IntegratorBase {
  vector y0,source0,source1;
  vector2 Y0,Src0,Src1;
  double a0,a1,a2;
  double b0,b1,b2;
  int init;
public:
  ABM3() : IntegratorBase(3) {}
  void Allocator() {
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
  bool new_y0;
  vector y0,source0;
  vector2 Y0,Src0;
  double halfdt;
public:
  PC(int order=2, bool fsal=false) : IntegratorBase(order,fsal) {}
  void Allocator() {
    Alloc(Y0,y0);
    Alloc0(Src0,source0);
    new_y0=true;
  }
  const char *Name() {return "Predictor-Corrector";}
  Solve_RC Solve();
	
  void TimestepDependence() {
    halfdt=0.5*dt;
  }
  
  void initialize0() {
    errmax=0.0;
    Set(y,Y[0]);
    Set(y0,Y0[0]);
    for(unsigned int i=0; i < ny; ++i)
      y0[i]=y[i];
  }
  
  void iSource() {
    Source(Src0,Y0,t);
  }

  void initialize() {
    errmax=0.0;
    if(new_y0) {
      swaparray(Y0,Y);
      Set(y,Y[0]);
      Set(y0,Y0[0]);
      if(first || !fsal()) {
	iSource();
	first=false;
      }
    }
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

class RK : public PC {
protected:  
  Array::array2<double> a,A;   // source coefficients for each stage
  Array::array1<double>::opt b,B,C; // b=error coefficients, C=time coefficients
  const unsigned int nstages;
  vector3 vSrc;
  vector2 vsource;
  unsigned int Astages;
public:

  RK(int order, int nstages, bool fsal=false) :
    PC(order,fsal), nstages(nstages), Astages(fsal ? nstages-1 : nstages) {}
    
  Var getvsource(unsigned int stage, unsigned int i) {
    return vsource[stage][i];
  }
  void setvsource(unsigned int stage, unsigned int i, Var value) {
    vsource[stage][i]=value;
  }

  unsigned int NStages() {return nstages;}
    
  void setnew_y0(bool flag) {new_y0=flag;}

  virtual void Source(const vector2& Src, const vector2& Y, double t) {
    Problem->Source(Src,Y,t);
  }
  
  void Source(unsigned int i) {
    Source(vSrc[i],Y,t+C[i]*dt);
  }
  
  void Allocator(const vector2& Y0, 
		 DynVector<unsigned int>* NY0,
		 const ivector& errmask0,size_t Align=0) {
    IntegratorBase::Allocator(Y0,NY0,errmask0,Align);
  }

  void Csum() {
    for(unsigned int s=0; s < Astages; ++s) {
      Real sum=0.0;
      rvector As=A[s];
      for(unsigned int k=0; k <= s; ++k)
	sum += As[k];
      C[s]=sum;
    }
    
    for(unsigned int k=0; k < nstages; k++)
      if(B[k] != 0.0) return;
    dynamic=0;
  }
  
  void Allocator() {
    Alloc(Y0,y0);
    vSrc.Allocate(nstages);
    vsource.Allocate(nstages);
    for(unsigned int s=0; s < nstages; s++)
      Alloc0(vSrc[s],vsource[s]);
    Src0.Dimension(vSrc[0]);
    new_y0=true;
    Csum();
  }
  
  bool fsal() {
    if(dynamic && FSAL) {
      swaparray(vSrc[0],vSrc[Astages]);
      Set(vsource[0],vSrc[0][0]);
      Set(vsource[Astages],vSrc[Astages][0]);
      return true;
    }
    return false;
  }
  
  virtual void Stage(unsigned int s, unsigned int start, unsigned int stop) {
    rvector as=a[s];
    for(unsigned int j=start; j < stop; j++) {
      Var sum=y0[j];
      for(unsigned int k=0; k <= s; k++)
	sum += as[k]*vsource[k][j];
      y[j]=sum;
    }
  }
  
  void Stage(unsigned int s, int start=0) {
    Stage(s,start,ny);
  }
  
  virtual void PStage(unsigned int s) {
    Stage(s);
  }

  void PredictorSource(unsigned int s) {
    double cs=C[s]*dt;
    Problem->BackTransform(Y,t+cs,cs,YI);
    Source(vSrc[s+1],Y,t+cs);
    if(Array::Active(YI)) {swaparray(YI,Y); Set(y,Y[0]);}
  }
  
  virtual void Predictor(unsigned int start, unsigned int stop) {
    for(unsigned int s=0; s < Astages-1; ++s) {
      Stage(s,start,stop);
      PredictorSource(s);
    }
  }
  
  virtual int Corrector(unsigned int start, unsigned int stop);
  
  void allocate() {
    A.Allocate(Astages,Astages);
    a.Allocate(Astages,Astages);
    A=0.0;
    
    Allocate(B,nstages);
    Allocate(b,nstages);
    for(unsigned int i=0; i < nstages; ++i)
      B[i]=0.0;
    
    Allocate(C,Astages);
  }
    
  virtual void TimestepDependence() {
    for(unsigned int s=0; s < Astages; ++s) {
      rvector as=a[s];
      rvector As=A[s];
      for(unsigned int k=0; k <= s; k++)
	as[k]=dt*As[k];
    }
    for(unsigned int k=0; k < nstages; k++)
      b[k]=dt*B[k];
  }
};
  
class RK1 : public RK {
public:
  const char *Name() {return "First-Order Runge-Kutta";}
  RK1() : RK(1,1) {
    allocate();
    A[0][0]=1.0;
  }
};

class RK2 : public RK {
public:
  const char *Name() {return "Second-Order Runge-Kutta";}
  
  RK2() : RK(2,2) {
    allocate();
    A[0][0]=0.5;
    
    A[1][1]=1.0;
    
    B[0]=1.0;
  }
};

class RKPC : public RK {
public:
  const char *Name() {return "Predictor-Corrector";}
  
  RKPC() : RK(2,2) {
    allocate();
    
    A[0][0]=1.0;
    
    A[1][0]=0.5;
    A[1][1]=0.5;
    
    B[0]=1.0;
  }
};

class RK3 : public RK {
public:
  const char *Name() {return "Third-Order Bogacki-Shampine Runge-Kutta";}
  
  RK3() : RK(3,4,true) {
    allocate();
    A[0][0]=0.5;
    
    A[1][1]=0.75;
    
    A[2][0]=2.0/9.0;
    A[2][1]=1.0/3.0;
    A[2][2]=4.0/9.0;
    
    B[0]=7.0/24.0;
    B[1]=0.25;
    B[2]=1.0/3.0;
    B[3]=0.125;
  }
};

class RK3C : public RK {
public:
  const char *Name() {return "Third-Order Classical Runge-Kutta";}
  
  RK3C() : RK(3,3) {
    allocate();
    A[0][0]=0.5;
    
    A[1][0]=-1.0; A[1][1]=2.0;
    
    A[2][0]=1.0/6.0;
    A[2][1]=2.0/3.0;
    A[2][2]=1.0/6.0;
    
    B[1]=1.0;
  }
};

class RK4 : public RK {
public:
  const char *Name() {return "Fourth-Order Runge-Kutta";}

  void Predictor(unsigned int start, unsigned int stop) {
    Real a00=a[0][0];
    for(unsigned int j=start; j < stop; j++)
      y[j]=y0[j]+a00*vsource[0][j];
    PredictorSource(0);
    
    Real a11=a[1][1];
    for(unsigned int j=start; j < stop; j++)
      y[j]=y0[j]+a11*vsource[1][j];
    PredictorSource(1);
    
    Real a22=a[2][2];
    for(unsigned int j=start; j < stop; j++)
      y[j]=y0[j]+a22*vsource[2][j];
    PredictorSource(2);
    
    if(dynamic) {
      rvector a3=a[3];
      Real a30=a3[0];
      Real a31=a3[1];
      for(unsigned int j=start; j < stop; j++)
	y[j]=y0[j]+a30*vsource[0][j]+a31*vsource[1][j];
      PredictorSource(3);
    }
  }

  int Corrector(unsigned int start, unsigned int stop) {
    if(dynamic) {
      rvector as=a[4];
      for(unsigned int j=start; j < stop; j++) {
	Var sum0=y0[j];
	Var sum=sum0+as[0]*vsource[0][j]+as[1]*vsource[1][j]+
	  as[2]*vsource[2][j]+as[3]*vsource[3][j];
	Var pred=sum0+b[0]*vsource[0][j]+b[1]*vsource[1][j]+
	  b[4]*vsource[4][j];
	if(!Array::Active(errmask) || errmask[j])
	  CalcError(sum0,sum,pred,sum);
	y[j]=sum;
      }
    } else {
      rvector as=a[4];
      for(unsigned int j=start; j < stop; j++) {
	y[j]=y0[j]+as[0]*vsource[0][j]+as[1]*vsource[1][j]+
	  as[2]*vsource[2][j]+as[3]*vsource[3][j];
      }
    }
    return 1;
  };
  
  RK4() : RK(4,5) {
    allocate();
    A[0][0]=0.5;

    A[1][1]=0.5;
    
    A[2][2]=1.0;
    
    A[3][0]=-1.0;
    A[3][1]=2.0;

    A[4][0]=1.0/6.0;
    A[4][1]=1.0/3.0;
    A[4][2]=1.0/3.0;
    A[4][3]=1.0/6.0;

    B[0]=1.0/6.0;
    B[1]=2.0/3.0;
    B[4]=1.0/6.0;

  }
};

class RK5 : public RK {
public:
  const char *Name() {return "Fifth-Order Runge-Kutta";}
  
  int Corrector(unsigned int start, unsigned int stop) {
    if(dynamic) {
      rvector as=a[5];
      for(unsigned int j=start; j < stop; j++) {
	Var sum0=y0[j];
	Var sum=sum0+as[0]*vsource[0][j]+as[2]*vsource[2][j]+
	  as[3]*vsource[3][j]+as[5]*vsource[5][j];
	Var pred=sum0+b[0]*vsource[0][j]+b[2]*vsource[2][j]+
	  +b[3]*vsource[3][j]+b[4]*vsource[4][j]+b[5]*vsource[5][j];
	if(!Array::Active(errmask) || errmask[j])
	  CalcError(sum0,sum,pred,sum);
	y[j]=sum;
      }
    } else {
      rvector as=a[5];
      for(unsigned int j=start; j < stop; j++) {
	y[j]=y0[j]+as[0]*vsource[0][j]+as[2]*vsource[2][j]+
	  as[3]*vsource[3][j]+as[5]*vsource[5][j];
      }
    }
    return 1;
  };
  
  RK5() : RK(5,6) {
    allocate();
    
    A[0][0]=0.2;
    A[1][0]=3.0/40.0; A[1][1]=9.0/40.0;
    
    A[2][0]=0.3; A[2][1]=-0.9; A[2][2]=1.2;
    
    A[3][0]=-11.0/54.0; A[3][1]=2.5; A[3][2]=-70.0/27.0; A[3][3]=35.0/27.0;
  
    A[4][0]=1631.0/55296.0; A[4][1]=175.0/512.0; A[4][2]=575.0/13824.0;
    A[4][3]=44275.0/110592.0; A[4][4]=253.0/4096.0;
  
    A[5][0]=37.0/378.0; A[5][2]=250.0/621.0; A[5][3]=125.0/594.0;
    A[5][5]=512.0/1771.0;
								  
    B[0]=2825.0/27648.0; B[2]=18575.0/48384.0; B[3]=13525.0/55296.0;
    B[4]=277.0/14336.0; B[5]=0.25;
  }
};

class Exact : public RK5 {
public:
  const char *Name() {return "Exact";}
  int Microfactor(){return 100;}
};

#endif
