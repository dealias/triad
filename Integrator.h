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
		double stepnoninverse0, double dtmin0, double dtmax0,
		int itmax0, int microsteps0, int verbose0, int dynamic0) {
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
    verbose=I.verbose;
    dynamic=I.dynamic;
  }
  void Integrate(double& t0, double tmax, double& dt0, double sample0,
		 int& iteration, unsigned long& nout);
  void ChangeTimestep(double dtnew);
	
  virtual void Source(const vector2& Src, const vector2& Y, double t) {
    Problem->Source(Src,Y,t);
  }
	
  virtual bool fsal() {return false;} // First Same As Last
    
  virtual void PSource(const vector2& Src, const vector2& Y, double t) {
    Source(Src,Y,t);
  }
	
  virtual void CSource(const vector2& Src, const vector2& Y, double t) {
    Source(Src,Y,t);
  }
  
  void CalcError(const Var& initial, const Var& norm, const Var& pred,
		 const Var& corr);
  Solve_RC CheckError();
  
  double Errmax() {
    return errmax;
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
  
  void Allocator(const vector2& Y0, DynVector<unsigned int>* NY0,
		 const ivector& mask);
  virtual void Allocator() {}
  virtual void Allocator(ProblemBase& problem) {
    SetProblem(problem);
    Allocator(problem.YVector(),problem.Sizes(),problem.ErrorMask());
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
    TimestepDependence();
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
  AB2() {order=2;}
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
  ABM3() {order=3;}
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
  PC() {order=2;}
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
  Array::array2<double> a,A;
  Array::array1<double> b,B,c;
  const unsigned int nstages;
  vector3 vSrc;
  vector2 vsource;
  bool FSAL;
public:
  RK(int nstages) : nstages(nstages) {}
  
  const unsigned int NStages() {return nstages;}
  const void Showy() {
    for (unsigned int i=0; i < ny; i++)
      cout << y[i] << ",";
    cout << endl;
  }
  void setnew_y0(bool flag) {new_y0=flag;}

  virtual void Source(const vector2& Src, const vector2& Y, double t) {
    Problem->Source(Src,Y,t);
  }
  
  void Source(unsigned int i) {
    //use vSrc[?] instead of vsource?
    Source(vSrc[i],Y,t+c[i]);
  }
  
  void Allocator(const vector2& Y0, 
		 DynVector<unsigned int>* NY0,
		 const ivector& errmask0) {
    IntegratorBase::Allocator(Y0,NY0,errmask0);
  }

  void Allocator() {
    Alloc(Y0,y0);
    vSrc.Allocate(nstages);
    vsource.Allocate(nstages);
    for(unsigned int s=0; s < nstages; s++)
      Alloc0(vSrc[s],vsource[s]);
    Src0.Dimension(vSrc[0]);
    new_y0=true;
  }
  
  void Stage(unsigned int s, unsigned int start=0, unsigned int stop) {
    rvector as=a[s];
    for(unsigned int j=start; j < stop; j++) {
      Var sum=y0[j];
      for(unsigned int k=0; k <= s; k++)
	sum += as[k]*vsource[k][j];
      y[j]=sum;
    }
  }
  
  void Stage(unsigned int s, unsigned int start=0) {
    Stage(s,0,ny);
  }
  
  void Predictor(unsigned int start, unsigned int stop) {
    for(unsigned int s=0; s < nstages-1; ++s) {
      Stage(s,start,stop);
      double cs=c[s];
      Problem->BackTransform(Y,t+cs,cs,YI);
      Source(vSrc[s+1],Y,t+cs);
      if(Array::Active(YI)) {swaparray(YI,Y); Set(y,Y[0]);}
    }
  }
  
  int Corrector(unsigned int start, unsigned int stop) {
    if(FSAL) msg(ERROR,"Too bad for you!");
    if(dynamic) {
      rvector as=a[nstages-1];
      for(unsigned int j=start; j < stop; j++) {
	Var sum0=y0[j];
	Var sum=sum0;
	Var pred=sum0;
	for(unsigned int k=0; k < nstages; k++) {
	  Var Skj=vsource[k][j];
	  sum += as[k]*Skj;
	  pred += b[k]*Skj;
	}
	if(!Array::Active(errmask) || errmask[j])
	  CalcError(sum0,sum,pred,sum);
	y[j]=sum;
      }
    } else Stage(nstages-1,start,stop);
    return 1;
  };
  
  void allocate(bool fsal=false) {
    FSAL=fsal;
    A.Allocate(nstages,nstages);
    a.Allocate(nstages,nstages);
    if(!fsal) {
      Allocate(B,nstages);
      Allocate(b,nstages);
    }
    Allocate(c,nstages);
    a=0.0;
  }
    
  void csum() {
    for(unsigned int s=0; s < nstages; ++s) {
      Real sum=0.0;
      for(unsigned int k=0; k <= s; ++k)
	sum += a[s][k];
      c[s]=sum;
    }
  }
  
  void TimestepDependence() {
    for(unsigned int s=0; s < nstages; ++s) {
      rvector as=a[s];
      rvector As=A[s];
      for(unsigned int k=0; k <= s; k++)
	as[k]=dt*As[k];
    }
    if(!FSAL)
      for(unsigned int k=0; k < nstages; k++)
	b[k]=dt*B[k];
    csum();
  }
};
  
class RK2p : public RK {
public:
  const char *Name() {return "Second-Order TEST Runge-Kutta";}
  
  RK2p() : RK(2) {
    allocate();
    A[0][0]=0.5;
    A[1][1]=1.0;
    
    B[0]=1.0;
    B[1]=0.0;
  }
};

class RK5p : public RK {
public:
  const char *Name() {return "Fifth-Order TEST Runge-Kutta";}
  
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
  
  RK5p() : RK(6) {
    allocate();
    order=5;
    
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


class RK2 : public PC {
public:
  const char *Name() {return "Second-Order Runge-Kutta";}
  void Predictor(unsigned int n0, unsigned int ny);
  int Corrector(unsigned int n0, unsigned int ny);
};

class RK3 : public PC {
protected:  
  vector source1,source2,source3;
  vector2 Src1,Src2,Src3;
  double a21;
  double a32;
  double b1,b2,b3;
  double B1,B2,B3,B4;
  double threefourthsdt;
public:
  RK3() {order=3;}
  const char *Name() {return "Third-Order Bogacki-Shampine Runge-Kutta";}
  
  void Allocator() {
    PC::Allocator();
    Alloc0(Src1,source1);
    Alloc0(Src2,source2);
    if(dynamic) Alloc0(Src3,source3);
  }
  
  void TimestepDependence();
  
  bool fsal() {
    if(dynamic) {
      swaparray(Src0,Src3);
      Set(source0,Src0[0]);
      Set(source3,Src3[0]);
      return true;
    }
    return false;
  }
  
  void Predictor(unsigned int n0, unsigned int ny);
  int Corrector(unsigned int n0, unsigned int ny);
};

class RK3C : public PC {
protected:  
  vector source1;
  vector2 Src1;
  double sixthdt;
public:
  RK3C() {order=3;}
  void Allocator() {
    PC::Allocator();
    Alloc0(Src1,source1);
  }
  const char *Name() {return "Third-Order Classical Runge-Kutta";}
  void TimestepDependence();
  void Predictor(unsigned int n0, unsigned int ny);
  int Corrector(unsigned int n0, unsigned int ny);
};

class RK4 : public RK3C {
protected:  
  vector source2,source3;
  vector2 Src2,Src3;
public:
  RK4() {order=4;}
  void Allocator() {
    RK3C::Allocator();
    Alloc0(Src2,source2);
    if(dynamic) Alloc0(Src3,source3);
  }
  const char *Name() {return "Fourth-Order Runge-Kutta";}
  void Predictor(unsigned int n0, unsigned int ny);
  int Corrector(unsigned int n0, unsigned int ny);
};

class RK5 : public RK4 {
protected:	
  vector source4;
  vector2 Src4;
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
    RK3C::Allocator();
    Alloc0(Src2,source2);
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
