// Allow for multiple B stages?

#ifndef __Integrator_h__
#define __Integrator_h__ 1

#include "parallel.h"

extern size_t threads;

class IntegratorBase {
protected:
  const char *abbrev;
  ProblemBase *Problem;
  size_t ny;
  vector y,source;
  vector2 Y,Src,YI,Yout; // arrays of dependent fields
  DynVector<size_t> NY; // number of variables in each field
  double t;
  double dt;
  double sample;
  double *errMax;
  double errmax;
  uvector errmask;
  double tolmax2,tolmin2;
  double stepfactor,stepinverse,stepnoninverse;
  double growfactor,shrinkfactor;
  double dtmin,dtmax;
  size_t itmax;
  size_t microsteps;
  Real microfactor;
  bool microprocess;
  size_t verbose;
  size_t dynamic;
  size_t order;
  double pgrow, pshrink;
  bool FSAL; // First Same As Last
  bool first;
  size_t align;

public:

  IntegratorBase(int order=0, bool fsal=false) : errMax(NULL), errmax(0.0),
                                                 order(order), FSAL(fsal) {}
  virtual ~IntegratorBase() {
    if(errMax)
      delete [] errMax;
  }

  void SetAbbrev(const char *abbrev0) {abbrev=abbrev0;}
  const char *Abbrev() {return abbrev;}
  void SetParam(double tolmax, double tolmin, double stepfactor0,
		double stepnoninverse0, double dtmin0, double dtmax0,
		size_t itmax0, size_t microsteps0, Real microfactor0,
		size_t verbose0, size_t dynamic0) {
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
    if(dynamic)
      errMax=new double[threads];
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
		 size_t& iteration, size_t& nout);
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
    double e=errMax[0];
    for(size_t t=1; t < threads; ++t)
      e=max(e,errMax[t]);
    return e;
  }

  size_t Order() {
    return order;
  }

  virtual bool fsal() {return FSAL;}

  virtual bool Exponential() {return false;}

  size_t Ny() {
    return ny;
  }

  // Add tolgood to inhibit time step adjustment.
  virtual void ExtrapolateTimestep () {
    errmax=Errmax();
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

  size_t Start(size_t field) {return Problem->Start(field);}
  size_t Stop(size_t field) {return Problem->Stop(field);}

  virtual void Allocator(const vector2& Y0, DynVector<size_t>* NY0,
			 const uvector& mask, size_t Align=0);
  virtual void Allocator() {}
  virtual void Allocator(ProblemBase& problem, size_t Align=0) {
    SetProblem(problem);
    align=Align;
    Allocator(problem.YVector(),problem.Sizes(),problem.ErrorMask(),align);
  }

  virtual const char *Name()=0;
  virtual Solve_RC Solve()=0;
  virtual size_t Microfactor() {return 1;}
  virtual void TimestepDependence() {}

  virtual void Unswap() {
    if(Yout.Size() && Yout != Y) set(Yout[0],Y[0],ny);
  }

  void SetTime(double t0, double dt0) {
    t=t0;
    dt=dt0;
  }

  void initError() {
    if(dynamic)
      for(size_t t=0; t < threads; ++t)
        errMax[t]=0.0;
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
  size_t thread=get_thread_num0();
  if(isfinite(norm0) && isfinite(diff)) {
    if(initial != 0.0 && pred != initial) {
      static const double epsilon=DBL_MIN/DBL_EPSILON;
      double error=max(abs2(diff)/(max(norm2(norm0),norm2(initial))+epsilon));
      errMax[thread]=max(errMax[thread],error);
    }
  } else errMax[thread]=HUGE_VAL;
}

inline Solve_RC IntegratorBase::CheckError()
{
  if(errmax > tolmax2) {
    return UNSUCCESSFUL;
  }
  if(errmax < tolmin2) {
    return ADJUST;
  }
  return SUCCESSFUL;
}

class Euler : public IntegratorBase {
public:
  const char *Name() {return "Euler";}
  Solve_RC Solve();
  virtual void Predictor(size_t n0, size_t ny) {
    PARALLELIF(
      ny > threshold,
      for(size_t j=n0; j < ny; j++)
        y[j] += dt*source[j];
      );
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
  size_t init;
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
  size_t init;
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
  PC(size_t order=2, bool fsal=false) : IntegratorBase(order,fsal) {}
  void Allocator() {
    Alloc(Y0,y0);
    Alloc0(Src0,source0);
    new_y0=true;
  }
  const char *Name() {return "Predictor-Corrector";}
  Solve_RC Solve();

  virtual bool isConservative() {return false;}

  void TimestepDependence() {
    halfdt=0.5*dt;
  }

  void iSource() {
    Source(Src0,Y0,t);
  }

  void initialize0() {
    initError();
    Set(y,Y[0]);
    Set(y0,Y0[0]);
    PARALLELIF(
      ny > threshold,
      for(size_t i=0; i < ny; ++i)
        y0[i]=y[i];
      );
  }

  void initialize() {
    initError();
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

  virtual void Predictor(size_t n0, size_t ny);
  virtual int Corrector(size_t n0, size_t ny);
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
      PARALLELIF(
        ny > threshold,
        for(size_t j=0; j < ny; j++)
          yp[j] = y[j];
        );
    }
    halfdt=0.5*dt;
  }
  void Predictor(size_t n0, size_t ny);
  int Corrector(size_t n0, size_t ny);
};

class SYM2 : public PC {
public:
  const char *Name() {return "Second-Order Symplectic";}
  void Predictor(size_t n0, size_t ny);
  int Corrector(size_t n0, size_t ny);
};

class RK : public PC {
protected:
  Array::array2<double> a,A;   // source coefficients for each stage
  Array::array1<double>::opt b,B,C; // b=error coefficients, C=step fractions
  const size_t nstages;
  vector3 vSrc;
  vector2 vsource;
  size_t Astages;
public:

  RK(size_t order, size_t nstages, bool fsal=false) :
    PC(order,fsal), nstages(fsal ? nstages : nstages-1), Astages(nstages-1) {}

  Var getvsource(size_t stage, size_t i) {
    return vsource[stage][i];
  }
  void setvsource(size_t stage, size_t i, Var value) {
    vsource[stage][i]=value;
  }

  size_t NStages() {return nstages;}

  virtual void Source(const vector2& Src, const vector2& Y, double t) {
    Problem->Source(Src,Y,t);
  }


  void Source(size_t i) {
    Source(vSrc[i],Y,t+C[i]*dt);
  }

  void Allocator(const vector2& YP,
		 DynVector<size_t>* NYP,
		 const uvector& errmask0,size_t Align=0) {
    IntegratorBase::Allocator(YP,NYP,errmask0,Align);
  }

  void Csum() {
    for(size_t s=0; s < Astages; ++s) {
      Real sum=0.0;
      rvector As=A[s];
      for(size_t k=0; k <= s; ++k)
	sum += As[k];
      C[s]=sum;
    }
    for(size_t k=0; k < nstages; k++)
      if(B[k] != 0.0) return;
  }

  void Allocator() {
    Alloc(Y0,y0);
    vSrc.Allocate(nstages);
    vsource.Allocate(nstages);
    for(size_t s=0; s < nstages; s++)
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

  virtual void Stage(size_t s, size_t start, size_t stop) {
    rvector as=a[s];
    PARALLELIF(
      (stop-start)*(s+1) > threshold,
      for(size_t j=start; j < stop; j++) {
        Var sum=y0[j];
        for(size_t k=0; k <= s; k++)
          sum += as[k]*vsource[k][j];
        y[j]=sum;
      });
  }

  void Stage(size_t s, size_t start=0) {
    Stage(s,start,ny);
  }

  virtual void PStage(size_t s) {
    Stage(s);
  }

  void PredictorSource(size_t s) {
    double cs=C[s]*dt;
    Problem->BackTransform(Y,t+cs,cs,YI);
    Source(vSrc[s+1],Y,t+cs);
    if(Array::Active(YI)) {swaparray(YI,Y); Set(y,Y[0]);}
  }

  virtual void Predictor(size_t start, size_t stop) {
    size_t laststage=Astages-1;
    for(size_t s=0; s < laststage; ++s) {
      Stage(s,start,stop);
      PredictorSource(s);
    }
  }

  virtual int Corrector(size_t start, size_t stop);

  void allocate() {
    A.Allocate(Astages,Astages);
    a.Allocate(Astages,Astages);
    A=0.0;

    Allocate(B,nstages);
    Allocate(b,nstages);
    for(size_t i=0; i < nstages; ++i)
      B[i]=0.0;

    Allocate(C,Astages);
  }

  virtual void TimestepDependence() {
    for(size_t s=0; s < Astages; ++s) {
      rvector as=a[s];
      rvector As=A[s];
      for(size_t k=0; k <= s; k++)
	as[k]=dt*As[k];
    }
    for(size_t k=0; k < nstages; k++)
      b[k]=dt*B[k];
  }
};

class RK1 : public RK {
public:
  const char *Name() {return "First-Order Runge-Kutta";}
  RK1() : RK(1,2) {
    allocate();
    A[0][0]=1.0;
  }
};

class RK2 : public RK {
public:
  const char *Name() {return "Second-Order Runge-Kutta";}

  RK2() : RK(2,3) {
    allocate();
    A[0][0]=0.5;

    A[1][1]=1.0;

    B[0]=1.0;
  }
};

class RKPC : public RK {
public:
  const char *Name() {
    return "Two-Stage Second-Order Heun Predictor-Corrector";
  }

  RKPC() : RK(2,3) {
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

class RK3K : public RK {
public:
  const char *Name() {return "Third-Order Kutta";}

  RK3K() : RK(3,4) {
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
  const char *Name() {return "Six-Stage Fourth-Order Runge-Kutta";}

  void Predictor(size_t start, size_t stop) {
    Real a00=a[0][0];
    PARALLELIF(
      stop-start > threshold,
      for(size_t j=start; j < stop; j++)
        y[j]=y0[j]+a00*vsource[0][j];
      );
    PredictorSource(0);

    Real a11=a[1][1];
    PARALLELIF(
      stop-start > threshold,
      for(size_t j=start; j < stop; j++)
        y[j]=y0[j]+a11*vsource[1][j];
      );
    PredictorSource(1);

    Real a22=a[2][2];
    PARALLELIF(
      stop-start > threshold,
      for(size_t j=start; j < stop; j++)
        y[j]=y0[j]+a22*vsource[2][j];
      );
    PredictorSource(2);

    if(dynamic) {
      rvector a3=a[3];
      Real a30=a3[0];
      Real a31=a3[1];
      PARALLELIF(
        stop-start > threshold,
        for(size_t j=start; j < stop; j++)
          y[j]=y0[j]+a30*vsource[0][j]+a31*vsource[1][j];
        );
      PredictorSource(3);
    }
  }

  int Corrector(size_t start, size_t stop) {
    if(dynamic) {
      rvector as=a[4];
      PARALLELIF(
        stop-start > threshold,
        for(size_t j=start; j < stop; j++) {
          Var sum0=y0[j];
          Var sum=sum0+as[0]*vsource[0][j]+as[1]*vsource[1][j]+
	  as[2]*vsource[2][j]+as[3]*vsource[3][j];
          Var pred=sum0+b[0]*vsource[0][j]+b[1]*vsource[1][j]+
	  b[4]*vsource[4][j];
          if(!Array::Active(errmask) || errmask[j])
            CalcError(sum0,sum,pred,sum);
          y[j]=sum;
        });
    } else {
      rvector as=a[4];
      PARALLELIF(
        stop-start > threshold,
        for(size_t j=start; j < stop; j++)
          y[j]=y0[j]+as[0]*vsource[0][j]+as[1]*vsource[1][j]+
            as[2]*vsource[2][j]+as[3]*vsource[3][j];
        );
    }
    return 1;
  };

  RK4() : RK(4,6) {
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

class RK43 : public RK {
public:
  const char *Name() {return "Robust Five-Stage Fourth-Order Runge-Kutta";}

  RK43() : RK(4,5,true) {
    allocate();
    A[0][0]=0.125;

    A[1][0]=-1.0/27.0;
    A[1][1]=10.0/27.0;

    A[2][0]=73.0/17.0;
    A[2][1]=-658.0/85.0;
    A[2][2]=378.0/85.0;

    A[3][0]=2.0/3.0;
    A[3][1]=-128.0/105.0;
    A[3][2]=1.35;
    A[3][3]=17.0/84.0;

    B[0]=11.0/3.0;
    B[1]= -704.0/105.0;
    B[2]=81.0/20.0;
    B[3]=-85.0/84.0;
    B[4]=1.0;
  }
};

class RK5 : public RK {
public:
  const char *Name() {return "Fifth-Order Cash-Karp Runge-Kutta";}

  int Corrector(size_t start, size_t stop) {
    if(dynamic) {
      rvector as=a[5];
      Real a0=as[0];
      Real a2=as[2];
      Real a3=as[3];
      Real a5=as[5];
      Real b0=b[0];
      Real b2=b[2];
      Real b3=b[3];
      Real b4=b[4];
      Real b5=b[5];
      vector vsource0=vsource[0];
      vector vsource2=vsource[2];
      vector vsource3=vsource[3];
      vector vsource4=vsource[4];
      vector vsource5=vsource[5];
      PARALLELIF(
        stop-start > threshold,
        for(size_t j=start; j < stop; j++) {
          Var sum0=y0[j];
          Var v0=vsource0[j];
          Var v2=vsource2[j];
          Var v3=vsource3[j];
          Var v5=vsource5[j];
          Var sum=sum0+a0*v0+a2*v2+a3*v3+a5*v5;
          Var pred=sum0+b0*v0+b2*v2+b3*v3+b4*vsource4[j]+b5*v5;
          if(!Array::Active(errmask) || errmask[j])
            CalcError(sum0,sum,pred,sum);
          y[j]=sum;
        });
    } else {
      rvector as=a[5];
      Real a0=as[0];
      Real a2=as[2];
      Real a3=as[3];
      Real a5=as[5];
      vector vsource0=vsource[0];
      vector vsource2=vsource[2];
      vector vsource3=vsource[3];
      vector vsource5=vsource[5];
      PARALLELIF(
        stop-start > threshold,
        for(size_t j=start; j < stop; j++)
          y[j]=y0[j]+a0*vsource0[j]+a2*vsource2[j]+
            a3*vsource3[j]+a5*vsource5[j];
        );
    }
    return 1;
  };

  RK5() : RK(5,7) {
    allocate();

    A[0][0]=0.2;

    A[1][0]=3.0/40.0;
    A[1][1]=9.0/40.0;

    A[2][0]=0.3;
    A[2][1]=-0.9;
    A[2][2]=1.2;

    A[3][0]=-11.0/54.0;
    A[3][1]=2.5;
    A[3][2]=-70.0/27.0;
    A[3][3]=35.0/27.0;

    A[4][0]=1631.0/55296.0;
    A[4][1]=175.0/512.0;
    A[4][2]=575.0/13824.0;
    A[4][3]=44275.0/110592.0;
    A[4][4]=253.0/4096.0;

    A[5][0]=37.0/378.0;
    A[5][2]=250.0/621.0;
    A[5][3]=125.0/594.0;
    A[5][5]=512.0/1771.0;

    B[0]=2825.0/27648.0;
    B[2]=18575.0/48384.0;
    B[3]=13525.0/55296.0;
    B[4]=277.0/14336.0;
    B[5]=0.25;
  }
};

class RK5DP : public RK {
public:
  const char *Name() {return "Fifth-Order Dormand-Prince Runge-Kutta";}

  RK5DP() : RK(5,7,true) {
    allocate();

    A[0][0]=0.2;

    A[1][0]=3.0/40.0;
    A[1][1]=9.0/40.0;

    A[2][0]=44.0/45.0;
    A[2][1]=-56.0/15.0;
    A[2][2]=32.0/9.0;

    A[3][0]=19372.0/6561.0;
    A[3][1]=-25360.0/2187.0;
    A[3][2]=64448.0/6561.0;
    A[3][3]=-212.0/729.0;

    A[4][0]=9017.0/3168.0;
    A[4][1]=-355.0/33.0;
    A[4][2]=46732.0/5247.0;
    A[4][3]=49.0/176.0;
    A[4][4]=-5103.0/18656.0;

    A[5][0]=35.0/384.0;
    A[5][2]=500.0/1113.0;
    A[5][3]=125.0/192.0;
    A[5][4]=-2187/6784.0;
    A[5][5]=11.0/84.0;

    B[0]=5179.0/57600.0;
    B[2]=7571.0/16695.0;
    B[3]=393.0/640.0;
    B[4]=-92097.0/339200.0;
    B[5]=187.0/2100.0;
    B[6]=1.0/40.0;
  }
};

class Exact : public RK5 {
public:
  const char *Name() {return "Exact";}
  size_t Microfactor(){return 100;}
};

#endif
