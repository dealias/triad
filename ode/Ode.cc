#include "options.h"
#include "kernel.h"

using namespace Array;

const char *method="Ode";
const char *integrator="PC";

typedef array1<Nu>::opt Nuvector;

// Global variables
Var A;
Var B;
Nuvector nu;

// Local vocabulary declarations and default values
static Nu nu0=1.0;

class OdeVocabulary : public VocabularyBase {
public:
  const char *Name() {return "Ode";}
  const char *Abbrev() {return "Ode";}
  OdeVocabulary();
};

class Ode : public ProblemBase {
public:
  Ode();
	
  void InitialConditions();
  void Initialize();
  void Output(int it);
	
  void NonLinearSource(const vector2& Src, const vector2& Y, double t);
  void LinearSource(const vector2& Src, const vector2& Y, double t);
  
  void Source(const vector2& Src, const vector2& Y, double t) {
    NonLinearSource(Src,Y,t);
    LinearSource(Src,Y,t);
  }
};

Ode *OdeProblem;

Ode::Ode()
{
  OdeProblem=this;
}

class Exponential {
protected:
  Nuvector nu_inv, expinv, onemexpinv;
public:
  void Allocate(unsigned int n) {
    Allocate1(expinv,n);
    Allocate1(onemexpinv,n);
    Allocate1(nu_inv,n);
    
    for(unsigned int j=0; j < n; j++) {
      if(nu[j] != 0.0) nu_inv[j]=1.0/nu[j];
    }
  }

  void TimestepDependence(double dt, unsigned int n) {
    for(unsigned int j=0; j < n; j++) {
      Nu temp=-expm1(-nu[j]*dt);
      expinv[j]=1.0-temp;
      if(nu[j] == 0.0) onemexpinv[j]=dt;
      else onemexpinv[j]=temp*nu_inv[j];
    }
  }
};

class IntegratingFactor {
protected:
  Nuvector expinv;
public:
  void Allocate(unsigned int n) {Allocate1(expinv,n);}
  void TimestepDependence(double dt, unsigned int n) {
    for(unsigned int j=0; j < n; j++) expinv[j]=exp(-nu[j]*dt);
  }
};
	
class Implicit : public IntegratorBase {
public:
  void Allocate() {IntegratorBase::Allocate();}
  const char *Name() {return "Implicit";}
  Solve_RC Solve();
};

class E_Euler : public Euler, public Exponential {
public:
  void Allocate() {Euler::Allocate(); Exponential::Allocate(ny);}
  const char *Name() {return "Exponential Euler";}
  Solve_RC Solve();
  void TimestepDependence() {
    Exponential::TimestepDependence(dt,ny);
  }
  void Source(const vector2& Src, const vector2& Y, double t) {
    OdeProblem->NonLinearSource(Src,Y,t);
  }
};

class I_Euler : public Euler, public IntegratingFactor {
public:
  void Allocate() {Euler::Allocate(); IntegratingFactor::Allocate(ny);}
  const char *Name() {return "Euler w/Integrating Factor";}
  Solve_RC Solve();
  void TimestepDependence() {
    IntegratingFactor::TimestepDependence(dt,ny);
  }
  void Source(const vector2& Src, const vector2& Y, double t) {
    OdeProblem->NonLinearSource(Src,Y,t);
  }
};

class RB1 : public Euler {
public:
  const char *Name() {return "First-Order Rosenbrock";}
  Solve_RC Solve();
};

class I_PC : public PC, public IntegratingFactor {
  double halfdt;
public:
  void Allocate() {PC::Allocate(); IntegratingFactor::Allocate(ny);}
  const char *Name() {return "Predictor-Corrector w/Integrating Factor";}
  void TimestepDependence() {
    PC::TimestepDependence();
    IntegratingFactor::TimestepDependence(dt,ny);
  }
  void Predictor();
  int Corrector();
  void Source(const vector2& Src, const vector2& Y, double t) {
    OdeProblem->NonLinearSource(Src,Y,t);
  }
};

class E_PC : public PC, public Exponential {
protected:
  double dtinv;
public:
  void Allocate() {PC::Allocate(); Exponential::Allocate(ny);}
  const char *Name() {return "Exponential Predictor-Corrector";}
  void TimestepDependence() {
    Exponential::TimestepDependence(dt,ny);
    dtinv=1.0/dt;
  }
  void Predictor();
  int Corrector();
  void Source(const vector2& Src, const vector2& Y, double t) {
    OdeProblem->NonLinearSource(Src,Y,t);
  }
};

class LE_PC : public E_PC {
public:
  const char *Name() {return "Linearized Exponential Predictor-Corrector";}
  void Predictor();
  int Corrector();
};

class CorrectC_PC {
 public:
  bool Correct(Real y0, Real& y,
	       Real source0, Real source, double dt);
  bool Correct(const Complex& y0, Complex& y,
	       const Complex& source0, const Complex& source, double dt);
};

class C_PC : public PC, public CorrectC_PC {
protected:  
public:
  const char *Name() {return "Conservative Predictor-Corrector";}
  int Corrector();
};

OdeVocabulary::OdeVocabulary()
{
  Vocabulary=this;
	
  VOCAB(nu0,(Var) 0.0,(Var) 0.0,"");
  VOCAB(A,(Var) 0.0,(Var) 0.0,"");
  VOCAB(B,(Var) 0.0,(Var) 0.0,"");
	
  METHOD(Ode);
	
  INTEGRATOR(I_Euler);
  INTEGRATOR(E_Euler);
  INTEGRATOR(RB1);
  INTEGRATOR(I_PC);
  INTEGRATOR(E_PC);
  INTEGRATOR(LE_PC);
  INTEGRATOR(C_PC);
  INTEGRATOR(Implicit);
}

OdeVocabulary Ode_Vocabulary;

ofstream fout;

void Ode::InitialConditions()
{
  NY[0]=1;
  Allocate();
  
  Allocate1(nu,ny);
	
  vector y=Y[0];
  y[0]=1.0;
  nu[0]=nu0;
	
  open_output(fout,dirsep,downcase(undashify(Integrator->Abbrev())));
}

void Ode::Initialize()
{
  fout << "# " << Integrator->Name() << endl;
}

void Ode::Output(int)
{
  fout << t << "\t" << abs(y[0]) << "\t" << endl;
}

void Ode::NonLinearSource(const vector2& Src, const vector2& Y, double)
{
  vector source=Src[0];
  vector y=Y[0];
//  source[0]=cos(y[0]);
//  source[0]=cos(t)*y[0];
  source[0]=-A*y[0]-B*y[0]*y[0];
}

void Ode::LinearSource(const vector2& Src, const vector2& Y, double)
{
  vector source=Src[0];
  vector y=Y[0];
  source[0] -= nu[0]*y[0];
}

Solve_RC E_Euler::Solve()
{
  Source(Src,Y,t);
  Problem->Transform(Y,t,dt,YI);
  for(unsigned int j=0; j < ny; j++)
    y[j]=expinv[j]*y[j]+onemexpinv[j]*source[j];
  Problem->BackTransform(Y,t+dt,dt,YI);
  return SUCCESSFUL;
}

Solve_RC I_Euler::Solve()
{
  Source(Src,Y,t);
  Problem->Transform(Y,t,dt,YI);
  for(unsigned int j=0; j < ny; j++)
    y[j]=(y[j]+dt*source[j])*expinv[j];
  Problem->BackTransform(Y,t+dt,dt,YI);
  return SUCCESSFUL;
}

Solve_RC RB1::Solve()
{
#if 0	 // ** JCB 
  Source(Src,Y,t);
  Problem->Transform(Y,t,dt,YI);
  for(int j=0; j < ny; j++)
    y[j] += dt*source[0]/(1.0-dt*(-sin(y[j])-nu[0]));
  Problem->BackTransform(Y0,t+dt,dt,YI);
#endif	
  return SUCCESSFUL;
}

void I_PC::Predictor()
{
  for(unsigned int j=0; j < ny; j++)
    y[j]=(y0[j]+dt*source0[j])*expinv[j];
}

int I_PC::Corrector()
{
  Source(Src,Y,t+dt);
  if(dynamic) {
    for(unsigned int j=0; j < ny; j++) {
      Var pred=y[j];
      y[j]=y0[j]*expinv[j]+halfdt*(source0[j]*expinv[j]+source[j]);
      if(!errmask || errmask[j]) 
	CalcError(y0[j]*expinv[j],y[j],pred,y[j]);
    }
    ExtrapolateTimestep();
  } else for(unsigned int j=0; j < ny; j++) {
    y[j]=y0[j]*expinv[j]+halfdt*(source0[j]*expinv[j]+source[j]);
  }
  return 1;
}

void E_PC::Predictor()
{
  for(unsigned int j=0; j < ny; j++)	
    y[j]=expinv[j]*y0[j]+onemexpinv[j]*source0[j];
}

int E_PC::Corrector()
{
  Source(Src,Y,t+dt);
  if(dynamic) {
    for(unsigned int j=0; j < ny; j++) {
      Var temp=0.5*(source0[j]+source[j]);
      y[j]=expinv[j]*y0[j]+onemexpinv[j]*temp;
      if(!errmask || errmask[j]) 
	CalcError(y0[j]*dtinv,temp,source0[j],temp);
    }
    ExtrapolateTimestep();
  } else for(unsigned int j=0; j < ny; j++) {
    source[j]=0.5*(source0[j]+source[j]);
    y[j]=expinv[j]*y0[j]+onemexpinv[j]*source[j];
  }
  return 1;
}

void LE_PC::Predictor()
{
  //	nu[0]=A+2.0*B*y0[0];
  nu[0]=(y0[0] != 0.0) ? -source0[0]/y0[0] : 0.0;
  source0[0] += nu[0]*y0[0];
  TimestepDependence();
  for(unsigned int j=0; j < ny; j++)
    y[j]=expinv[j]*y0[j]+onemexpinv[j]*source0[j];
}

int LE_PC::Corrector()
{
  Source(Src,Y,t+dt);
  source[0] += nu[0]*y0[0];
  if(dynamic) {
    for(unsigned int j=0; j < ny; j++) {
      source[j]=0.5*(source0[j]+source[j]);
      y[j]=expinv[j]*y0[j]+onemexpinv[j]*source[j];
    }
    for(unsigned int j=0; j < ny; j++)
      if(!errmask || errmask[j]) 
	CalcError(y0[j]*dtinv,source[j],source0[j],source[j]);
    ExtrapolateTimestep();
  } else for(unsigned int j=0; j < ny; j++)
    y[j]=expinv[j]*y0[j]+onemexpinv[j]*0.5*(source0[j]+source[j]);
		
  return 1;
}

inline bool CorrectC_PC::Correct(Real y0, Real& y,
				 Real source0, Real source,
				 double dt)
{
  Real discr=y0*y0+dt*(y0*source0+y*source);
  if(discr >= 0.0) y=sgn(y)*sqrt(discr);
  else {
    if(hybrid) y=y0+0.5*dt*(source0+source);
    else return false;
  }
  return true;
}

inline bool CorrectC_PC::Correct(const Complex& y0, Complex& y,
				 const Complex& source0, const Complex& source,
				 double dt)
{
  if(!Correct(y0.re,y.re,source0.re,source.re,dt)) return false;
  return Correct(y0.im,y.im,source0.im,source.im,dt);
}

int C_PC::Corrector()
{
  Source(Src,Y,t+dt);
  if(dynamic) {
    for(unsigned int j=0; j < ny; j++) {
      Var pred=y[j];
      if(!Correct(y0[j],y[j],source0[j],source[j],dt)) return 0;
      if(!errmask || errmask[j]) CalcError(y0[j],y[j],pred,y[j]);
    }
    ExtrapolateTimestep();
  } else for(unsigned int j=0; j < ny; j++)
    if(!Correct(y0[j],y[j],source0[j],source[j],dt)) return 0;
  
  return 1;
}

Solve_RC Implicit::Solve()
{
  Source(Src,Y,t);
  cout << endl;
  
  for(unsigned int j=0; j < ny; j++) {
#if 0		
    // Implicit midpoint rule		
    Var a=0.25*dt*B;
    Var b=dt*(B*y[j]+0.5*A)+1.0;
    Var c=(dt*(0.25*B*y[j]+0.5*A)-1.0)*y[j];
    if(a != 0.0) {
      y[j]=(-b+sqrt(b*b-4.0*a*c))/(2.0*a);
    }
    else y[j]=-c/b;
#endif
		
#if 0
    // Backward Euler
    Var a=dt*B;
    Var b=1.0+dt*A;
    Var c=-y[j];
    if(a != 0.0) y[j]=(-b+sqrt(b*b-4.0*a*c))/(2.0*a);
    else y[j]=-c/b;
#endif
		
		
#if 1
    // Linearized Backward Euler
    y[j]=y[j]+dt*source[j]/(1.0-dt*(-A-2.0*B*y[j]));
    //		Complex s=-A-2.0*B*y[j];
    //		y[j]=exp(s*dt)*(y[j]+dt*(source[j]-s*y[j])/(1.0-dt*(-A-2.0*B*y[j]-s)));
#endif
		
#if 0
    // Linearized Trapezoidal
    y[j]=y[j]+dt*source[j]/(1.0-0.5*dt*(-A-2.0*B*y[j]));
#endif		
    cout << y[j] << endl;
  }
	
  return SUCCESSFUL;
}

