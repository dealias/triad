#include "options.h"
#include "kernel.h"

const char *method="Ode";
const char *integrator="PC";

// Global variables
Var A;
Var B;
Nu *nu;

// Local vocabulary declarations and default values
static Nu nu0=1.0;

class OdeVocabulary : public VocabularyBase {
public:
  char *Name() {return "Ode";}
  char *Abbrev() {return "Ode";}
  OdeVocabulary();
};

class Ode : public ProblemBase {
public:
  Ode();
	
  void InitialConditions();
  void Initialize();
  void Output(int it);
	
  void NonLinearSrc(Var *, Var *, double);
  void LinearSrc(Var *, Var *, double);
	
  void Source(Var *src, Var *Y, double t) {
    NonLinearSrc(src,Y,t);
    LinearSrc(src,Y,t);
  }
};

Ode *OdeProblem;

Ode::Ode()
{
  OdeProblem=this;
}

class Exponential {
protected:
  Nu *nu_inv;
  Nu *expinv,*onemexpinv;
public:
  void Allocate(int n) {
    expinv=new Nu[n];
    onemexpinv=new Nu[n];
	
    nu_inv=new Nu[n];

    for(int j=0; j < n; j++) {
      if(nu[j] != 0.0) nu_inv[j]=1.0/nu[j];
    }
  }

  void TimestepDependence(double dt, int n) {
    for(int j=0; j < n; j++) {
      Nu temp=-expm1(-nu[j]*dt);
      expinv[j]=1.0-temp;
      if(nu[j] == 0.0) onemexpinv[j]=dt;
      else onemexpinv[j]=temp*nu_inv[j];
    }
  }
};

class IntegratingFactor {
protected:
  Nu *expinv;
public:
  void Allocate(int n) {expinv=new Nu[n];}
  void TimestepDependence(double dt, int n) {
    for(int j=0; j < n; j++) expinv[j]=exp(-nu[j]*dt);
  }
};
	
class Implicit : public IntegratorBase {
public:
  void Allocate(int n) {IntegratorBase::Allocate(n);}
  const char *Name() {return "Implicit";}
  Solve_RC Solve(double, double&);
};

class E_Euler : public Euler, public Exponential {
public:
  void Allocate(int n) {Euler::Allocate(n); Exponential::Allocate(n);}
  char *Name() {return "Exponential Euler";}
  Solve_RC Solve(double, double);
  void TimestepDependence(double dt) {
    Exponential::TimestepDependence(dt,ny);
  }
  void Source(Var *src, Var *Y, double t) {
    OdeProblem->NonLinearSrc(src,Y,t);
  }
};

class I_Euler : public Euler, public IntegratingFactor {
public:
  void Allocate(int n) {Euler::Allocate(n); IntegratingFactor::Allocate(n);}
  char *Name() {return "Euler w/Integrating Factor";}
  Solve_RC Solve(double, double);
  void TimestepDependence(double dt) {
    IntegratingFactor::TimestepDependence(dt,ny);
  }
  void Source(Var *src, Var *Y, double t) {
    OdeProblem->NonLinearSrc(src,Y,t);
  }
};

class RB1 : public Euler {
public:
  char *Name() {return "First-Order Rosenbrock";}
  Solve_RC Solve(double, double);
};

class I_PC : public PC, public IntegratingFactor {
public:
  void Allocate(int n) {PC::Allocate(n); IntegratingFactor::Allocate(n);}
  char *Name() {return "Predictor-Corrector w/Integrating Factor";}
  void TimestepDependence(double dt) {
    IntegratingFactor::TimestepDependence(dt,ny);
  }
  void Predictor(double, double, unsigned int, unsigned int);
  int Corrector(double, int, unsigned int, unsigned int);
  void Source(Var *src, Var *Y, double t) {
    OdeProblem->NonLinearSrc(src,Y,t);
  }
};

class E_PC : public PC, public Exponential {
protected:
  double dtinv;
public:
  void Allocate(int n) {PC::Allocate(n); Exponential::Allocate(n);}
  char *Name() {return "Exponential Predictor-Corrector";}
  void TimestepDependence(double dt) {
    Exponential::TimestepDependence(dt,ny);
    dtinv=1.0/dt;
  }
  void Predictor(double, double, unsigned int, unsigned int);
  int Corrector(double, int, unsigned int, unsigned int);
  void Source(Var *src, Var *Y, double t) {
    OdeProblem->NonLinearSrc(src,Y,t);
  }
};

class LE_PC : public E_PC {
public:
  char *Name() {return "Linearized Exponential Predictor-Corrector";}
  void Predictor(double, double, unsigned int, unsigned int);
  int Corrector(double, int, unsigned int, unsigned int);
};

OdeVocabulary::OdeVocabulary()
{
  Vocabulary=this;
	
  VOCAB(nu0,(Var) 0.0,(Var) 0.0);
  VOCAB(A,(Var) 0.0,(Var) 0.0);
  VOCAB(B,(Var) 0.0,(Var) 0.0);
	
  METHOD(Ode);
	
  INTEGRATOR(I_Euler);
  INTEGRATOR(E_Euler);
  INTEGRATOR(RB1);
  INTEGRATOR(I_PC);
  INTEGRATOR(E_PC);
  INTEGRATOR(LE_PC);
  INTEGRATOR(Implicit);
}

OdeVocabulary Ode_Vocabulary;

ofstream fout;

void Ode::InitialConditions()
{
  ny=1;
  y=new Var[ny];
  nu=new Var[ny];
	
  y[0]=1.0;
  nu[0]=nu0;
	
  dynamic=0;
  open_output(fout,dirsep,downcase(undashify(Integrator->Abbrev())));
}

void Ode::Initialize()
{
  fout << "# " << Integrator->Name() << endl;
}

void Ode::Output(int)
{
  fout << t << "\t" << y[0].re << "\t" << endl;
}

void Ode::NonLinearSrc(Var *source, Var *, double t)
{
  //	source[0]=cos(y[0]);
  //	source[0]=cos(t)*y[0];
  source[0]=-A*y[0]-B*y[0]*y[0];
}

void Ode::LinearSrc(Var *source, Var *y, double)
{
  source[0] -= nu[0]*y[0];
}

Solve_RC E_Euler::Solve(double t, double dt)
{
  Source(source,y0,t);
  Problem->Transform(y0,t,dt,yi);
  for(unsigned int j=0; j < ny; j++)
    y0[j]=expinv[j]*y0[j]+onemexpinv[j]*source[j];
  Problem->BackTransform(y0,t+dt,dt,yi);
  return SUCCESSFUL;
}

Solve_RC I_Euler::Solve(double t, double dt)
{
  Source(source,y0,t);
  Problem->Transform(y0,t,dt,yi);
  for(unsigned int j=0; j < ny; j++)
    y0[j]=(y0[j]+dt*source[j])*expinv[j];
  Problem->BackTransform(y0,t+dt,dt,yi);
  return SUCCESSFUL;
}

Solve_RC RB1::Solve(double t, double dt)
{
#if 0	 // ** JCB 
  Source(source,y0,t);
  Problem->Transform(y0,t,dt,yi);
  for(int j=0; j < ny; j++)
    y0[j]=y0[j]+dt*source[0]/(1.0-dt*(-sin(y0[j])-nu[0]));
  Problem->BackTransform(y0,t+dt,dt,yi);
#endif	
  return SUCCESSFUL;
}

void I_PC::Predictor(double t, double dt, unsigned int start,
		     unsigned int stop)
{
  for(unsigned int j=start; j < stop; j++)
    y1[j]=(y0[j]+dt*source0[j])*expinv[j];
  Source(source,y1,t+dt);
}

int I_PC::Corrector(double dt, int dynamic, unsigned int start,
		    unsigned int stop)
{
  const double halfdt=0.5*dt;
  if(dynamic) for(unsigned int j=start; j < stop; j++) {
    Var pred=y[j];
    y[j]=y0[j]*expinv[j]+halfdt*(source0[j]*expinv[j]+source[j]);
    if(!errmask || errmask[j]) 
      CalcError(y0[j]*expinv[j],y[j],pred,y[j]);
  } else for(unsigned int j=start; j < stop; j++) {
    y[j]=y0[j]*expinv[j]+halfdt*(source0[j]*expinv[j]+source[j]);
  }
  return 1;
}

void E_PC::Predictor(double t, double, unsigned int start, unsigned int stop)
{
  for(unsigned int j=start; j < stop; j++)
    y1[j]=expinv[j]*y0[j]+onemexpinv[j]*source0[j];
  Source(source,y1,t+dt);
}

int E_PC::Corrector(double, int dynamic, unsigned int start, unsigned int stop)
{
  unsigned int j;
  if(dynamic) {
    for(j=start; j < stop; j++) {
      source[j]=0.5*(source0[j]+source[j]);
      y[j]=expinv[j]*y0[j]+onemexpinv[j]*source[j];
    }
    for(j=start; j < stop; j++)
      if(!errmask || errmask[j]) 
	CalcError(y0[j]*dtinv,source[j],source0[j],source[j]);
    ExtrapolateTimestep();
  } else for(j=start; j < stop; j++)
    y[j]=expinv[j]*y0[j]+onemexpinv[j]*0.5*(source0[j]+source[j]);
		
  return 1;
}

void LE_PC::Predictor(double t, double dt, unsigned int start,
		      unsigned int stop)
{
  //	nu[0]=A+2.0*B*y0[0];
  nu[0]=(y0[0] != 0.0) ? -source0[0]/y0[0] : 0.0;
  source0[0] += nu[0]*y0[0];
  TimestepDependence(dt);
  for(unsigned int j=start; j < stop; j++)
    y1[j]=expinv[j]*y0[j]+onemexpinv[j]*source0[j];
  Source(source,y1,t+dt);
}

int LE_PC::Corrector(double, int dynamic, unsigned int start,
		     unsigned int stop)
{
  unsigned int j;
  source[0] += nu[0]*y0[0];
  if(dynamic) {
    for(j=start; j < stop; j++) {
      source[j]=0.5*(source0[j]+source[j]);
      y[j]=expinv[j]*y0[j]+onemexpinv[j]*source[j];
    }
    for(j=start; j < stop; j++)
      if(!errmask || errmask[j]) 
	CalcError(y0[j]*dtinv,source[j],source0[j],source[j]);
    ExtrapolateTimestep();
  } else for(j=start; j < stop; j++)
    y[j]=expinv[j]*y0[j]+onemexpinv[j]*0.5*(source0[j]+source[j]);
		
  return 1;
}

Solve_RC Implicit::Solve(double t, double& dt)
{
  Source(source,y0,t);
  cout << endl;
  for(unsigned int j=0; j < ny; j++) {
#if 0		
    // Implicit midpoint rule		
    Var a=0.25*dt*B;
    Var b=dt*(B*y0[j]+0.5*A)+1.0;
    Var c=(dt*(0.25*B*y0[j]+0.5*A)-1.0)*y0[j];
    if(a != 0.0) {
      y0[j]=(-b+sqrt(b*b-4.0*a*c))/(2.0*a);
    }
    else y0[j]=-c/b;
#endif
		
#if 0
    // Backward Euler
    Var a=dt*B;
    Var b=1.0+dt*A;
    Var c=-y0[j];
    if(a != 0.0) y0[j]=(-b+sqrt(b*b-4.0*a*c))/(2.0*a);
    else y0[j]=-c/b;
#endif
		
		
#if 1
    // Linearized Backward Euler
    y0[j]=y0[j]+dt*source[j]/(1.0-dt*(-A-2.0*B*y0[j]));
    //		Complex s=-A-2.0*B*y0[j];
    //		y0[j]=exp(s*dt)*(y0[j]+dt*(source[j]-s*y0[j])/(1.0-dt*(-A-2.0*B*y0[j]-s)));
#endif
		
#if 0
    // Linearized Trapezoidal
    y0[j]=y0[j]+dt*source[j]/(1.0-0.5*dt*(-A-2.0*B*y0[j]));
#endif		
    cout << y0[j] << endl;
  }
	
  return SUCCESSFUL;
}

