#ifndef __Ode_h__
#define __Ode_h__ 1

// Specialized ODE integrators

class Exponential0 {
protected:
  Nuvector nu_inv, expinv, onemexpinv;
public:
  void Allocator(unsigned int n) {
    Allocate(expinv,n);
    Allocate(onemexpinv,n);
    Allocate(nu_inv,n);
    
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
  void Allocator(unsigned int n) {Allocate(expinv,n);}
  void TimestepDependence(double dt, unsigned int n) {
    for(unsigned int j=0; j < n; j++) expinv[j]=exp(-nu[j]*dt);
  }
};
	
class Implicit : public IntegratorBase {
public:
  void Allocator() {IntegratorBase::Allocator();}
  const char *Name() {return "Implicit";}
  Solve_RC Solve(const Duvector& Fields=AllFields);
};

class E_Euler : public Euler, public Exponential0 {
public:
  void Allocator() {Euler::Allocator(); Exponential0::Allocator(ny);}
  const char *Name() {return "Exponential Euler";}
  Solve_RC Solve(const Duvector& Fields=AllFields);
  void TimestepDependence() {
    Exponential0::TimestepDependence(dt,ny);
  }
  void Source(const vector2& Src, const vector2& Y, double t) {
    OdeProblem->NonLinearSource(Src,Y,t);
  }
};

class I_Euler : public Euler, public IntegratingFactor {
public:
  void Allocator() {Euler::Allocator(); IntegratingFactor::Allocator(ny);}
  const char *Name() {return "Euler w/Integrating Factor";}
  Solve_RC Solve(const Duvector& Fields=AllFields);
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
  Solve_RC Solve(const Duvector& Fields=AllFields);
};

class I_PC : public PC, public IntegratingFactor {
  double halfdt;
public:
  void Allocator() {PC::Allocator(); IntegratingFactor::Allocator(ny);}
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

#if 0
class LE_PC : public PC, public Exponential0 {
protected:
  double dtinv;
public:
  void Allocator() {PC::Allocator(); Exponential0::Allocator(ny);}
  const char *Name() {return "Linearized Exponential Predictor-Corrector";}
  void TimestepDependence() {
    Exponential0::TimestepDependence(dt,ny);
    dtinv=1.0/dt;
  }
  void Predictor();
  int Corrector();
  void Source(const vector2& Src, const vector2& Y, double t) {
    OdeProblem->NonLinearSource(Src,Y,t);
  }
};
#endif

Solve_RC E_Euler::Solve(const Duvector& Fields)
{
  unsigned int nFields=Fields.Size();
  Source(Src,Y,t);
  Problem->Transform(Y,t,dt,YI);
  for(unsigned int s=0; s < nFields; s++) {
    unsigned int f=Fields[s]; vector y=Y[f]; vector source=Src[f];
    unsigned int ny=y.Size();
    for(unsigned int j=0; j < ny; j++)
      y[j]=expinv[j]*y[j]+onemexpinv[j]*source[j];
  }
  Problem->BackTransform(Y,t+dt,dt,YI);
  return SUCCESSFUL;
}

Solve_RC I_Euler::Solve(const Duvector& Fields)
{
  unsigned int nFields=Fields.Size();
  Source(Src,Y,t);
  Problem->Transform(Y,t,dt,YI);
  for(unsigned int s=0; s < nFields; s++) {
    unsigned int f=Fields[s]; vector y=Y[f]; vector source=Src[f];
    unsigned int ny=y.Size();
    for(unsigned int j=0; j < ny; j++)
      y[j]=(y[j]+dt*source[j])*expinv[j];
  }
  Problem->BackTransform(Y,t+dt,dt,YI);
  return SUCCESSFUL;
}

Solve_RC RB1::Solve(const Duvector& Fields)
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
  for(unsigned int s=0; s < nFields; s++) {
    unsigned int f=Fields[s]; vector y0=Y0[f]; vector source0=Src0[f];
    vector y=Y[f];
    unsigned int ny=y.Size();
    for(unsigned int j=0; j < ny; j++)
      y[j]=(y0[j]+dt*source0[j])*expinv[j];
  }
}

int I_PC::Corrector()
{
  Source(Src,Y,t+dt);
  for(unsigned int s=0; s < nFields; s++) {
    unsigned int f=Fields[s]; vector y0=Y0[f]; vector source0=Src0[f];
    vector y=Y[f]; vector source=Src[f];
    unsigned int ny=y.Size();
    if(dynamic) {
      for(unsigned int j=0; j < ny; j++) {
	Var pred=y[j];
	y[j]=y0[j]*expinv[j]+halfdt*(source0[j]*expinv[j]+source[j]);
	if(!Active(errmask) || errmask[j])
	  CalcError(y0[j]*expinv[j],y[j],pred,y[j]);
      }
  } else for(unsigned int j=0; j < ny; j++)
    y[j]=y0[j]*expinv[j]+halfdt*(source0[j]*expinv[j]+source[j]);
  }
  return 1;
}

#if 0
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
  for(unsigned int s=0; s < nFields; s++) {
    unsigned int f=Fields[s]; vector y0=Y0[f]; vector source0=Src0[f];
    vector y=Y[f]; vector source=Src[f];
    unsigned int ny=y.Size();
    if(dynamic) {
      for(unsigned int j=0; j < ny; j++) {
	source[j]=0.5*(source0[j]+source[j]);
	y[j]=expinv[j]*y0[j]+onemexpinv[j]*source[j];
      }
      for(unsigned int j=0; j < ny; j++)
	if(!Active(errmask) || errmask[j])
	  CalcError(y0[j]*dtinv,source[j],source0[j],source[j]);
    } else for(unsigned int j=0; j < ny; j++)
      y[j]=expinv[j]*y0[j]+onemexpinv[j]*0.5*(source0[j]+source[j]);
  }
		
  return 1;
}
#endif

Solve_RC Implicit::Solve(const Duvector& Fields)
{
  unsigned int nFields=Fields.Size();
  Source(Src,Y,t);
  cout << endl;
  
  for(unsigned int s=0; s < nFields; s++) {
    unsigned int f=Fields[s]; vector y=Y[f];
    unsigned int ny=y.Size();
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
  }
	
  return SUCCESSFUL;
}

#endif
