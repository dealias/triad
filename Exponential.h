#ifndef __Exponential_h__
#define __Exponential_h__ 1

#include "options.h"
#include "kernel.h"
#include "expratio.h"

// Exponential integrators

typedef Array::array1<Nu>::opt Nuvector;

class Exponential {
protected:
  Nuvector coeff0,coeff1;
  unsigned Nspecial;
public:
  virtual ~Exponential() {}
  
  unsigned NSpecial();
  void Source(const vector2& Src, const vector2& Y, double t);
  Nu LinearCoeff(unsigned int i);
  
  void Allocator() {
    Nspecial=NSpecial();
    Allocate(coeff0,Nspecial);
    Allocate(coeff1,Nspecial);
  }
};

class E_PC : public PC, public Exponential {
protected:
  Real dtinv;
public:
  void Allocator() {PC::Allocator(); Exponential::Allocator();}
  const char *Name() {return "Exponential Predictor-Corrector";}
  void TimestepDependence() {
    PC::TimestepDependence();
    dtinv=1.0/dt;
    for(unsigned int j=0; j < Nspecial; j++) {
      Nu nuk=LinearCoeff(j);
      Nu x=-nuk*dt;
      Nu temp=expratio1(x); // (e^x-1)/x
      coeff0[j]=x*temp+1.0;
      coeff1[j]=temp*dt;
    }
  }
  void Predictor();
  int Corrector();
};

class E_RK2 : public PC, public Exponential {
protected:
  Nuvector coeff2;
public:
  void Allocator() {
    PC::Allocator();
    Exponential::Allocator();
    Allocate(coeff2,Nspecial);
  }
  const char *Name() {return "Second-Order Exponential Runge-Kutta";}
  void TimestepDependence() {
    PC::TimestepDependence();
    for(unsigned int j=0; j < Nspecial; j++) {
      Nu nuk=LinearCoeff(j);
      Real x=-nuk*dt;
      Nu temp2=expratio2(x);
      Nu temp=x*temp2+1.0;
      coeff0[j]=x*temp+1.0;
      coeff1[j]=temp*dt;
      coeff2[j]=temp2*dt;
      }
    }

  void Predictor();
  int Corrector();
};

class E_RK3 : public PC, public Exponential {
protected:
  Nuvector coeff0h,coeff1h,coeffA,coeffB,coeffC;
  Nuvector coeff2;
  vector source1,source2;
  vector2 Src1,Src2;
public:
  E_RK3() {order=3;}
  void Allocator0() {
    PC::Allocator();
    Exponential::Allocator();
    Allocate(coeff0h,Nspecial);
    Allocate(coeff1h,Nspecial);
    Allocate(coeffA,Nspecial);
    Allocate(coeffB,Nspecial);
    Allocate(coeffC,Nspecial);
    Alloc0(Src1,source1);
    Alloc0(Src2,source2);
    
  }
  void Allocator() {
    Allocator0();
    if(dynamic) Allocate(coeff2,Nspecial);
  }
  const char *Name() {return "Third-Order Exponential Runge-Kutta";}
  void TimestepDependence() {
    PC::TimestepDependence();
    for(unsigned int j=0; j < Nspecial; j++) {
      Nu nuk=LinearCoeff(j);
      Real x=-nuk*dt;
      Nu temph=expratio1(0.5*x);
      Nu temp3=expratio3(x);
      Nu temp=x*x*temp3+0.5*x+1.0;
      coeff0[j]=x*temp+1.0;
      coeff0h[j]=0.5*x*temph+1.0;
      coeff1[j]=temp*dt;
      coeff1h[j]=temph*halfdt;
      coeffA[j]=((x*x-3.0*x+4.0)*temp3+0.5*(x-1))*dt;
      coeffB[j]=2.0*((x-2)*temp3+0.5)*dt;
      coeffC[j]=(-(x-4)*temp3-0.5)*dt;
      if(Active(coeff2)) coeff2[j]=(x*temp3+0.5)*dt;
    }
  }
  void Predictor();
  int Corrector();
};

class E_RK4 : public E_RK3 {
protected:
  vector source3;
  vector2 Src3;
public:
  E_RK4() {order=4;}
  void Allocator() {
    E_RK3::Allocator0();
    if(dynamic) Alloc0(Src3,source3);
  }
  const char *Name() {return "Fourth-Order Exponential Runge-Kutta";}
  void Predictor();
  int Corrector();
};

#endif
