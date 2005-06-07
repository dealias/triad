#ifndef __Exponential_h__
#define __Exponential_h__ 1

#include "options.h"
#include "kernel.h"
#include "expratio.h"

// Exponential integrators

typedef Array::Array1<Nu>::opt NuVector;

class Exponential {
protected:
  unsigned int start,stop;
  unsigned int startN,stopN;
  NuVector coeff0,coeff1;
public:
  virtual ~Exponential() {}
  
  virtual void Allocator(unsigned int) {};
  
  void Allocator() {
    unsigned int ny=stop-start;
    Allocate(coeff0,ny,start,0);
    Allocate(coeff1,ny,start,0);
    Allocator(ny);
  }
  
};

template<class T>
class E_PC : public PC, public Exponential {
protected:
  T *parent;
public:
  E_PC(T *parent) : parent(parent) {}
  
  const char *Name() {return "Exponential Predictor-Corrector";}
  
  void Allocator() {
    PC::Allocator();
    parent->IndexLimits(start,stop,startN,stopN);
    Exponential::Allocator();
  }
  
  void Source(const vector2& Src, const vector2& Y, double t) {
    parent->ExponentialSource(Src,Y,t);
  }
  
  void PSource(const vector2&, const vector2&, double) {}
  void CSource(const vector2&, const vector2&, double) {}
  
  inline void TimestepDependence(),Predictor(),Corrector();
  
  virtual void Predictor(unsigned int, unsigned int) {
    E_PC<T>::Predictor();
    PC::Predictor(startN,stopN);
  }
  
  virtual int Corrector(unsigned int, unsigned int) {
    E_PC<T>::Corrector();
    return PC::Corrector(startN,stopN);
  }
};

template<class T>
inline void E_PC<T>::TimestepDependence() {
  PC::TimestepDependence();
  for(unsigned int j=start; j < stop; j++) {
    Nu nuk=parent->LinearCoeff(j);
    Nu x=-nuk*dt;
    Nu temp=expratio1(x); // (e^x-1)/x
    coeff0[j]=x*temp+1.0;
    coeff1[j]=temp*dt;
  }
}

template<class T>
inline void E_PC<T>::Predictor() {
  for(unsigned int j=start; j < stop; j++)	
    y[j]=coeff0[j]*y0[j]+coeff1[j]*source0[j];
}
  
template<class T>
inline void E_PC<T>::Corrector() {
  Source(Src,Y,t+dt);
  if(dynamic) {
    for(unsigned int j=start; j < stop; j++) {
      Var temp=coeff0[j]*y0[j];
      Var val=temp+coeff1[j]*0.5*(source0[j]+source[j]);
      if(!Active(errmask) || errmask[j])
	CalcError(y0[j],val,temp+coeff1[j]*source0[j],val);
      y[j]=val;
    }
  } else {
    for(unsigned int j=start; j < stop; j++)
      y[j]=coeff0[j]*y0[j]+coeff1[j]*0.5*(source0[j]+source[j]);
  }
}
  
template<class T>
class E_RK2 : public RK2, public Exponential {
protected:
  T *parent;
  NuVector coeff0h,coeff1h;
public:
  E_RK2(T *parent) : parent(parent) {}
  
  const char *Name() {return "Second-Order Exponential Runge-Kutta";}
  
  void Allocator() {
    RK2::Allocator();
    parent->IndexLimits(start,stop,startN,stopN);
    Exponential::Allocator();
  }
  
  void Allocator(unsigned int ny) {
    Allocate(coeff0h,ny,start,0);
    Allocate(coeff1h,ny,start,0);
  }
  
  void Source(const vector2& Src, const vector2& Y, double t) {
    parent->ExponentialSource(Src,Y,t);
  }
  
  void PSource(const vector2&, const vector2&, double) {}
  void CSource(const vector2&, const vector2&, double) {}
  
  inline void TimestepDependence(),Predictor(),Corrector();
  
  virtual void Predictor(unsigned int, unsigned int) {
    E_RK2<T>::Predictor();
    RK2::Predictor(startN,stopN);
  }
  
  virtual int Corrector(unsigned int, unsigned int) {
    E_RK2<T>::Corrector();
    return RK2::Corrector(startN,stopN);
  }
};

template<class T>
inline void E_RK2<T>::TimestepDependence() {
  PC::TimestepDependence();
  for(unsigned int j=start; j < stop; j++) {
    Nu nuk=parent->LinearCoeff(j);
    Real x=-nuk*dt;
    Nu temph=0.5*expratio1(0.5*x);
    coeff0h[j]=x*temph+1.0;
    coeff1h[j]=temph*dt;
    Nu temp=expratio1(x);
    coeff0[j]=x*temp+1.0;
    coeff1[j]=temp*dt;
  }
}

template<class T>
inline void E_RK2<T>::Predictor()
{
  for(unsigned int j=start; j < stop; j++)
    y[j]=coeff0h[j]*y0[j]+coeff1h[j]*source0[j];
  Problem->BackTransform(Y,t+halfdt,halfdt,YI);
}

template<class T>
inline void E_RK2<T>::Corrector()
{
  Source(Src,Y,t+halfdt);
  if(dynamic) {
    for(unsigned int j=start; j < stop; j++) {
      Var val=coeff0[j]*y0[j]+coeff1[j]*source[j];
      if(!Active(errmask) || errmask[j])
	CalcError(y0[j],val,coeff0[j]*y0[j]+coeff1[j]*source0[j],val);
      val=y[j];
    }
  } else {
    for(unsigned int j=start; j < stop; j++)
      y[j]=coeff0[j]*y0[j]+coeff1[j]*source[j];
  }
}

class RK3_Exponential : public Exponential {
public:  
  NuVector coeff0h,coeff1h,coeffA,coeffB,coeffC;
  
  void Allocator(unsigned int ny, unsigned int start) {
    Allocate(coeff0h,ny,start,0);
    Allocate(coeff1h,ny,start,0);
    Allocate(coeffA,ny,start,0);
    Allocate(coeffB,ny,start,0);
    Allocate(coeffC,ny,start,0);
  }
  
  virtual Nu LinearCoeff(unsigned int j)=0;
  
  void TimestepDependence(Real dt) {
    for(unsigned int j=start; j < stop; j++) {
      Nu nuk=LinearCoeff(j);
      Real x=-nuk*dt;
      Nu temph=0.5*expratio1(0.5*x);
      coeff0h[j]=x*temph+1.0;
      coeff1h[j]=temph*dt;
      Nu temp3=expratio3(x);
      Nu temp=x*x*temp3+0.5*x+1.0;
      coeff0[j]=x*temp+1.0;
      coeff1[j]=temp*dt;
      coeffA[j]=((x*x-3.0*x+4.0)*temp3+0.5*(x-1.0))*dt;
      coeffB[j]=2.0*((x-2.0)*temp3+0.5)*dt;
      coeffC[j]=(-(x-4.0)*temp3-0.5)*dt;
    }
  }
};

template<class T>
class E_RK3 : public RK3, public RK3_Exponential {
protected:
  T *parent;
public:
  E_RK3(T *parent) : parent(parent) {order=3;}
  
  const char *Name() {return "Third-Order Exponential Runge-Kutta";}
  
  void Allocator() {
    RK3::Allocator();
    parent->IndexLimits(start,stop,startN,stopN);
    Exponential::Allocator();
  }
  
  void Allocator(unsigned int ny) {
    RK3_Exponential::Allocator(ny,start);
  }
  
  Nu LinearCoeff(unsigned int j) {
    return parent->LinearCoeff(j);
  };
  
  void Source(const vector2& Src, const vector2& Y, double t) {
    parent->ExponentialSource(Src,Y,t);
  }
  
  void PSource(const vector2&, const vector2&, double) {}
  void CSource(const vector2&, const vector2&, double) {}
  
  void TimestepDependence() {
    RK3::TimestepDependence();
    RK3_Exponential::TimestepDependence(dt);
  }

  inline void Predictor(),Corrector();
  
  virtual void Predictor(unsigned int, unsigned int) {
    E_RK3<T>::Predictor();
    RK3::Predictor(startN,stopN);
  }
  
  virtual int Corrector(unsigned int, unsigned int) {
    E_RK3<T>::Corrector();
    return RK3::Corrector(startN,stopN);
  }
};

template<class T>
void E_RK3<T>::Predictor()
{
  for(unsigned int j=start; j < stop; j++)	
    y[j]=coeff0h[j]*y0[j]+coeff1h[j]*source0[j];
  Source(Src1,Y,t+halfdt);
  for(unsigned int j=start; j < stop; j++)	
    y[j]=coeff0[j]*y0[j]+coeff1[j]*(2.0*source1[j]-source0[j]);
}

template<class T>
void E_RK3<T>::Corrector()
{
  Source(Src,Y,t+dt);
  if(dynamic) {
    for(unsigned int j=start; j < stop; j++) {
      Var pred=coeff0[j]*y0[j]+coeff1[j]*source1[j];
      Var val=coeff0[j]*y0[j]+coeffA[j]*source0[j]+2.0*coeffB[j]*source1[j]
	+coeffC[j]*source[j];
      if(!Active(errmask) || errmask[j])
	CalcError(y0[j],val,pred,val);
      y[j]=val;
    }
  } else {
    for(unsigned int j=start; j < stop; j++)
      y[j]=coeff0[j]*y0[j]+coeffA[j]*source0[j]+2.0*coeffB[j]*source1[j]
	+coeffC[j]*source[j];
  }
}

template<class T>
class E_RK4 : public RK4, public RK3_Exponential {
protected:
  T *parent;
  vector source3;
  vector2 Src3;
public:
  E_RK4(T *parent) : parent(parent) {order=4;}
  
  const char *Name() {return "Fourth-Order Exponential Runge-Kutta";}
  
  void Allocator() {
    RK4::Allocator();
    parent->IndexLimits(start,stop,startN,stopN);
    Exponential::Allocator();
  }
  
  void Allocator(unsigned int ny) {
    RK3_Exponential::Allocator(ny,start);
    if(dynamic) Alloc0(Src3,source3);
  }
  
  Nu LinearCoeff(unsigned int j) {
    return parent->LinearCoeff(j);
  };
  
  void Source(const vector2& Src, const vector2& Y, double t) {
    parent->ExponentialSource(Src,Y,t);
  }
  
  void PSource(const vector2&, const vector2&, double) {}
  void CSource(const vector2&, const vector2&, double) {}
  
  void TimestepDependence() {
    RK4::TimestepDependence();
    RK3_Exponential::TimestepDependence(dt);
  }

  inline void Predictor(),Corrector();
  
  virtual void Predictor(unsigned int, unsigned int) {
    E_RK4<T>::Predictor();
    RK4::Predictor(startN,stopN);
  }
  
  virtual int Corrector(unsigned int, unsigned int) {
    E_RK4<T>::Corrector();
    return RK4::Corrector(startN,stopN);
  }
};

template<class T>
inline void E_RK4<T>::Predictor()
{
  for(unsigned int j=start; j < stop; j++)	
    y[j]=coeff0h[j]*y0[j]+coeff1h[j]*source0[j];
  Source(Src1,Y,t+halfdt);
  for(unsigned int j=start; j < stop; j++)	
    source[j]=coeff0h[j]*y0[j]+coeff1h[j]*source1[j];
  Source(Src2,Src,t+halfdt);
  for(unsigned int j=start; j < stop; j++)	
    y[j]=coeff0h[j]*y[j]+coeff1h[j]*(2.0*source2[j]-source0[j]);
}

template<class T>
inline void E_RK4<T>::Corrector()
{
  Source(Src,Y,t+dt);
  
  if(dynamic) {
    for(unsigned int j=start; j < stop; j++)	
      y[j]=coeff0[j]*y0[j]+coeff1[j]*(2.0*source1[j]-source0[j]);
    Source(Src3,Y,t+dt);
    for(unsigned int j=start; j < stop; j++) {
      Var val=coeff0[j]*y0[j]+coeffA[j]*source0[j]
	+coeffB[j]*(source1[j]+source2[j])+coeffC[j]*source[j];
      if(!Active(errmask) || errmask[j]) {
	Var pred=coeff0[j]*y0[j]+coeffA[j]*source0[j]+2.0*coeffB[j]*source1[j]
	  +coeffC[j]*source3[j];
	CalcError(y0[j],val,pred,val);
	y[j]=val;
      }
    }
  } else {
    for(unsigned int j=start; j < stop; j++)
      y[j]=coeff0[j]*y0[j]+coeffA[j]*source0[j]
	+coeffB[j]*(source1[j]+source2[j])+coeffC[j]*source[j];
  }
}

#endif
