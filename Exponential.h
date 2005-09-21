#ifndef __Exponential_h__
#define __Exponential_h__ 1

#include "options.h"
#include "kernel.h"
#include "phi.h"

// Exponential integrators

typedef Array::Array1<Nu>::opt NuVector;

class Exponential {
protected:
  unsigned int start,stop;
  unsigned int startN,stopN;
  unsigned int stop0,start0; // Unused
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
    parent->IndexLimits(start,stop,startN,stop0,start0,stopN);
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
    Nu ph1=phi1(x); // (e^x-1)/x
    coeff0[j]=x*ph1+1.0;
    coeff1[j]=ph1*dt;
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
      if(!Array::Active(errmask) || errmask[j])
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
  vector source1;
  vector2 Src2;
  NuVector coeff0h,coeff1h,b0,b1;
public:
  E_RK2(T *parent) : parent(parent) {}
  
  const char *Name() {return "Second-Order Exponential Runge-Kutta";}
  
  void Allocator() {
    RK2::Allocator();
    parent->IndexLimits(start,stop,startN,stop0,start0,stopN);
    Exponential::Allocator();
  }
  
  void Allocator(unsigned int ny) {
    Allocate(coeff0h,ny,start,0);
    Allocate(coeff1h,ny,start,0);
    Allocate(b0,ny,start,0);
    Allocate(b1,ny,start,0);
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
    Nu x=-nuk*dt;
    Nu xh=0.5*x;
    Nu ph1H=phi1(xh)*halfdt;
    coeff0h[j]=-nuk*ph1H+1;
    coeff1h[j]=ph1H;
    Nu ph2=phi2(x)*dt;
    Nu ph1=x*ph2+dt;
    coeff0[j]=-nuk*ph1+1;
    coeff1[j]=ph1;
    b0[j]=ph1-2*ph2;
    b1[j]=2*ph2;
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
      Var val=coeff0[j]*y0[j]+b0[j]*source0[j]+b1[j]*source[j];
      if(!Array::Active(errmask) || errmask[j])
	CalcError(y0[j],val,coeff0[j]*y0[j]+coeff1[j]*source0[j],val);
      y[j]=val;
    }
  } else {
    for(unsigned int j=start; j < stop; j++)
      y[j]=coeff0[j]*y0[j]+b0[j]*source0[j]+b1[j]*source[j];
  }
}

template<class T>
class E_RK3 : public RK3 {
protected:
  unsigned int start,stop;
  unsigned int startN,stopN;
  unsigned int stop0,start0;
  T *parent;
  NuVector coeff0,coeff0h,coeff0i;
  NuVector a21;
  NuVector a31,a32;
  NuVector b1,b2,b3;
  NuVector B1,B2,B3,B4;
public:
  E_RK3(T *parent) : parent(parent) {order=3;}
  
  const char *Name() {
    return "Third-Order Exponential Bogacki-Shampine Runge-Kutta";
  }
  
  void Allocator() {
    RK3::Allocator();
    parent->IndexLimits(start,stop,startN,stop0,start0,stopN);

    unsigned int ny=stop-start;
    Allocate(coeff0,ny,start,0);
    Allocate(coeff0h,ny,start,0);
    Allocate(coeff0i,ny,start,0);
    Allocate(a21,ny,start,0);
    Allocate(a31,ny,start,0);
    Allocate(a32,ny,start,0);
    Allocate(b1,ny,start,0);
    Allocate(b2,ny,start,0);
    Allocate(b3,ny,start,0);
    if(dynamic) {
      Allocate(B1,ny,start,0);
      Allocate(B2,ny,start,0);
      Allocate(B3,ny,start,0);
      Allocate(B4,ny,start,0);
    }
  }
  
  void TimestepDependence() {
    if(startN < stopN) RK3::TimestepDependence();
    halfdt=0.5*dt;
    threefourthsdt=0.75*dt;
    for(unsigned int j=start; j < stop; j++) {
      Nu nuk=parent->LinearCoeff(j);
      Nu x=-nuk*dt;
      
      Nu ph2=phi2(x)*dt;
      Nu ph1=x*ph2+dt;
      
      Nu xh=0.5*x;
      Nu ph2h=phi2(xh)*dt;
      Nu ph1H=0.5*xh*ph2h+halfdt;
      
      Nu xi=0.75*x;
      Nu ph2i=phi2(xi)*dt;
      Nu ph1I=0.75*xi*ph2i+threefourthsdt;
      
      coeff0h[j]=-nuk*ph1H+1;
      coeff0i[j]=-nuk*ph1I+1;
      coeff0[j]=-nuk*ph1+1;
      a21[j]=ph1H;
      Nu a32j=9.0/8.0*ph2i+3.0/8.0*ph2h;
      a31[j]=ph1I-a32j;
      a32[j]=a32j;
      
      Nu b2j=ph1/3.0;
      Nu b3j=4.0/3.0*ph2-2.0/9.0*ph1;
      
      b1[j]=ph1-b2j-b3j;
      b2[j]=b2j;
      b3[j]=b3j;
      if(dynamic) {
	B1[j]=ph1-17.0/12.0*ph2;
	B2[j]=0.5*ph2;
	B3[j]=2.0/3.0*ph2;
	B4[j]=0.25*ph2;
      }
    }
  }
  
  void Source(const vector2& Src, const vector2& Y, double t) {
    parent->ExponentialSource(Src,Y,t);
  }
  
  void PSource(const vector2&, const vector2&, double) {}
  void CSource(const vector2&, const vector2&, double) {}
  
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
inline void E_RK3<T>::Predictor()
{
  for(unsigned int j=start; j < stop; j++)	
    y[j]=coeff0h[j]*y0[j]+a21[j]*source0[j];
  Source(Src1,Y,t+halfdt);
  for(unsigned int j=start; j < stop; j++)	
    y[j]=coeff0i[j]*y0[j]+a31[j]*source0[j]+a32[j]*source1[j];
}

template<class T>
inline void E_RK3<T>::Corrector()
{
  Source(Src,Y,t+threefourthsdt);
  
  for(unsigned int j=start; j < stop; j++)	
    y[j]=coeff0[j]*y0[j]+b1[j]*source0[j]+b2[j]*source1[j]+b3[j]*source[j];
  
  if(dynamic) {
    Source(Src3,Y,t+dt);
    for(unsigned int j=start; j < stop; j++)
      if(!Array::Active(errmask) || errmask[j])
	CalcError(y0[j],y[j],coeff0[j]*y0[j]+B1[j]*source0[j]+
		  B2[j]*source1[j]+B3[j]*source[j]+B4[j]*source3[j],y[j]);
  }
}

class RK3C_Exponential : public Exponential {
protected:
  NuVector coeff0h,coeff1h,coeffA,coeffB,coeffC,B0,B1;
public:  
  void Allocator(unsigned int ny, unsigned int start, bool extra=false) {
    Allocate(coeff0h,ny,start,0);
    Allocate(coeff1h,ny,start,0);
    Allocate(coeffA,ny,start,0);
    Allocate(coeffB,ny,start,0);
    Allocate(coeffC,ny,start,0);
    if(extra) {
      Allocate(B0,ny,start,0);
      Allocate(B1,ny,start,0);
    }
  }
  
  virtual Nu LinearCoeff(unsigned int j)=0;
  
  void TimestepDependence(Real dt, bool extra=false) {
    for(unsigned int j=start; j < stop; j++) {
      Nu nuk=LinearCoeff(j);
      Nu x=-nuk*dt;
      Nu xh=0.5*x;
      Real halfdt=0.5*dt;
      Nu ph1h=phi1(xh)*halfdt;
      coeff0h[j]=-nuk*ph1h+1;
      coeff1h[j]=ph1h;
      Nu ph3=phi3(x)*dt;
      Nu ph2=x*ph3+0.5*dt;
      Nu ph1=x*ph2+dt;
      coeff0[j]=-nuk*ph1+1;
      coeff1[j]=ph1;
      coeffA[j]=ph1-3.0*ph2+4.0*ph3;
      coeffB[j]=2.0*(ph2-2.0*ph3);
      coeffC[j]=4.0*ph3-ph2;
      if(extra) {
	B0[j]=ph1-2.0*ph2;
	B1[j]=2.0*ph2;
      }
    }
  }
};

template<class T>
class E_RK3C : public RK3C, public RK3C_Exponential {
protected:
  T *parent;
public:
  E_RK3C(T *parent) : parent(parent) {order=3;}
  
  const char *Name() {return "Third-Order Exponential Classical Runge-Kutta";}
  
  void Allocator() {
    RK3C::Allocator();
    parent->IndexLimits(start,stop,startN,stop0,start0,stopN);
    Exponential::Allocator();
  }
  
  void Allocator(unsigned int ny) {
    RK3C_Exponential::Allocator(ny,start,dynamic);
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
    RK3C::TimestepDependence();
    RK3C_Exponential::TimestepDependence(dt,dynamic);
  }

  inline void Predictor(),Corrector();
  
  virtual void Predictor(unsigned int, unsigned int) {
    E_RK3C<T>::Predictor();
    RK3C::Predictor(startN,stopN);
  }
  
  virtual int Corrector(unsigned int, unsigned int) {
    E_RK3C<T>::Corrector();
    return RK3C::Corrector(startN,stopN);
  }
};

template<class T>
void E_RK3C<T>::Predictor()
{
  for(unsigned int j=start; j < stop; j++)	
    y[j]=coeff0h[j]*y0[j]+coeff1h[j]*source0[j];
  Source(Src1,Y,t+halfdt);
  for(unsigned int j=start; j < stop; j++)	
    y[j]=coeff0[j]*y0[j]+coeff1[j]*(2.0*source1[j]-source0[j]);
}

template<class T>
void E_RK3C<T>::Corrector()
{
  Source(Src,Y,t+dt);
  if(dynamic) {
    for(unsigned int j=start; j < stop; j++) {
      Var val=coeff0[j]*y0[j]+coeffA[j]*source0[j]+2.0*coeffB[j]*source1[j]
	+coeffC[j]*source[j];
      if(!Array::Active(errmask) || errmask[j])
	CalcError(y0[j],val,coeff0[j]*y0[j]+B0[j]*source0[j]+B1[j]*source1[j],
		  val);
      y[j]=val;
    }
  } else {
    for(unsigned int j=start; j < stop; j++)
      y[j]=coeff0[j]*y0[j]+coeffA[j]*source0[j]+2.0*coeffB[j]*source1[j]
	+coeffC[j]*source[j];
  }
}

template<class T>
class E_RK4 : public RK4, public RK3C_Exponential {
protected:
  T *parent;
public:
  E_RK4(T *parent) : parent(parent) {order=4;}
  
  const char *Name() {
    return "Fourth-Order Exponential Classical Runge-Kutta";
  }
  
  void Allocator() {
    RK4::Allocator();
    parent->IndexLimits(start,stop,startN,stop0,start0,stopN);
    Exponential::Allocator();
  }
  
  void Allocator(unsigned int ny) {
    RK3C_Exponential::Allocator(ny,start);
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
    RK3C_Exponential::TimestepDependence(dt);
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
      if(!Array::Active(errmask) || errmask[j]) {
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

#if 0
// Needs to be derived from a RK4HO class
template<class T>
class E_RK4HO : public RK4 {
protected:
  unsigned int start,stop;
  unsigned int startN,stopN;
  T *parent;
  NuVector coeff0,coeff0h;
  NuVector a21;
  NuVector a31,a32;
  NuVector a41,a42;
  NuVector a51,a52,a54;
  NuVector b1,b4,b5;
  NuVector B1,B2,B4;
public:
  E_RK4HO(T *parent) : parent(parent) {order=4;}
  
  const char *Name() {
    return "Fourth-Order Exponential Runge-Kutta (Hochbruck-Ostermann)";
  }
  
  void Allocator() {
    RK4::Allocator();
    if(!dynamic) Alloc0(Src3,source3);
    parent->IndexLimits(start,stop,startN,stop0,start0,stopN);

    unsigned int ny=stop-start;
    Allocate(coeff0,ny,start,0);
    Allocate(coeff0h,ny,start,0);
    Allocate(a21,ny,start,0);
    Allocate(a31,ny,start,0);
    Allocate(a32,ny,start,0);
    Allocate(a41,ny,start,0);
    Allocate(a42,ny,start,0);
    Allocate(a51,ny,start,0);
    Allocate(a52,ny,start,0);
    Allocate(a54,ny,start,0);
    Allocate(b1,ny,start,0);
    Allocate(b4,ny,start,0);
    Allocate(b5,ny,start,0);
    if(dynamic) {
      Allocate(B1,ny,start,0);
      Allocate(B2,ny,start,0);
      Allocate(B4,ny,start,0);
    }
  }
  
  void TimestepDependence() {
    RK4::TimestepDependence();
    for(unsigned int j=start; j < stop; j++) {
      Nu nuk=parent->LinearCoeff(j);
      halfdt=0.5*dt;
      Nu x=-nuk*dt;
      Nu xh=0.5*x;
      Nu ph3h=phi3(xh)*dt;
      Nu ph2h=xh*ph3h+halfdt;
      Nu ph1H=0.5*xh*ph2h+halfdt;
      Nu ph3=phi3(x)*dt;
      Nu ph2=x*ph3+halfdt;
      Nu ph1=x*ph2+dt;
      
      coeff0h[j]=-nuk*ph1H+1;
      coeff0[j]=-nuk*ph1+1;
      a21[j]=ph1H;
      a31[j]=ph1H-ph2h;
      a32[j]=ph2h;
      a41[j]=ph1-2*ph2;
      a42[j]=ph2;
      Nu A52=0.5*(ph2h-ph3h)-ph3+0.25*ph2;
      a51[j]=ph1H-A52-0.25*ph2;
      a52[j]=A52;
      a54[j]=0.25*ph2h-A52;
      b1[j]=ph1-3*ph2+4*ph3;
      b4[j]=4*ph3-ph2;
      b5[j]=4*ph2-8*ph3;
      Real twothirdsdt=2.0/3.0*dt;
      if(dynamic) {
	B1[j]=ph1-3*ph2+twothirdsdt;
	B2[j]=2*ph2-twothirdsdt;
	B4[j]=twothirdsdt-ph2;
      }
    }
  }
  
  void Source(const vector2& Src, const vector2& Y, double t) {
    parent->ExponentialSource(Src,Y,t);
  }
  
  void PSource(const vector2&, const vector2&, double) {}
  void CSource(const vector2&, const vector2&, double) {}
  
  inline void Predictor(),Corrector();
  
  virtual void Predictor(unsigned int, unsigned int) {
    E_RK4HO<T>::Predictor();
    RK4::Predictor(startN,stopN);
  }
  
  virtual int Corrector(unsigned int, unsigned int) {
    E_RK4HO<T>::Corrector();
    return RK4::Corrector(startN,stopN);
  }
};

template<class T>
inline void E_RK4HO<T>::Predictor()
{
  for(unsigned int j=start; j < stop; j++)	
    y[j]=coeff0h[j]*y0[j]+a21[j]*source0[j];
  Source(Src1,Y,t+halfdt);
  for(unsigned int j=start; j < stop; j++)	
    y[j]=coeff0h[j]*y0[j]+a31[j]*source0[j]+a32[j]*source1[j];
  Source(Src2,Y,t+halfdt);
  for(unsigned int j=start; j < stop; j++)	
    y[j]=coeff0[j]*y0[j]+a41[j]*source0[j]+a42[j]*(source1[j]+source2[j]);
  Source(Src3,Y,t+dt);
  for(unsigned int j=start; j < stop; j++)	
    y[j]=coeff0h[j]*y0[j]+a51[j]*source0[j]+a52[j]*(source1[j]+source2[j])
      +a54[j]*source3[j];
}

template<class T>
inline void E_RK4HO<T>::Corrector()
{
  Source(Src,Y,t+halfdt);
  
  if(dynamic) {
    for(unsigned int j=start; j < stop; j++) {
      Var val=coeff0[j]*y0[j]+b1[j]*source0[j]+b4[j]*source3[j]+
	b5[j]*source[j];
      if(!Array::Active(errmask) || errmask[j]) {
	CalcError(y0[j],val,coeff0[j]*y0[j]+B1[j]*source0[j]+
		  B2[j]*(source1[j]+source2[j])+B4[j]*source3[j],val);
      }
      y[j]=val;
    }
  } else {
    for(unsigned int j=start; j < stop; j++)
      y[j]=coeff0[j]*y0[j]+b1[j]*source0[j]+b4[j]*source3[j]+b5[j]*source[j];
  }
}
#endif

template<class T>
void ExponentialIntegrators(Table<IntegratorBase> *t, T *parent) {
  new entry<E_PC<T>,IntegratorBase,T>("E_PC",t,parent);
  new entry<E_RK2<T>,IntegratorBase,T>("E_RK2",t,parent);
  new entry<E_RK3<T>,IntegratorBase,T>("E_RK3",t,parent);
  new entry<E_RK3C<T>,IntegratorBase,T>("E_RK3C",t,parent);
  new entry<E_RK4<T>,IntegratorBase,T>("E_RK4",t,parent);
//  new entry<E_RK4HO<T>,IntegratorBase,T>("E_RK4HO",t,parent);
}

#endif
