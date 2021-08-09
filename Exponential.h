#ifndef __Exponential_h__
#define __Exponential_h__ 1

#include "options.h"
#include "kernel.h"
#include "phi.h"

// Exponential integrators

typedef Array::Array1<Nu>::opt NuVector;

template<class T>
class E_RK : public RK {
protected:  
  T *parent;
  unsigned int start,stop;
  unsigned int startN,stopN;
  Array::Array3<Nu> e;
  Array::Array2<Nu> f,phi0;
public:
  E_RK(T *parent, int order, int nstages, bool fsal=false) : 
    RK(order,nstages,fsal), parent(parent) {}

  void Allocator() {
    RK::Allocator();
    unsigned int stop0,start0; // Unused
    parent->IndexLimits(start,stop,startN,stop0,start0,stopN);
    this->start=start;
    this->stop=stop;
    this->startN=startN;
    this->stopN=stopN;
    unsigned int ny=stop-start;
    phi0.Allocate(Astages,ny,0,(int) start);
    phi0=0.0;
    e.Allocate(Astages,ny,Astages,0,start,0);
    e=0.0;
    f.Allocate(ny,nstages,start,0);
    f=0.0;
  }
  
  void Source(const vector2& Src, const vector2& Y, double t) {
    parent->ExponentialSource(Src,Y,t);
  }
  
  void Stage(unsigned int s, unsigned int start, unsigned int stop) {
    NuVector phi0s=phi0[s];
    Array::Array2<Nu> es=e[s];
#pragma omp parallel for num_threads(threads)
    for(unsigned int j=start; j < stop; j++) {
      Var sum=phi0s[j]*y0[j];
      NuVector esj=es[j];
      for(unsigned int k=0; k <= s; k++)
	sum += esj[k]*vsource[k][j];
      y[j]=sum;
    }
  }
  
  virtual void PStage(unsigned int s) {
    Stage(s,start,stop);
    RK::Stage(s,startN,stopN);
  }
  
  void Predictor(unsigned int, unsigned int) {
    for(unsigned int s=0; s < Astages-1; ++s) {
      PStage(s);
      PredictorSource(s);
    }
  }

  int Corrector(unsigned int, unsigned int) {
    if(dynamic) {
      NuVector phi0s=phi0[Astages-1];
      if(FSAL) {
	Stage(Astages-1,start,stop);
	RK::Corrector(startN,stopN);
#pragma omp parallel for num_threads(threads)
	for(unsigned int j=start; j < stop; j++) {
	  NuVector fj=f[j];
	  Var sum0=y0[j];
	  Var pred=phi0s[j]*sum0;
	  for(unsigned int k=0; k < nstages; k++)
	    pred += fj[k]*vsource[k][j];
	  if(!Array::Active(errmask) || errmask[j])
	    CalcError(sum0,y[j],pred,y[j]);
	}
      } else {
	Array::Array2<Nu> es=e[Astages-1];
#pragma omp parallel for num_threads(threads)
	for(unsigned int j=start; j < stop; j++) {
	  NuVector esj=es[j];
	  NuVector fj=f[j];
	  Var sum0=y0[j];
	  Var sum=phi0s[j]*sum0;
	  Var pred=sum;
	  for(unsigned int k=0; k < Astages; k++) {
	    Var Skj=vsource[k][j];
	    sum += esj[k]*Skj;
	    pred += fj[k]*Skj;
	  }
	  if(!Array::Active(errmask) || errmask[j])
	    CalcError(sum0,sum,pred,sum);
	  y[j]=sum;
	}
	RK::Corrector(startN,stopN);
      }
    } else {
      Stage(Astages-1,start,stop);
      RK::Stage(Astages-1,startN,stopN);
    }
    return 1;
  }
};

template<class T>
class E_Euler : public E_RK<T> {
protected:
public:
  const char *Name() {return "Exponential Euler";}
  E_Euler(T *parent) : E_RK<T>(parent,1,1) {
    RK::allocate();
    this->A[0][0]=1.0;
  }
  
  inline void TimestepDependence() {
    if(this->startN < this->stopN) RK::TimestepDependence();
    for(unsigned int j=this->start; j < this->stop; j++) {
      Nu nuk=this->parent->LinearCoeff(j);
      Nu x=-nuk*this->dt;
      Nu ph1=phi1(x); // (e^x-1)/x
      this->phi0[0][j]=x*ph1+1.0;
      this->e[0][j][0]=ph1*this->dt;
    }
  }
};

template<class T>
class I_Euler : public E_RK<T> {
protected:
public:
  const char *Name() {return "Integrating Factor Euler";}
  I_Euler(T *parent) : E_RK<T>(parent,1,1) {
    RK::allocate();
    this->A[0][0]=1.0;
  }
  
  inline void TimestepDependence() {
    if(this->startN < this->stopN) RK::TimestepDependence();
    for(unsigned int j=this->start; j < this->stop; j++) {
      Nu nuk=this->parent->LinearCoeff(j);
      Nu ph0=exp(-nuk*this->dt);
      this->phi0[0][j]=ph0;
      this->e[0][j][0]=ph0*this->dt;
    }
  }
};

template<class T>
class E_PC : public E_RK<T> {
protected:
public:
  const char *Name() {return "Exponential Predictor-Corrector";}
  
  E_PC(T *parent) : E_RK<T>(parent,2,2) {
    RK::allocate();
    this->A[0][0]=1.0;
    
    this->A[1][0]=0.5;
    this->A[1][1]=0.5;
    
    this->B[0]=1.0;
  }
  
  inline void TimestepDependence() {
    if(this->startN < this->stopN) RK::TimestepDependence();
    for(unsigned int j=this->start; j < this->stop; j++) {
      Nu nuk=this->parent->LinearCoeff(j);
      Nu x=-nuk*this->dt;
      Nu ph1=phi1(x); // (e^x-1)/x
      this->phi0[1][j]=this->phi0[0][j]=x*ph1+1.0;
      ph1 *= this->dt;
      this->e[0][j][0]=ph1;
      this->f[j][0]=ph1;
      this->e[1][j][1]=this->e[1][j][0]=0.5*ph1;
    }
  }
};

template<class T>
class E_RK2 : public E_RK<T> {
protected:
public:
  const char *Name() {return "Second-Order Exponential Runge-Kutta";}
  
  E_RK2(T *parent) : E_RK<T>(parent,2,2) {
    RK::allocate();
    
    this->A[0][0]=0.5;
    this->A[1][1]=1.0;
    
    this->B[0]=1.0;
  }
  
  inline void TimestepDependence() {
    if(this->startN < this->stopN) RK::TimestepDependence();
    for(unsigned int j=this->start; j < this->stop; j++) {
      Nu nuk=this->parent->LinearCoeff(j);
      Nu x=-nuk*this->dt;
      Nu xh=0.5*x;
      Nu ph1H=0.5*phi1(xh)*this->dt;
      this->phi0[0][j]=-nuk*ph1H+1;
      this->e[0][j][0]=ph1H;
      Nu ph2=phi2(x)*this->dt;
      Nu ph1=x*ph2+this->dt;
      this->phi0[1][j]=-nuk*ph1+1;
      this->f[j][0]=ph1;
      this->e[1][j][0]=ph1-2*ph2;
      this->e[1][j][1]=2*ph2;
    }
  }
};

template<class T>
class E_RK3 : public E_RK<T> {
protected:
public:
  const char *Name() {
    return "Third-Order Exponential Bogacki-Shampine Runge-Kutta";
  }
  
  E_RK3(T *parent) : E_RK<T>(parent,3,4,true) {
    RK::allocate();
    
    this->A[0][0]=0.5;
    
    this->A[1][1]=0.75;
    
    this->A[2][0]=2.0/9.0;
    this->A[2][1]=1.0/3.0;
    this->A[2][2]=4.0/9.0;
    
    this->B[0]=7.0/24.0;
    this->B[1]=0.25;
    this->B[2]=1.0/3.0;
    this->B[3]=0.125;
  }
  
  inline void TimestepDependence() {
    if(this->startN < this->stopN) RK::TimestepDependence();
    for(unsigned int j=this->start; j < this->stop; j++) {
      Nu nuk=this->parent->LinearCoeff(j);
      Nu x=-nuk*this->dt;
      Nu xh=0.5*x;
      Nu xi=0.75*x;
      
      Nu ph1=phi1(x)*this->dt;
      Nu ph2=phi2(x)*this->dt;
      
      Nu ph1h=phi1(xh)*this->dt;
      Nu ph2h=phi2(xh)*this->dt;
      
      Nu ph1i=phi1(xi)*this->dt;
      Nu ph2i=phi2(xi)*this->dt;
      
      this->phi0[0][j]=exp(xh);
      this->phi0[1][j]=exp(xi);
      this->phi0[2][j]=exp(x);
      
      this->e[0][j][0]=0.5*ph1h;
      
      Nu a11j=9.0/8.0*ph2i+3.0/8.0*ph2h;
      this->e[1][j][0]=0.75*ph1i-a11j;
      this->e[1][j][1]=a11j;
      
      Nu a21j=ph1/3.0;
      Nu a22j=4.0/3.0*ph2-2.0/9.0*ph1;
      
      this->e[2][j][0]=2.0*a21j-a22j;
      this->e[2][j][1]=a21j;
      this->e[2][j][2]=a22j;
      
      this->f[j][0]=ph1-17.0/12.0*ph2;
      this->f[j][1]=0.5*ph2;
      this->f[j][2]=2.0/3.0*ph2;
      this->f[j][3]=0.25*ph2;
    }
  }
};

template<class T>
class E_RK3BP : public E_RK<T> {
protected:
public:
  const char *Name() {
    return "Third-Order Exponential Bogacki-Shampine Runge-Kutta BP";
  }
  
  E_RK3BP(T *parent) : E_RK<T>(parent,3,4,true) {
    RK::allocate();
    
    this->A[0][0]=0.5;
    
    this->A[1][1]=0.75;
    
    this->A[2][0]=2.0/9.0;
    this->A[2][1]=1.0/3.0;
    this->A[2][2]=4.0/9.0;
    
    this->B[0]=7.0/24.0;
    this->B[1]=0.25;
    this->B[2]=1.0/3.0;
    this->B[3]=0.125;
  }
  
  inline void TimestepDependence() {
    if(this->startN < this->stopN) RK::TimestepDependence();
    for(unsigned int j=this->start; j < this->stop; j++) {
      Nu nuk=this->parent->LinearCoeff(j);
      Nu x=-nuk*this->dt;
      Nu xh=0.5*x;
      Nu xi=0.75*x;
      
      Nu ph1=phi1(x)*this->dt;
      Nu ph2=phi2(x)*this->dt;
      Nu ph3=phi3(x)*this->dt;
      
      Nu ph1h=phi1(xh)*this->dt;
      
      Nu ph1i=phi1(xi)*this->dt;
      Nu ph2i=phi2(xi)*this->dt;
      
      this->phi0[0][j]=exp(xh);
      this->phi0[1][j]=exp(xi);
      this->phi0[2][j]=exp(x);
      
      this->e[0][j][0]=0.5*ph1h;
      
      this->e[1][j][0]=0.75*ph1i-1.5*ph2i;
      this->e[1][j][1]=1.5*ph2i;
      
      this->e[2][j][0]=ph1-2.0*ph2+4.0/3.0*ph3;
      this->e[2][j][1]=2.0*ph2-4.0*ph3;
      this->e[2][j][2]=8.0*ph3/3;
      
      this->f[j][0]=ph1-17.0/12.0*ph2;
      this->f[j][1]=0.5*ph2;
      this->f[j][2]=2.0/3.0*ph2;
      this->f[j][3]=0.25*ph2;
    }
  }
};

template<class T>
void ExponentialIntegrators(Table<IntegratorBase> *t, T *parent) {
  new entry<E_Euler<T>,IntegratorBase,T>("E_Euler",t,parent);
  new entry<I_Euler<T>,IntegratorBase,T>("I_Euler",t,parent);
  new entry<E_PC<T>,IntegratorBase,T>("E_PC",t,parent);
  new entry<E_RK2<T>,IntegratorBase,T>("E_RK2",t,parent);
  new entry<E_RK3<T>,IntegratorBase,T>("E_RK3",t,parent);
  new entry<E_RK3BP<T>,IntegratorBase,T>("E_RK3BP",t,parent);
}

#endif
