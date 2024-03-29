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
  size_t start,stop;
  size_t startN,stopN;
  Array::Array3<Nu> e;
  Array::Array2<Nu> f,phi0;
public:
  E_RK(T *parent, size_t order, size_t nstages, bool fsal=false) :
    RK(order,nstages,fsal), parent(parent) {}

  bool Exponential() {return true;}

  void Allocator() {
    RK::Allocator();
    size_t stop0,start0; // Unused
    parent->IndexLimits(start,stop,startN,stop0,start0,stopN);
    this->start=start;
    this->stop=stop;
    this->startN=startN;
    this->stopN=stopN;
    size_t ny=stop-start;
    phi0.Allocate(Astages,ny,0,start);
    phi0=0.0;
    e.Allocate(Astages,ny,Astages,0,start,0);
    e=0.0;
    f.Allocate(ny,nstages,start,0);
    f=0.0;
  }

  void Source(const vector2& Src, const vector2& Y, double t) {
    parent->ExponentialSource(Src,Y,t);
  }

  void Stage(size_t s, size_t start, size_t stop) {
    NuVector phi0s=phi0[s];
    Array::Array2<Nu> es=e[s];
    PARALLELIF(
      stop-start > threshold,
      for(size_t j=start; j < stop; j++) {
        Var sum=phi0s[j]*y0[j];
        NuVector esj=es[j];
        for(size_t k=0; k <= s; k++)
          sum += esj[k]*vsource[k][j];
        y[j]=sum;
      });
  }

  virtual void PStage(size_t s) {
    Stage(s,start,stop);
    RK::Stage(s,startN,stopN);
  }

  void Predictor(size_t, size_t) {
    for(size_t s=0; s < Astages-1; ++s) {
      PStage(s);
      PredictorSource(s);
    }
  }

  int Corrector(size_t, size_t) {
    if(dynamic) {
      NuVector phi0s=phi0[Astages-1];
      if(FSAL) {
	Stage(Astages-1,start,stop);
	RK::Corrector(startN,stopN);
        PARALLELIF(
          stop-start > threshold,
          for(size_t j=start; j < stop; j++) {
            NuVector fj=f[j];
            Var sum0=y0[j];
            Var pred=phi0s[j]*sum0;
            for(size_t k=0; k < nstages; k++)
              pred += fj[k]*vsource[k][j];
            if(!Array::Active(errmask) || errmask[j])
              CalcError(sum0,y[j],pred,y[j]);
          });
      } else {
	Array::Array2<Nu> es=e[Astages-1];
        PARALLELIF(
          stop-start > threshold,
          for(size_t j=start; j < stop; j++) {
            NuVector esj=es[j];
            NuVector fj=f[j];
            Var sum0=y0[j];
            Var sum=phi0s[j]*sum0;
            Var pred=sum;
            for(size_t k=0; k < Astages; k++) {
              Var Skj=vsource[k][j];
              sum += esj[k]*Skj;
              pred += fj[k]*Skj;
            }
            if(!Array::Active(errmask) || errmask[j])
              CalcError(sum0,sum,pred,sum);
            y[j]=sum;
          });
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
  E_Euler(T *parent) : E_RK<T>(parent,1,2) {
    RK::allocate();
    this->A[0][0]=1.0;
  }

  inline void TimestepDependence() {
    if(this->startN < this->stopN) RK::TimestepDependence();
    PARALLELIF(
      this->stop-this->start > threshold,
      for(size_t j=this->start; j < this->stop; j++) {
        Nu nuk=this->parent->LinearCoeff(j);
        Nu x=-nuk*this->dt;
        Nu ph1=phi1(x); // (e^x-1)/x
        this->phi0[0][j]=x*ph1+1.0;
        this->e[0][j][0]=ph1*this->dt;
      });
  }
};

template<class T>
class I_Euler : public E_RK<T> {
protected:
public:
  const char *Name() {return "Integrating Factor Euler";}
  I_Euler(T *parent) : E_RK<T>(parent,1,2) {
    RK::allocate();
    this->A[0][0]=1.0;
  }

  inline void TimestepDependence() {
    if(this->startN < this->stopN) RK::TimestepDependence();
    PARALLELIF(
      this->stop-this->start > threshold,
      for(size_t j=this->start; j < this->stop; j++) {
        Nu nuk=this->parent->LinearCoeff(j);
        Nu ph0=exp(-nuk*this->dt);
        this->phi0[0][j]=ph0;
        this->e[0][j][0]=ph0*this->dt;
      });
  }
};

template<class T>
class E_PC : public E_RK<T> {
protected:
public:
  const char *Name() {return "Exponential Predictor-Corrector";}

  E_PC(T *parent) : E_RK<T>(parent,2,3) {
    RK::allocate();

    this->A[0][0]=1.0;

    this->A[1][0]=0.5;
    this->A[1][1]=0.5;

    this->B[0]=1.0;
  }

  inline void TimestepDependence() {
    if(this->startN < this->stopN) RK::TimestepDependence();
    PARALLELIF(
      this->stop-this->start > threshold,
      for(size_t j=this->start; j < this->stop; j++) {
        Nu nuk=this->parent->LinearCoeff(j);
        Nu x=-nuk*this->dt;
        Nu ph1=phi1(x); // (e^x-1)/x
        this->phi0[1][j]=this->phi0[0][j]=x*ph1+1.0;
        ph1 *= this->dt;
        this->e[0][j][0]=ph1;
        this->f[j][0]=ph1;
        this->e[1][j][1]=this->e[1][j][0]=0.5*ph1;
      });
  }
};

template<class T>
class E_RK2 : public E_RK<T> {
protected:
public:
  const char *Name() {return "Second-Order Exponential Runge-Kutta";}

  E_RK2(T *parent) : E_RK<T>(parent,2,3) {
    RK::allocate();

    this->A[0][0]=0.5;
    this->A[1][1]=1.0;

    this->B[0]=1.0;
  }

  inline void TimestepDependence() {
    if(this->startN < this->stopN) RK::TimestepDependence();
    PARALLELIF(
      this->stop-this->start > threshold,
      for(size_t j=this->start; j < this->stop; j++) {
        Nu nuk=this->parent->LinearCoeff(j);
        Nu x=-nuk*this->dt;
        Nu xh=0.5*x;
        Nu ph1H=0.5*phi1(xh)*this->dt;
        this->phi0[0][j]=-nuk*ph1H+1.0;
        this->e[0][j][0]=ph1H;
        Nu ph2=phi2(x)*this->dt;
        Nu ph1=x*ph2+this->dt;
        this->phi0[1][j]=-nuk*ph1+1.0;
        this->f[j][0]=ph1;
        this->e[1][j][0]=ph1-2.0*ph2;
        this->e[1][j][1]=2.0*ph2;
      });
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
    PARALLELIF(
      this->stop-this->start > threshold,
      for(size_t j=this->start; j < this->stop; j++) {
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
      });
  }
};

template<class T>
class E_RK32ZB : public E_RK<T> {
protected:
public:
  const char *Name() {
    return "Robust Third-Order Exponential Runge-Kutta";
  }

  E_RK32ZB(T *parent) : E_RK<T>(parent,3,4,true) {
    RK::allocate();

    this->A[0][0]=0.5;

    this->A[1][1]=0.75;

    this->A[2][0]=2.0/9.0;
    this->A[2][1]=1.0/3.0;
    this->A[2][2]=4.0/9.0;

    this->B[0]=2101.0/2520.0;
    this->B[1]=-179.0/252.0;
    this->B[2]=3.0/35.0;
    this->B[3]=1993.0/2520.0;
  }

  inline void TimestepDependence() {
    if(this->startN < this->stopN) RK::TimestepDependence();
    PARALLELIF(
      this->stop-this->start > threshold,
      for(size_t j=this->start; j < this->stop; j++) {
        Nu nuk=this->parent->LinearCoeff(j);
        Nu x=-nuk*this->dt;
        Nu x1=0.5*x;
        Nu x2=0.75*x;

        this->phi0[0][j]=exp(x1);
        this->phi0[1][j]=exp(x2);
        this->phi0[2][j]=exp(x);

        Nu w1=phi1(x)*this->dt;
        Nu w2=phi2(x)*this->dt;
        Nu w3=phi3(x)*this->dt;

        Nu w1c1=phi1(x1)*this->dt;
        Nu w2c1=phi2(x1)*this->dt;
        Nu w3c1=phi3(x1)*this->dt;

        Nu w1c2=phi1(x2)*this->dt;
        Nu w2c2=phi2(x2)*this->dt;

        this->e[0][j][0]=0.5*w1c1;

        Nu a11j=9.0/8.0*w2c2+3.0/8.0*w2c1;
        this->e[1][j][0]=0.75*w1c2-a11j;
        this->e[1][j][1]=a11j;

        Nu a21=0.75*w2-0.25*w3;
        Nu a22=5.0/6.0*w2+1.0/6.0*w3;

        this->e[2][j][0]=w1-a21-a22;
        this->e[2][j][1]=a21;
        this->e[2][j][2]=a22;

        this->f[j][0]=29.0/18.0*w1+7.0/6.0*w1c2+9.0/14.0*w1c1+0.75*w2+2.0/7.0*w2c2+1.0/12.0*w2c1-8083.0/420.0*w3+11.0/30.0*w3c1;
        this->f[j][1]=-1.0/9.0*w1-1.0/6.0*w1c2-0.5*w2-1.0/7.0*w2c2-1.0/3.0*w2c1+1.0/6.0*w3+1.0/6.0*w3c1;
        this->f[j][2]=2.0/3.0*w1-0.5*w1c2-1.0/7.0*w1c1+1.0/3.0*w2-1.0/7.0*w2c2-0.2*w3c1;
        this->f[j][3]=-7.0/6.0*w1-0.5*w1c2-0.5*w1c1-7.0/12.0*w2+0.25*w2c1+2671.0/140.0*w3-1.0/3.0*w3c1;
      });
  }
};

template<class T>
class E_RK43ZB : public E_RK<T> {
protected:
public:
  const char *Name() {
    return "Robust Fourth-Order Five-Stage Embedded Exponential Runge-Kutta";
  }

  E_RK43ZB(T *parent) : E_RK<T>(parent,4,6) {
    RK::allocate();

    this->A[0][0]=1.0/6.0;

    this->A[1][0]=-0.5;
    this->A[1][1]=1.0;

    this->A[2][0]=-2.5;
    this->A[2][1]=3.5;
    this->A[2][2]=-0.5;

    this->A[3][0]=1.0;
    this->A[3][1]=-1.5;
    this->A[3][2]=2.0;
    this->A[3][3]=-0.5;

    this->A[4][0]=1.0/6.0;
    this->A[4][1]=0.0;
    this->A[4][2]=5.0/6.0;
    this->A[4][3]=-1.0/6.0;
    this->A[4][4]=1.0/6.0;

    for(size_t i=0; i < 4; ++i)
      this->B[i]=this->A[3][i];
    this->B[4]=0.0;
  }

  inline void TimestepDependence() {
    const double sixth=1.0/6.0;
    if(this->startN < this->stopN)
      RK::TimestepDependence();
    PARALLELIF(
      this->stop-this->start > threshold,
      for(size_t j=this->start; j < this->stop; j++) {
        Nu nuk=this->parent->LinearCoeff(j);
        Nu x=-nuk*this->dt;
        Nu x1=sixth*x;
        Nu x2=0.5*x;

        // Automate this with expfactors?
        this->phi0[0][j]=exp(x1);
        Nu e2=exp(x2);
        this->phi0[1][j]=e2;
        this->phi0[2][j]=e2;
        Nu e=exp(x);
        this->phi0[3][j]=e;
        this->phi0[4][j]=e;

        Nu w1_1=phi1(x1)*this->dt;
        Nu w2_1=phi2(x1)*this->dt;

        Nu w1_2=phi1(x2)*this->dt;
        Nu w2_2=phi2(x2)*this->dt;
        Nu w3_2=phi3(x2)*this->dt;

        Nu w1=phi1(x)*this->dt;
        Nu w2=phi2(x)*this->dt;
        Nu w3=phi3(x)*this->dt;

        this->e[0][j][0]=sixth*w1_1;

        Nu a1_1=1.5*w2_2+0.5*w2_1;

        this->e[1][j][0]=0.5*w1_2-a1_1;
        this->e[1][j][1]=a1_1;

        Nu a2_1=19.0/60.0*w1+0.5*w1_2+0.5*w1_1+2.0*w2_2+13.0/6.0*w2_1+0.6*w3_2;
        Nu a2_2=-19.0/180.0*w1-sixth*(w1_2+w1_1+w2_2)+1.0/9.0*w2_1-0.2*w3_2;
        this->e[2][j][0]=0.5*w1_2-a2_1-a2_2;
        this->e[2][j][1]=a2_1;
        this->e[2][j][2]=a2_2;

        Nu a3_3=w2+w2_2-6.0*w3-3.0*w3_2;
        Nu a3_1=3.0*w2-4.5*w2_2-2.5*w2_1+6.0*a3_3+a2_1;
        Nu a3_2=6.0*w3+3.0*w3_2-2.0*a3_3+a2_2;

        // Low-order approximation
        this->f[j][0]=w1-a3_1-a3_2-a3_3;
        this->f[j][1]=a3_1;
        this->f[j][2]=a3_2;
        this->f[j][3]=a3_3;
        this->f[j][4]=0;

        for(size_t i=0; i < 4; ++i)
          this->e[3][j][i]=this->f[j][i];

        // High-order approximation
        this->e[4][j][0]=w1-(67.0/9.0)*w2+(52.0/3.0)*w3;
        this->e[4][j][1]=8.0*w2-24.0*w3;
        this->e[4][j][2]=(26.0/3.0)*w3-(11.0/9.0)*w2;
        this->e[4][j][3]=(7.0/9.0)*w2-(10.0/3.0)*w3;
        this->e[4][j][4]=(4.0/3.0)*w3-(1.0/9.0)*w2;
      });
  }
};

template<class T>
void ExponentialIntegrators(Table<IntegratorBase> *t, T *parent) {
  new entry<E_Euler<T>,IntegratorBase,T>("E_Euler",t,parent);
  new entry<I_Euler<T>,IntegratorBase,T>("I_Euler",t,parent);
  new entry<E_PC<T>,IntegratorBase,T>("E_PC",t,parent);
  new entry<E_RK2<T>,IntegratorBase,T>("E_RK2",t,parent);
  new entry<E_RK3<T>,IntegratorBase,T>("E_RK3",t,parent);
  new entry<E_RK32ZB<T>,IntegratorBase,T>("E_RK32ZB",t,parent);
  new entry<E_RK43ZB<T>,IntegratorBase,T>("E_RK43ZB",t,parent);
}

#endif
