#ifndef __Conservative_h__
#define __Conservative_h__ 1

// Conservative integrators

template<class T>
class C_RK : public RK {
protected:
  T *parent;
  Array::Array2<Var> ys;
  unsigned int start,stop;
  unsigned int startT,stopT;
  unsigned int startM,stopM;
public:
  C_RK(T *parent, int Order, int nstages, bool FSAL=false) :
    RK(Order,nstages,FSAL), parent(parent) {}

  bool isConservative() {return true;}

  void Source(const vector2& Src, const vector2& Y, double t) {
    parent->ConservativeSource(Src,Y,t);
  }

  void Allocator() {
    RK::Allocator();
    parent->IndexLimits(start,stop,startT,stopT,startM,stopM);
    unsigned int ny=stop-start;
    ys.Allocate(nstages-2,ny,0,start);
  }

  virtual void Predictor(unsigned int, unsigned int) {
    unsigned int laststage=Astages-1;
    for(unsigned int s=0; s < laststage; ++s) {
      RK::Stage(s,start,stop);
      if(s < Astages-2) {
        Array::Array1<Var>::opt yss=ys[s];
#pragma omp parallel for num_threads(threads)
	for(unsigned int j=start; j < stop; j++)
	  yss[j]=y[j];
      }
      RK::Stage(s,startT,stopT);
      PredictorSource(s);
    }
  }

  inline Real signedsquareRoot(Real x, Real ref) {
    return sgn(ref)*sqrt(x);
  }

  inline Complex signedsquareRoot(const Complex& x, const Complex & ref) {
    return Complex(signedsquareRoot(x.re,ref.re),signedsquareRoot(x.im,ref.im));
  }

  inline bool nonnegative(Real x) {
    return x >= 0.0;
  }

  inline bool nonnegative(const Complex& x) {
    return x.re >= 0.0 && x.im >= 0.0;
  }

  int Corrector(unsigned int, unsigned int) {
    unsigned int laststage=Astages-1;
    if(dynamic) {
      if(FSAL) {
	msg(ERROR,"Sorry; FSAL not yet implemented");
      } else {
	rvector as=a[laststage];
        bool cont=true;
#pragma omp parallel for num_threads(threads)
	for(unsigned int j=start; j < stop; j++) {
          if(cont) {
            Var y0j=y0[j];
            Var Skj=vsource[0][j];
            Var temp=y0j+as[0]*Skj;
            Var yS=product(y0j,Skj);
            Var discr=as[0]*yS;
            unsigned int k;
            for(k=1; k < laststage; k++) {
              Var Skj=vsource[k][j];
              temp += as[k]*Skj;
              Var yS=product(ys[k-1][j],Skj);
              discr += as[k]*yS;
            }
            Skj=vsource[k][j];
            temp += as[k]*Skj;
            yS=product(y[j],Skj);
            discr=product(y0j,y0j)+2.0*(discr+as[k]*yS);
            Var val;
            if(nonnegative(discr))
              val=signedsquareRoot(discr,temp);
            else {
              if(hybrid) {
                val=temp;
                discr=product(temp,temp);
              } else {cont=false; continue;}
            }
            if(!Array::Active(errmask) || errmask[j])
              CalcError(y0j,val,y[j],val);
            y[j]=val;
          }
	}
        if(!cont) return 0;
      }
    } else {
      rvector as=a[laststage];
        bool cont=true;
#pragma omp parallel for num_threads(threads)
      for(unsigned int j=start; j < stop; j++) {
        if(cont) {
          Var y0j=y0[j];
          Var Skj=vsource[0][j];
          Var temp=y0j+as[0]*Skj;
          Var discr=as[0]*product(y0j,Skj);
          unsigned int k;
          for(k=1; k < laststage; k++) {
            Var Skj=vsource[k][j];
            temp += as[k]*Skj;
            discr += as[k]*product(ys[k-1][j],Skj);
          }
          Skj=vsource[k][j];
          temp += as[k]*Skj;
          discr=product(y0j,y0j)+2.0*(discr+as[k]*product(y[j],Skj));
          if(nonnegative(discr))
            y[j]=signedsquareRoot(discr,temp);
          else {
            if(hybrid)
              y[j]=temp;
            cont=false;
          }
        }
      }
      if(!cont) return 0;
    }

    RK::Corrector(startT,stopT);

    // Average the conservative moments with a trapezoidal rule
    parent->NonConservativeSource(vSrc[0],Y0,t);
    parent->NonConservativeSource(vSrc[1],Y,t+this->dt);
    ::vector source0=vsource[0];
    ::vector source=vsource[1];
    double halfdt=0.5*this->dt;

    if(dynamic) {
#pragma omp parallel for num_threads(threads)
      for(unsigned int j=startM; j < stopM; j++) {
	Var val=y0[j]+halfdt*(source0[j]+source[j]);
	if(!Array::Active(this->errmask) || this->errmask[j])
	  CalcError(y0[j],val,y0[j]+this->dt*source0[j],val);
	y[j]=val;
      }
    } else {
#pragma omp parallel for num_threads(threads)
      for(unsigned int j=startM; j < stopM; j++)
	y[j]=y0[j]+halfdt*(source0[j]+source[j]);
    }
    return 1;
  }

};

template<class T>
class C_PC : public C_RK<T> {
protected:
public:
  const char *Name() {return "Conservative Predictor-Corrector";}

  C_PC(T *parent) : C_RK<T>(parent,2,2) {
    RK::allocate();
    this->A[0][0]=1.0;

    this->A[1][0]=0.5;
    this->A[1][1]=0.5;

    this->B[0]=1.0;
  }
};

template<class T>
class C_RK2 : public C_RK<T> {
protected:
public:
  const char *Name() {return "Conservative Second-Order Runge-Kutta";}

  C_RK2(T *parent) : C_RK<T>(parent,2,2) {
    RK::allocate();

    this->A[0][0]=0.5;

    this->A[1][1]=1.0;

    this->B[0]=1.0;
  }
};

template<class T>
class C_RK4 : public C_RK<T> {
protected:
public:
  const char *Name() {return "Conservative Fourth-Order Runge-Kutta";}

  C_RK4(T *parent) : C_RK<T>(parent,4,5) {
    RK::allocate();

    this->A[0][0]=0.5;

    this->A[1][1]=0.5;

    this->A[2][2]=1.0;

    this->A[3][0]=-1.0;
    this->A[3][1]=2.0;

    this->A[4][0]=1.0/6.0;
    this->A[4][1]=1.0/3.0;
    this->A[4][2]=1.0/3.0;
    this->A[4][3]=1.0/6.0;

    this->B[0]=1.0/6.0;
    this->B[1]=2.0/3.0;
    this->B[4]=1.0/6.0;
  }
};

template<class T>
class C_RK5 : public C_RK<T> {
protected:
public:
  const char *Name() {return "Conservative Fifth-Order Runge-Kutta";}

  C_RK5(T *parent) : C_RK<T>(parent,5,6) {
    RK::allocate();

    this->A[0][0]=0.2;
    this->A[1][0]=3.0/40.0; this->A[1][1]=9.0/40.0;

    this->A[2][0]=0.3; this->A[2][1]=-0.9; this->A[2][2]=1.2;

    this->A[3][0]=-11.0/54.0; this->A[3][1]=2.5; this->A[3][2]=-70.0/27.0;
    this->A[3][3]=35.0/27.0;

    this->A[4][0]=1631.0/55296.0; this->A[4][1]=175.0/512.0;
    this->A[4][2]=575.0/13824.0;
    this->A[4][3]=44275.0/110592.0; this->A[4][4]=253.0/4096.0;

    this->A[5][0]=37.0/378.0; this->A[5][2]=250.0/621.0;
    this->A[5][3]=125.0/594.0;
    this->A[5][5]=512.0/1771.0;

    this->B[0]=2825.0/27648.0; this->B[2]=18575.0/48384.0;
    this->B[3]=13525.0/55296.0;
    this->B[4]=277.0/14336.0; this->B[5]=0.25;
  }
};
template<class T>
void ConservativeIntegrators(Table<IntegratorBase> *t, T *parent) {
  new entry<C_PC<T>,IntegratorBase,T>("C_PC",t,parent);
  new entry<C_RK2<T>,IntegratorBase,T>("C_RK2",t,parent);
  new entry<C_RK4<T>,IntegratorBase,T>("C_RK4",t,parent);
  new entry<C_RK5<T>,IntegratorBase,T>("C_RK5",t,parent);
}

#endif
