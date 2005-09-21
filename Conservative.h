#ifndef __Conservative_h__
#define __Conservative_h__ 1

// Conservative integrators

class Conservative {
protected:
  unsigned int start,stop;
  unsigned int startT,stopT;
  unsigned int startM,stopM;
public:
  virtual ~Conservative() {}
};

template<class T>
class C_PC : public PC, public Conservative {
  T *parent;
public:
  C_PC(T *parent) : parent(parent) {}
  const char *Name() {return "Conservative Predictor-Corrector";}
  
  void Allocator() {
    PC::Allocator();
    parent->IndexLimits(start,stop,startT,stopT,startM,stopM);
  }
  
  inline bool Correct(Real y0, Real& y, Real source0, Real source);
  inline bool Correct(const Complex& y0, Complex& y,
		      const Complex& source0, const Complex& source);
  
  void PSource(const vector2& Src, const vector2& Y, double t) {
    parent->ConservativeSource(Src,Y,t);
  }
  
  void CSource(const vector2& Src, const vector2& Y, double t) {
    parent->NonConservativeSource(Src,Y,t+dt);
  }
  
  inline void Predictor(unsigned int, unsigned int) {
    PC::Predictor(start,stopT);
  }
  
  inline int Corrector();
  
  virtual int Corrector(unsigned int, unsigned int) {
    int rc=C_PC<T>::Corrector();
    // Average the conservative moments with a trapezoidal rule
    if(rc) PC::Corrector(startT,stopM);
    return rc;
  }
};

template<class T>
inline bool C_PC<T>::Correct(Real y0, Real& y, Real source0, Real source)
{
  Real discr=y0*y0+dt*(y0*source0+y*source);
  if(discr >= 0.0) y=sgn(y)*sqrt(discr);
  else {
    if(hybrid) y=y0+0.5*dt*(source0+source);
    else return false;
  }
  return true;
}

template<class T>
inline bool C_PC<T>::Correct(const Complex& y0, Complex& y,
			     const Complex& source0, const Complex& source)
{
  if(!Correct(y0.re,y.re,source0.re,source.re)) return false;
  return Correct(y0.im,y.im,source0.im,source.im);
}

template<class T>
inline int C_PC<T>::Corrector()
{
  parent->ConservativeSource(Src,Y,t+dt);
  if(dynamic) {
    for(unsigned int j=start; j < stop; j++) {
      Var val=y[j];
      if(!Correct(y0[j],val,source0[j],source[j])) return 0;
      if(!Array::Active(errmask) || errmask[j])
	CalcError(y0[j],val,y[j],val);
      y[j]=val;
    }
  } else for(unsigned int j=start; j < stop; j++)
    if(!Correct(y0[j],y[j],source0[j],source[j])) return 0;
  
  return 1;
}

template<class T>
class C_RK2 : public RK2, public Conservative {
  T *parent;
public:
  C_RK2(T *parent) : parent(parent) {}
  const char *Name() {return "Conservative Second-Order Runge-Kutta";}
  
  void Allocator() {
    RK2::Allocator();
    parent->IndexLimits(start,stop,startT,stopT,startM,stopM);
  }
  
  inline bool Correct(Real y0, Real& y, Real source);
  inline bool Correct(const Complex& y0, Complex& y, const Complex& source);
  
  void PSource(const vector2& Src, const vector2& Y, double t) {
    parent->ConservativeSource(Src,Y,t);
  }
  
  void CSource(const vector2& Src, const vector2& Y, double t) {}
  
  inline void Predictor(unsigned int, unsigned int) {
    RK2::Predictor(start,stopT);
  }
  
  inline int Corrector();
  
  virtual int Corrector(unsigned int, unsigned int) {
    int rc=C_RK2<T>::Corrector();
    if(rc) {
      RK2::Corrector(startT,stopT);
      // Average the conservative moments with a trapezoidal rule
      parent->NonConservativeSource(Src,Y,t+dt);
      PC::Corrector(startM,stopM);
    }
    return rc;
  }
};

template<class T>
inline bool C_RK2<T>::Correct(Real y0, Real& y, Real source)
{
  Real temp=dt*source;
  Real discr=y0*y0+2.0*y*temp;
  if(discr >= 0.0) y=sgn(y0+temp)*sqrt(discr);
  else {
    if(hybrid) y=y0+temp;
    else return false;
  }
  return true;
}

template<class T>
inline bool C_RK2<T>::Correct(const Complex& y0, Complex& y,
			      const Complex& source)
{
  if(!Correct(y0.re,y.re,source.re)) return false;
  return Correct(y0.im,y.im,source.im);
}

template<class T>
inline int C_RK2<T>::Corrector()
{
  parent->ConservativeSource(Src,Y,t+halfdt);
  if(dynamic) {
    for(unsigned int j=start; j < stop; j++) {
      Var val=y[j];
      if(!Correct(y0[j],val,source[j])) return 0;
      if(!Array::Active(errmask) || errmask[j])
	CalcError(y0[j],val,y0[j]+dt*source0[j],val);
      y[j]=val;
    }
  } else for(unsigned int j=start; j < stop; j++) {
      if(!Correct(y0[j],y[j],source[j])) return 0;
  }
  
  return 1;
}

template<class T>
class C_RK4 : public RK4, public Conservative {
  T *parent;
  double sixthdt2;
public:
  C_RK4(T *parent) : parent(parent) {}
  const char *Name() {return "Conservative Fourth-Order Runge-Kutta";}
  
  void Allocator() {
    RK4::Allocator();
    parent->IndexLimits(start,stop,startT,stopT,startM,stopM);
  }
  
  void TimestepDependence();
  
  inline bool Correct(Real y0, Real& y, Real source0, Real source1,
		      Real source2, Real source);
  inline bool Correct(const Complex& y0, Complex& y, const Complex& source0,
		      const Complex& source1, const Complex& source2, 
		      const Complex& source);
  
  void PSource(const vector2& Src, const vector2& Y, double t) {
    parent->ConservativeSource(Src,Y,t);
  }
  
  void CSource(const vector2& Src, const vector2& Y, double t) {}
  
  inline void Predictor(unsigned int, unsigned int) {
    RK4::Predictor(start,stopT);
  }
  
  inline int Corrector();
  
  virtual int Corrector(unsigned int, unsigned int) {
    int rc=C_RK4<T>::Corrector();
    if(rc) {
      RK4::Corrector(startT,stopT);
      // Average the conservative moments with a trapezoidal rule
      parent->NonConservativeSource(Src,Y,t+dt);
      PC::Corrector(startM,stopM);
    }
    return rc;
  }
};

template<class T>
void C_RK4<T>::TimestepDependence()
{
  PC::TimestepDependence();
  RK4::TimestepDependence();
  sixthdt2=sixthdt*dt;
}

template<class T>
inline bool C_RK4<T>::Correct(Real y0, Real& y, Real Source0, Real Source1,
			      Real Source2, Real Source)
{			
  Real temp=sixthdt*(Source0+2.0*(Source1+Source2)+Source);
  Real discr=y0*y0+2.0*(y0*temp+sixthdt2*(Source1*(Source0+Source2)+
					  Source2*Source));
  if(discr >= 0.0) {
    y=sgn(y0+temp)*sqrt(discr);
  } else {
    if(hybrid) y=y0+temp;
    else return false;
  }
  return true;
}


template<class T>
inline bool C_RK4<T>::Correct(const Complex& y0, Complex& y, 
			      const Complex& source0, const Complex& source1,
			      const Complex& source2, const Complex& source)
{			
  if(!Correct(y0.re,y.re,source0.re,source1.re,source2.re,source.re))
    return false;
  return Correct(y0.im,y.im,source0.im,source1.im,source2.im,source.im);
}

template<class T>
inline int C_RK4<T>::Corrector()
{
  parent->ConservativeSource(Src,Y,t+dt);
  if(dynamic) {
    for(unsigned int j=start; j < stop; j++)
      y[j]=y0[j]+dt*(2.0*source1[j]-source0[j]);
    parent->ConservativeSource(Src3,Y,t+dt);
    for(unsigned int j=start; j < stop; j++) {
      Var val;
      if(!Correct(y0[j],val,source0[j],source1[j],source2[j],source[j]))
	return 0;
      if(!Array::Active(errmask) || errmask[j]) 
	CalcError(y0[j],val,
		  y0[j]+sixthdt*(source0[j]+4.0*source1[j]+source3[j]),val);
      y[j]=val;
    }
  } else for(unsigned int j=start; j < stop; j++) {
    if(!Correct(y0[j],y[j],source0[j],source1[j],source2[j],source[j]))
      return 0;
  }
  return 1;
}

template<class T>
class C_RK5 : public RK5, public Conservative {
  T *parent;
  vector y2,y3,y4;
  vector2 Y2,Y3,Y4;
public:
  C_RK5(T *parent) : parent(parent) {}
  const char *Name() {return "Conservative Fifth-Order Runge-Kutta";}
  
  void Allocator() {
    RK5::Allocator();
    Alloc(Y2,y2);
    Alloc(Y3,y3);
    Alloc(Y4,y4);
    parent->IndexLimits(start,stop,startT,stopT,startM,stopM);
  }
  
  void TimestepDependence();
  
  inline bool Correct(Real y0, Real& y2, Real &y3,
		      Real y4, Real& y,
		      Real source0, Real source2, 
		      Real source3, Real source4,
		      Real source, Real& pred, Real& corr);
  inline bool Correct(const Complex& y0, Complex& y2, Complex &y3,
		      const Complex& y4, Complex& y,
		      const Complex& source0, const Complex& source2, 
		      const Complex& source3, const Complex& source4,
		      const Complex& source, Complex& pred, Complex& corr);
  
  void PSource(const vector2& Src, const vector2& Y, double t) {
    parent->ConservativeSource(Src,Y,t);
  }
  
  void CSource(const vector2& Src, const vector2& Y, double t) {}
  
  inline void Predictor(unsigned int start, unsigned int stop);
  
  inline int Corrector();
  
  virtual int Corrector(unsigned int, unsigned int) {
    int rc=C_RK5<T>::Corrector();
    if(rc) {
      RK5::Corrector(startT,stopT);      
      // Average the conservative moments with a trapezoidal rule
      parent->NonConservativeSource(Src,Y,t+dt);
      PC::Corrector(startM,stopM);
    }
    return rc;
  }
};

template<class T>
void C_RK5<T>::TimestepDependence()
{
  PC::TimestepDependence();
  RK5::TimestepDependence();
}

template<class T>
void C_RK5<T>::Predictor(unsigned int start, unsigned int stop)
{
  for(unsigned int j=start; j < stopT; j++)
    y[j]=y0[j]+b10*source0[j];
  PSource(Src,Y,t+a1);
  for(unsigned int j=start; j < stopT; j++)
    y2[j]=y0[j]+b20*source0[j]+b21*source[j];
  PSource(Src2,Y2,t+a2);
  for(unsigned int j=start; j < stopT; j++)
    y3[j]=y0[j]+b30*source0[j]+b31*source[j]+b32*source2[j];
  PSource(Src3,Y3,t+a3);
  for(unsigned int j=start; j < stopT; j++)
    y4[j]=y0[j]+b40*source0[j]+b41*source[j]+b42*source2[j]+b43*source3[j];
  PSource(Src4,Y4,t+a4);
  for(unsigned int j=start; j < stopT; j++) 
    y[j]=y0[j]+b50*source0[j]+b51*source[j]+b52*source2[j]+b53*source3[j]+
      b54*source4[j];
}

template<class T>
inline bool C_RK5<T>::Correct(Real y0, Real& y2, Real& y3,
			      Real y4, Real& y,
			      Real source0, Real source2, 
			      Real source3, Real source4,
			      Real source, Real& pred, Real& discr)
{
  pred=y0*y0+2.0*(d0*y0*source0+d2*y2*source2+d3*y3*source3+d4*y4*source4+
		  d5*y*source);
  discr=y0*y0+2.0*(c0*y0*source0+c2*y2*source2+c3*y3*source3+c5*y*source);
	
  Real temp=y0+c0*source0+c2*source2+c3*source3+c5*source;
  if(discr >= 0.0)
    y=sgn(temp)*sqrt(discr);
  else {
    if(hybrid) {
      y=temp;
      discr=y*y;
    } else return false;
  }
  return true;
}

template<class T>
inline bool C_RK5<T>::Correct(const Complex& y0, Complex& y2, Complex& y3,
			      const Complex& y4, Complex& y,
			      const Complex& source0, const Complex& source2, 
			      const Complex& source3, const Complex& source4,
			      const Complex& source,
			      Complex& pred, Complex& corr)
{
  if(!Correct(y0.re,y2.re,y3.re,y4.re,y.re,source0.re,source2.re,
	      source3.re,source4.re,source.re,pred.re,corr.re)) return false;
  return Correct(y0.im,y2.im,y3.im,y4.im,y.im,source0.im,source2.im,
		 source3.im,source4.im,source.im,pred.im,corr.im);
}

template<class T>
inline int C_RK5<T>::Corrector()
{
  parent->ConservativeSource(Src,Y,t+a5);
  if(dynamic) {
    for(unsigned int j=start; j < stop; j++) {
      Var pred,corr;
      Var val=y[j];
      if(!Correct(y0[j],y2[j],y3[j],y4[j],val,source0[j],source2[j],
		  source3[j],source4[j],source[j],pred,corr)) return 0;
      if(!Array::Active(errmask) || errmask[j])
//	CalcError(y0[j]*y0[j],corr,pred,corr);
	CalcError(y0[j]*y0[j],val,y[j],val);
      y[j]=val;
    }
  } else {
    for(unsigned int j=start; j < stop; j++) {
      Var pred,corr;
      if(!Correct(y0[j],y2[j],y3[j],y4[j],y[j],source0[j],source2[j],
		  source3[j],source4[j],source[j],pred,corr)) return 0;
    }
  }
  return 1;
}

template<class T>
void ConservativeIntegrators(Table<IntegratorBase> *t, T *parent) {
  new entry<C_PC<T>,IntegratorBase,T>("C_PC",t,parent);
  new entry<C_RK2<T>,IntegratorBase,T>("C_RK2",t,parent);
  new entry<C_RK4<T>,IntegratorBase,T>("C_RK4",t,parent);
  new entry<C_RK5<T>,IntegratorBase,T>("C_RK5",t,parent);
}

#endif
