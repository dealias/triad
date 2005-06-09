#ifndef __Conservative_h__
#define __Conservative_h__ 1

// Conservative integrators

class CorrectC_PC {
public:
  bool Correct(Real y0, Real& y,
	       Real source0, Real source, double dt);
  bool Correct(const Complex& y0, Complex& y,
	       const Complex& source0, const Complex& source, double dt);
};

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

class Conservative {
protected:
  unsigned int start,stop;
  unsigned int startN,stopN;
public:
  virtual ~Conservative() {}
};

template<class T>
class C_PC : public PC, public CorrectC_PC, public Conservative {
  T *parent;
public:
  C_PC(T *parent) : parent(parent) {}
  const char *Name() {return "Conservative Predictor-Corrector";}
  
  void Allocator() {
    PC::Allocator();
    parent->IndexLimits(start,stop,startN,stopN);
  }
  
// Average derived conserved quantities with trapezoidal rule
  void CSource(const vector2& Src, const vector2& Y, double t) {
    parent->NonConservativeSource(Src,Y,t);
  }
  
  inline int Corrector();
  
  virtual int Corrector(unsigned int, unsigned int) {
    int rc=C_PC<T>::Corrector();
    if(rc) rc=PC::Corrector(startN,stopN);
    return rc;
  }
};

template<class T>
inline int C_PC<T>::Corrector()
{
  parent->ConservativeSource(Src,Y,t+dt);
  if(dynamic) {
    for(unsigned int j=start; j < stop; j++) {
      Var pred=y[j];
      if(!Correct(y0[j],y[j],source0[j],source[j],dt)) return 0;
      if(!Array::Active(errmask) || errmask[j])
	CalcError(y0[j],y[j],pred,y[j]);
    }
  } else for(unsigned int j=start; j < stop; j++)
    if(!Correct(y0[j],y[j],source0[j],source[j],dt)) return 0;
  
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
    parent->IndexLimits(start,stop,startN,stopN);
  }
  
  inline bool Correct(Real y0, Real& y, Real source);
  inline bool Correct(const Complex& y0, Complex& y, const Complex& source);
  
  void CSource(const vector2& Src, const vector2& Y, double t) {
    parent->NonConservativeSource(Src,Y,t);
  }
  
  inline int Corrector();
  
  virtual int Corrector(unsigned int, unsigned int) {
    int rc=C_RK2<T>::Corrector();
    if(rc) rc=PC::Corrector(startN,stopN);
    return rc;
  }
};

template<class T>
inline bool C_RK2<T>::Correct(Real Y0, Real& Y, Real Source)
{
  Real temp=dt*Source;
  Real discr=Y0*Y0+2.0*Y*temp;
  if(discr >= 0.0) Y=sgn(Y0+temp)*sqrt(discr);
  else {
    if(hybrid) Y=Y0+temp;
    else return false;
  }
  return true;
}

template<class T>
inline bool C_RK2<T>::Correct(const Complex& Y0, Complex& Y,
			   const Complex& Source)
{
  if(!Correct(Y0.re,Y.re,Source.re)) return false;
  return Correct(Y0.im,Y.im,Source.im);
}

template<class T>
inline int C_RK2<T>::Corrector()
{
  parent->ConservativeSource(Src,Y,t+halfdt);
  if(dynamic) {
    for(unsigned int j=start; j < stop; j++) {
      if(!Correct(y0[j],y[j],source[j])) return 0;
      if(!Array::Active(errmask) || errmask[j])
	CalcError(y0[j],y[j],y0[j]+dt*source0[j],y[j]);
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
    parent->IndexLimits(start,stop,startN,stopN);
  }
  
  void TimestepDependence();
  
  inline bool Correct(Real y0, Real& y, Real source0, Real source1,
		      Real source2, Real source);
  inline bool Correct(const Complex& y0, Complex& y, const Complex& source0,
		      const Complex& source1, const Complex& source2, 
		      const Complex& source);
  
  void CSource(const vector2& Src, const vector2& Y, double t) {
    parent->NonConservativeSource(Src,Y,t);
  }
  
  inline int Corrector();
  
  virtual int Corrector(unsigned int, unsigned int) {
    int rc=C_RK4<T>::Corrector();
    if(rc) rc=PC::Corrector(startN,stopN);
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
inline bool C_RK4<T>::Correct(Real Y0, Real& Y, Real Source0, Real Source1,
			      Real Source2, Real Source)
{			
  Real temp=sixthdt*(Source0+2.0*(Source1+Source2)+Source);
  Real discr=Y0*Y0+2.0*(Y0*temp+sixthdt2*(Source1*(Source0+Source2)+
					  Source2*Source));
  if(discr >= 0.0) {
    Y=sgn(Y0+temp)*sqrt(discr);
  } else {
    if(hybrid) Y=Y0+temp;
    else return false;
  }
  return true;
}


template<class T>
inline bool C_RK4<T>::Correct(const Complex& Y0, Complex& Y, 
			      const Complex& Source0, const Complex& Source1,
			      const Complex& Source2, const Complex& Source)
{			
  if(!Correct(Y0.re,Y.re,Source0.re,Source1.re,Source2.re,Source.re))
    return false;
  return Correct(Y0.im,Y.im,Source0.im,Source1.im,Source2.im,Source.im);
}

template<class T>
inline int C_RK4<T>::Corrector()
{
  parent->ConservativeSource(Src,Y,t+dt);
  if(dynamic) {
    for(unsigned int j=start; j < stop; j++) {
      Var pred=y[j]; // Fix me
      if(!Correct(y0[j],y[j],source0[j],source1[j],source2[j],source[j]))
	return 0;
      if(!Array::Active(errmask) || errmask[j]) CalcError(y0[j],y[j],pred,y[j]);
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
    parent->IndexLimits(start,stop,startN,stopN);
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
  
  void CSource(const vector2& Src, const vector2& Y, double t) {
    parent->NonConservativeSource(Src,Y,t);
  }
  
  void Predictor(unsigned int start, unsigned int stop);
  inline int Corrector();
  
  virtual int Corrector(unsigned int, unsigned int) {
    int rc=C_RK5<T>::Corrector();
    if(rc) rc=PC::Corrector(startN,stopN);
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
  for(unsigned int j=start; j < stop; j++)
    y[j]=y0[j]+b10*source0[j];
  Source(Src,Y,t+a1);
  for(unsigned int j=start; j < stop; j++)
    y2[j]=y0[j]+b20*source0[j]+b21*source[j];
  Source(Src2,Y2,t+a2);
  for(unsigned int j=start; j < stop; j++)
    y3[j]=y0[j]+b30*source0[j]+b31*source[j]+b32*source2[j];
  Source(Src3,Y3,t+a3);
  for(unsigned int j=start; j < stop; j++)
    y4[j]=y0[j]+b40*source0[j]+b41*source[j]+b42*source2[j]+b43*source3[j];
  Source(Src4,Y4,t+a4);
  for(unsigned int j=start; j < stop; j++) 
    y[j]=y0[j]+b50*source0[j]+b51*source[j]+b52*source2[j]+b53*source3[j]+
      b54*source4[j];
}

template<class T>
inline bool C_RK5<T>::Correct(Real Y0, Real& Y2, Real& Y3,
			      Real Y4, Real& Y,
			      Real Source0, Real Source2, 
			      Real Source3, Real Source4,
			      Real Source, Real& pred, Real& discr)
{
  pred=Y0*Y0+2.0*(d0*Y0*Source0+d2*Y2*Source2+d3*Y3*Source3+d4*Y4*Source4+
		  d5*Y*Source);
  discr=Y0*Y0+2.0*(c0*Y0*Source0+c2*Y2*Source2+c3*Y3*Source3+c5*Y*Source);
	
  if(discr >= 0.0)
    Y=sgn(Y0+c0*Source0+c2*Source2+c3*Source3+c5*Source)*sqrt(discr);
  else {
    if(hybrid) {
      Y=Y0+c0*Source0+c2*Source2+c3*Source3+c5*Source;
      discr=Y*Y;
    } else return false;
  }
  return true;
}

template<class T>
inline bool C_RK5<T>::Correct(const Complex& Y0, Complex& Y2, Complex& Y3,
			      const Complex& Y4, Complex& Y,
			      const Complex& Source0, const Complex& Source2, 
			      const Complex& Source3, const Complex& Source4,
			      const Complex& Source,
			      Complex& pred, Complex& corr)
{
  if(!Correct(Y0.re,Y2.re,Y3.re,Y4.re,Y.re,Source0.re,Source2.re,
	      Source3.re,Source4.re,Source.re,pred.re,corr.re)) return false;
  return Correct(Y0.im,Y2.im,Y3.im,Y4.im,Y.im,Source0.im,Source2.im,
		 Source3.im,Source4.im,Source.im,pred.im,corr.im);
}

template<class T>
inline int C_RK5<T>::Corrector()
{
  parent->ConservativeSource(Src,Y,t+a5);
  if(dynamic) {
    for(unsigned int j=start; j < stop; j++) {
      Var pred,corr;
      if(!Correct(y0[j],y2[j],y3[j],y4[j],y[j],source0[j],source2[j],
		 source3[j],source4[j],source[j],pred,corr)) return 0;
      if(!Array::Active(errmask) || errmask[j])
	CalcError(y0[j]*y0[j],corr,pred,corr);
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

#endif
