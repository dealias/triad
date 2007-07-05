#ifndef __Helmholtz_h__
#define __Helmholtz_h__ 1

namespace Array {
  
// Solve the Helmholtz equation
//
// A\Del^2 u+Bu=f

  template<class T>
  class Helmholtz2 : public Grid2<T> {
  public:	
    virtual T A()=0;
    virtual T B()=0;
    void Defect(const Array2<T>& d0, const Array2<T>& u, const Array2<T>& f);
    void GaussSeidel(const Array2<T>&, const Array2<T>&, int, int, int, int);
    void Smooth(const Array2<T>& u, const Array2<T>& f);
    virtual void BoundaryConditions(const Array2<T>& u)=0;
  };

  template<class T>
  void Helmholtz2<T>::Defect(const Array2<T>& d0, const Array2<T>& u,
			     const Array2<T>& f)
  {
    T coeffx2=this->hx2inv*A();
    T coeffy2=this->hy2inv*A();
    T coeffx2y2=-2.0*(coeffx2+coeffy2)+B();
    for(int i=this->i1; i <= this->i2; i++) {
      typename Array1<T>::opt di=d0[i];
      typename Array1<T>::opt fi=f[i];
      typename Array1<T>::opt um=u[i-1];
      typename Array1<T>::opt ui=u[i];
      typename Array1<T>::opt up=u[i+1];
      for(int j=this->j1; j <= this->j2; j++) {
	di[j]=coeffx2*(um[j]+up[j])+coeffx2y2*ui[j]+
	  coeffy2*(ui[j-1]+ui[j+1])-fi[j];
      }
    }
  }

  template<class T>
  void Helmholtz2<T>::GaussSeidel(const Array2<T>& u, const Array2<T>& f,
				  int i0, int j0, int istep, int jstep)
  {
    T coeffx2=this->hx2inv*A();
    T coeffy2=this->hy2inv*A();
    T factor=1.0/(2.0*(coeffx2+coeffy2)-B());
    for(int i=this->i1+i0; i <= this->i2; i += istep) {
      typename Array1<T>::opt fi=f[i];
      typename Array1<T>::opt um=u[i-1];
      typename Array1<T>::opt ui=u[i];
      typename Array1<T>::opt up=u[i+1];
      for(int j=this->j1+j0; j <= this->j2; j += jstep) {
	ui[j]=factor*(coeffx2*(um[j]+up[j])+coeffy2*(ui[j-1]+ui[j+1])-
		      fi[j]);
      }
    }
  }

  template<class T>
  void Helmholtz2<T>::Smooth(const Array2<T>& u, const Array2<T>& f)
  {
    RedBlack(u,f);
  }

}

#endif
