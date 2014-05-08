#ifndef __Poisson_h__
#define __Poisson_h__ 1

namespace Array {

  template<class T>
  class Poisson2 : public Grid2<T> {
  public:	
    Poisson2() {}
    void Defect(const Array2<T>& d0, const Array2<T>& u, const Array2<T>& f);
    void GaussSeidel(const Array2<T>&, const Array2<T>&, int, int, int, int);
    void Smooth(const Array2<T>& u, const Array2<T>& f);
    virtual inline void BoundaryConditions(const Array2<T>& u)=0;
  };

  template<class T>
  void Poisson2<T>::Defect(const Array2<T>& d0, const Array2<T>& u,
			   const Array2<T>& f)
  {
    T coeffx2y2=-2.0*(this->hx2inv+this->hy2inv);
    for(int i=this->i1; i <= this->i2; i++) {
      typename Array1<T>::opt di=d0[i];
      typename Array1<T>::opt fi=f[i];
      typename Array1<T>::opt um=u[i-1];
      typename Array1<T>::opt ui=u[i];
      typename Array1<T>::opt up=u[i+1];
      for(int j=this->j1; j <= this->j2; j++) {
	di[j]=this->hx2inv*(um[j]+up[j])+coeffx2y2*ui[j]+
	  this->hy2inv*(ui[j-1]+ui[j+1])-fi[j];
      }
    }
  }

  template<class T>
  void Poisson2<T>::GaussSeidel(const Array2<T>& u, const Array2<T>& f,
				int i0, int j0, int istep, int jstep)
  {
    T factor=0.5/(this->hx2inv+this->hy2inv);
    for(int i=this->i1+i0; i <= this->i2; i += istep) {
      typename Array1<T>::opt fi=f[i];
      typename Array1<T>::opt um=u[i-1];
      typename Array1<T>::opt ui=u[i];
      typename Array1<T>::opt up=u[i+1];
      for(int j=this->j1+j0; j <= this->j2; j += jstep) {
	ui[j]=factor*(this->hx2inv*(um[j]+up[j])+
		      this->hy2inv*(ui[j-1]+ui[j+1])-fi[j]);
      }
    }
  }

  template<class T>
  void Poisson2<T>::Smooth(const Array2<T>& u, const Array2<T>& f)
  {
    this->RedBlack(u,f);
  }

}

#endif
