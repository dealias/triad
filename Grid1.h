#ifndef __Grid1_h__
#define __Grid1_h__ 1

#include "MultiGrid.h"

namespace Array {
  
template<class T>
class Grid1 : public Grid<Array1<T>,T> {
 protected:	
  // number of points in previous and current levels and incl. boundaries
  // index limits
  int nx, nx1bc, nxbc, rx, offx, ox, i1, i1p, i2, i2p;
  typename Array1<Real>::opt x;
  Real hx, hxinv, hx2, hx2inv;
 public:
  Grid1() {radix=2; dimension=1;}
  virtual ~Grid1() {};

  virtual Limits XMeshRange()=0;
  Real X(int i) {return x[i];}
  Real *X() {return x;}
  int Nx() {return nx;}
  int Nxbc() {return nxbc;}
  int Nx1bc() {return nx1bc;}
  int I1() {return i1;}
  int I2() {return i2;}
  int Ox() {return ox;}
  Real Hx() {return hx;}
	
  void Allocate(int allocate=1) {
    Mesh(x,XMeshRange(),nx,nx1bc,nxbc,hx,hxinv,hx2,hx2inv,rx,
	 offx,ox,i1,i1p,i2,i2p);
    if(!allocate) return;
    d.Allocate(nxbc,ox);
    if(level > 0) {
      v.Allocate(nx1bc,ox);
      if(nonlinear) v2.Allocate(nx1bc,ox);
    }
  }

  virtual void Defect(const Array1<T>& d0, const Array1<T>& u, const
		      Array1<T>& f)=0; 
  virtual void Smooth(const Array1<T>& u, const Array1<T>& f)=0;
	
  virtual void GaussSeidel(const Array1<T>&, const Array1<T>&,
			   int, int) {};
	
  virtual void Restrict(const Array1<T>& r, const Array1<T>& u) {
    if(&r != &u) XDirichlet(r,u,1);
    typename Array1<T>::opt u0=u+offx;
    for(int i=i1; i <= i2p; i++) {
      r[i]=0.5*(0.5*(u0[rx*i-1]+u0[rx*i+1])+u0[rx*i]);
    }
  }
	
  virtual void SubtractProlongation(const Array1<T>& u,
				    const Array1<T>& v0) { 
    typename Array1<T>::opt u0=u+offx;
    for(int i=i1p; i <= i2p; i++) {
      u0[rx*i] -= v0[i];
      u0[rx*i+1] -= 0.5*(v0[i]+v0[i+1]);
    }
  }
	
  virtual inline void L0inv(const Array1<T>& u,
			    const Array1<T>& f) {};
	
  void Jacobi(const Array1<T>& u, const Array1<T>& f, Real omegah2) {
    Defect(d,u,f);
    for(int i=i1; i <= i2; i++) u[i] -= omegah2*d[i];
  }
	
  void Lexicographical(const Array1<T>& u, const Array1<T>& f) {
    GaussSeidel(u,f,0,1);
  }
	
  void RedBlack(const Array1<T>& u, const Array1<T>& f) {
    GaussSeidel(u,f,0,2);
    GaussSeidel(u,f,1,2);
  }
	
  void Sum2(const Array1<T>& u, T& s) {
    for(int i=i1; i <= i2; i++) s += abs2(u[i]);
  }

  virtual inline void BoundaryConditions(const Array1<T>& u)=0;
	
  void XDirichlet(const Array1<T>&) {}
	
  void XDirichlet2(const Array1<T>&) {}
	
  void XDirichlet(const Array1<T>& u, T b0, T b1) {
    if(homogeneous) return;
    u[i1-1]=b0;
    u[i2+1]=b1;
  }
	
  void XDirichlet2(const Array1<T>& u, T b0, T b1) {
    if(homogeneous) return;
    u[i1-2]=2.0*b0-u[i1];
    u[i1-1]=b0;
    u[i2+1]=b1;
    u[i2+2]=2.0*b1-u[i2];
  }
	
  void XDirichlet(const Array1<T>& u, const Array1<T>& b, int contract=0) {
    if(homogeneous) return;
    int I2;
    if(contract) I2=i2p;
    else I2=i2;
    u[i1-1]=b[i1-1];
    u[I2+1]=b[i2+1];
  }
	
  void XDirichlet2(const Array1<T>& u, const Array1<T>& b,
		   int contract=0) {
    if(homogeneous) return;
    int I2;
    if(contract) I2=i2p;
    else I2=i2;
    u[i1-2]=b[i1-2];
    u[i1-1]=b[i1-1];
    u[I2+1]=b[i2+1];
    u[I2+2]=b[i2+2];
  }
	
  void XNeumann(const Array1<T>& u) {
    u[i1-1]=u[i1+1];
    u[i2+1]=u[i2-1];
  }
	
  void XNeumann2(const Array1<T>& u) {
    u[i1-2]=u[i1+2];
    u[i1-1]=u[i1+1];
    u[i2+1]=u[i2-1];
    u[i2+2]=u[i2-2];
  }
	
  void XConstant(const Array1<T>& u) {
    u[i1-1]=u[i1];
    u[i2+1]=u[i2];
  }
	
  void XConstant2(const Array1<T>& u) {
    u[i1-1]=u[i1-2]=u[i1];
    u[i2+2]=u[i2+1]=u[i2];
  }
	
  void XMixedA(const Array1<T>& u) {
    u[i1-1]=u[i1+1];
  }
	
  void XMixedA2(const Array1<T>& u) {
    u[i1-2]=u[i1+2];
    u[i1-1]=u[i1+1];
  }
	
  void XMixedB(const Array1<T>& u) {
    u[i2+1]=u[i2-1];
  }
	
  void XMixedB2(const Array1<T>& u) {
    u[i2+1]=u[i2-1];
    u[i2+2]=u[i2-2];
  }
	
  void XPeriodic(const Array1<T>& u) {
    u[i1-1]=u[i2];
    u[i2+1]=u[i1];
  }
	
  void XPeriodic2(const Array1<T>& u) {
    u[i1-2]=u[i2-1];
    u[i1-1]=u[i2];
    u[i2+1]=u[i1];
    u[i2+2]=u[i1+1];
  }
};

}

#endif
