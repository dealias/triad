#ifndef __Grid2_h__
#define __Grid2_h__ 1

#include "MultiGrid.h"

namespace Array {
  
template<class T>
class Grid2 : public Grid<Array2<T>,T> {
 protected:	
  // number of points in previous and current levels and incl. boundaries
  int nx, nx1bc, nxbc, rx, offx, ox, i1, i1p, i2, i2p;
  int ny, ny1bc, nybc, ry, offy, oy, j1, j1p, j2, j2p;
  Array1<Real>::opt x,y;
  Real hx, hxinv, hx2, hx2inv;
  Real hy, hyinv, hy2, hy2inv;
  Real hxyinv;
 public:
  Grid2() {radix=2; dimension=2;}
  virtual ~Grid2() {};
	
  virtual Limits XMeshRange()=0;
  virtual Limits YMeshRange()=0;
  Real X(int i) {return x[i];}
  Real Y(int i) {return y[i];}
  Real *X() {return x;}
  Real *Y() {return y;}
  int Nx() {return nx;}
  int Ny() {return ny;}
  int Nxbc() {return nxbc;}
  int Nybc() {return nybc;}
  int Nx1bc() {return nx1bc;}
  int Ny1bc() {return ny1bc;}
  int I1() {return i1;}
  int I2() {return i2;}
  int J1() {return j1;}
  int J2() {return j2;}
  int Ox() {return ox;}
  int Oy() {return oy;}
  Real Hx() {return hx;}
  Real Hy() {return hy;}

  void Allocate(int allocate=1) {
    Mesh(x,XMeshRange(),nx,nx1bc,nxbc,hx,hxinv,hx2,hx2inv,rx,
	 offx,ox,i1,i1p,i2,i2p);
    Mesh(y,YMeshRange(),ny,ny1bc,nybc,hy,hyinv,hy2,hy2inv,ry,
	 offy,oy,j1,j1p,j2,j2p);
    hxyinv=1.0/(hx*hy);
    if(!allocate) return;
    d.Allocate(nxbc,nybc,ox,oy);
    if(level > 0) {
      v.Allocate(nx1bc,ny1bc,ox,oy);
      if(nonlinear) v2.Allocate(nx1bc,ny1bc,ox,oy);
    }
  }

  virtual void Defect(const Array2<T>& d0, const Array2<T>& u,
		      const Array2<T>& f)=0;
  virtual void Smooth(const Array2<T>& u, const Array2<T>& f)=0;
	
  virtual void GaussSeidel(const Array2<T>&, const Array2<T>&, int, int,
			   int, int) {}; 
  virtual void XGaussSeidel(const Array2<T>&, const Array2<T>&, int, int) {};
  virtual void YGaussSeidel(const Array2<T>&, const Array2<T>&, int, int)	{};

  virtual void Restrict(const Array2<T>& r, const Array2<T>& u) {
    if(&r != &u) {
      XDirichlet(r,u,1);
      YDirichlet(r,u,1);
    }
    for(int i=i1; i <= i2p; i++) {
      int i0=rx*i+offx;
      Array1<T>::opt ri=r[i];
      Array1<T>::opt um=u[i0-1]+offy, uz=u[i0]+offy, up=u[i0+1]+offy;
      for(int j=j1; j <= j2p; j++) {
	ri[j]=0.25*(0.5*(0.5*(um[ry*j-1]+um[ry*j+1]+up[ry*j-1]+
			      up[ry*j+1])+
			 um[ry*j]+uz[ry*j-1]+uz[ry*j+1]+up[ry*j])+
		    uz[ry*j]);
      }
    }
  }
	
  virtual void SubtractProlongation(const Array2<T>& u,
				    const Array2<T>& v0) {
    for(int i=i1p; i <= i2p; i++) {
      int i0=rx*i+offx;
      Array1<T>::opt vz=v0[i], vp=v0[i+1], uz=u[i0]+offy, up=u[i0+1]+offy;
      for(int j=j1p; j <= j2p; j++) {
	uz[ry*j] -= vz[j];
	up[ry*j] -= 0.5*(vz[j]+vp[j]);
	uz[ry*j+1] -= 0.5*(vz[j]+vz[j+1]);
	up[ry*j+1] -= 0.25*(vz[j]+vp[j]+vz[j+1]+vp[j+1]);
      }
    }
  }
	
  virtual inline void L0inv(const Array2<T>& u, const Array2<T>& f) {};
	
  void Jacobi(const Array2<T>& u, const Array2<T>& f, Real omegah2) {
    Defect(d,u,f);
    for(int i=i1; i <= i2; i++) {
      Array1<T>::opt di=d[i], ui=u[i];
      for(int j=j1; j <= j2; j++) {
	ui[j] -= omegah2*di[j];
      }
    }
  }
	
  void Lexicographical(const Array2<T>& u, const Array2<T>& f) {
    GaussSeidel(u,f,0,0,1,1);
  }
	
  void XLexicographical(const Array2<T>& u, const Array2<T>& f) {
    XGaussSeidel(u,f,0,1);
  }
	
  void YLexicographical(const Array2<T>& u, const Array2<T>& f) {
    YGaussSeidel(u,f,0,1);
  }
	
  void RedBlack(const Array2<T>& u, const Array2<T>& f) {
    GaussSeidel(u,f,0,0,2,2);
    GaussSeidel(u,f,1,1,2,2);
		
    GaussSeidel(u,f,0,1,2,2);
    GaussSeidel(u,f,1,0,2,2);
  }
	
  void XZebra(const Array2<T>& u, const Array2<T>& f) {
    XGaussSeidel(u,f,0,2);
    XGaussSeidel(u,f,1,2);
  }
	
  void YZebra(const Array2<T>& u, const Array2<T>& f) {
    YGaussSeidel(u,f,0,2);
    YGaussSeidel(u,f,1,2);
  }
	
  void Sum2(const Array2<T>& u, T& s) {
    for(int i=i1; i <= i2; i++) {
      Array1<T>::opt ui=u[i];
      for(int j=j1; j <= j2; j++) s += abs2(ui[j]);
    }
  }
	
  virtual inline void BoundaryConditions(const Array2<T>& u)=0;
	
  void XDirichlet(const Array2<T>&) {}
	
  void XDirichlet2(const Array2<T>&) {}
		
  void XDirichlet(const Array2<T>& u, T b0, T b1) {
    if(homogeneous) return;
    Array1<T>::opt u0=u[i1-1], unx1=u[i2+1];
    int nybco=nybc+oy;
    for(int j=oy; j < nybco; j++) {
      u0[j]=b0;
      unx1[j]=b1;
    }
  }
	
  void XDirichlet2(const Array2<T>& u, T b0, T b1) {
    if(homogeneous) return;
    Array1<T>::opt um1=u[i1-2], u0=u[i1-1], unx1=u[i2+1], unx2=u[i2+2];
    int nybco=nybc+oy;
    for(int j=oy; j < nybco; j++) {
      um1[j]=b0;
      u0[j]=b0;
      unx1[j]=b1;
      unx2[j]=b1;
    }
  }
	
  void XDirichlet(const Array2<T>& u, const Array2<T>& b, int contract=0) {
    if(homogeneous) return;
    int I2,nybco;
    if(contract) {I2=i2p; nybco=ny1bc+oy;} 
    else {I2=i2; nybco=nybc+oy;} 
    Array1<T>::opt u0=u[i1-1], unx1=u[I2+1];
    Array1<T>::opt b0=b[i1-1], bnx1=b[i2+1];
    for(int j=oy; j < nybco; j++) {
      u0[j]=b0[j];
      unx1[j]=bnx1[j];
    }
  }
	
  void XDirichlet2(const Array2<T>& u, const Array2<T>& b, int contract=0) {
    if(homogeneous) return;
    int I2,nybco;
    if(contract) {I2=i2p; nybco=ny1bc+oy;} 
    else {I2=i2; nybco=nybc+oy;} 
    Array1<T>::opt um1=u[i1-2], u0=u[i1-1], unx1=u[I2+1], unx2=u[I2+2];
    Array1<T>::opt bm1=b[i1-2], b0=b[i1-1], bnx1=b[i2+1], bnx2=b[i2+2];
    for(int j=oy; j < nybco; j++) {
      um1[j]=bm1[j];
      u0[j]=b0[j];
      unx1[j]=bnx1[j];
      unx2[j]=bnx2[j];
    }
  }
	
  void XNeumann(const Array2<T>& u) {
    Array1<T>::opt u0=u[i1-1], u2=u[i1+1];
    Array1<T>::opt unxm1=u[i2-1], unx1=u[i2+1];
    int nybco=nybc+oy;
    for(int j=oy; j < nybco; j++) {
      u0[j]=u2[j];
      unx1[j]=unxm1[j];
    }
  }
	
  void XNeumann2(const Array2<T>& u) {
    Array1<T>::opt um1=u[i1-2], u0=u[i1-1], u2=u[i1+1], u3=u[i1+2];
    Array1<T>::opt unxm2=u[i2-2], unxm1=u[i2-1], unx1=u[i2+1], unx2=u[i2+2];
    int nybco=nybc+oy;
    for(int j=oy; j < nybco; j++) {
      um1[j]=u3[j];
      u0[j]=u2[j];
      unx1[j]=unxm1[j];
      unx2[j]=unxm2[j];
    }
  }
	
  void XConstant(const Array2<T>& u) {
    Array1<T>::opt u0=u[i1-1], u1=u[i1];
    Array1<T>::opt unx=u[i2], unx1=u[i2+1];
    int nybco=nybc+oy;
    for(int j=oy; j < nybco; j++) {
      u0[j]=u1[j];
      unx1[j]=unx[j];
    }
  }
	
  void XConstant2(const Array2<T>& u) {
    Array1<T>::opt um1=u[i1-2], u0=u[i1-1], u1=u[i1];
    Array1<T>::opt unx=u[i2], unx1=u[i2+1], unx2=u[i2+2];
    int nybco=nybc+oy;
    for(int j=oy; j < nybco; j++) {
      u0[j]=um1[j]=u1[j];
      unx2[j]=unx1[j]=unx[j];
    }
  }
	
  void XMixedA(const Array2<T>& u) {
    Array1<T>::opt u0=u[i1-1], u2=u[i1+1];
    int nybco=nybc+oy;
    for(int j=oy; j < nybco; j++) {
      u0[j]=u2[j];
    }
  }
	
  void XMixedA2(const Array2<T>& u) {
    Array1<T>::opt um1=u[i1-2], u0=u[i1-1], u2=u[i1+1], u3=u[i1+2];
    int nybco=nybc+oy;
    for(int j=oy; j < nybco; j++) {
      um1[j]=u3[j];
      u0[j]=u2[j];
    }
  }
	
  void XMixedB(const Array2<T>& u) {
    Array1<T>::opt unxm1=u[i2-1], unx1=u[i2+1];
    int nybco=nybc+oy;
    for(int j=oy; j < nybco; j++) {
      unx1[j]=unxm1[j];
    }
  }
	
  void XMixedB2(const Array2<T>& u) {
    Array1<T>::opt unxm2=u[i2-2], unxm1=u[i2-1], unx1=u[i2+1], unx2=u[i2+2];
    int nybco=nybc+oy;
    for(int j=oy; j < nybco; j++) {
      unx1[j]=unxm1[j];
      unx2[j]=unxm2[j];
    }
  }
	
  void XPeriodic(const Array2<T>& u) {
    Array1<T>::opt u0=u[i1-1], u1=u[i1];
    Array1<T>::opt unx=u[i2], unx1=u[i2+1];
    int nybco=nybc+oy;
    for(int j=oy; j < nybco; j++) {
      u0[j]=unx[j];
      unx1[j]=u1[j];
    }
  }
	
  void XPeriodic2(const Array2<T>& u) {
    Array1<T>::opt um1=u[i1-2], u0=u[i1-1], u1=u[i1], u2=u[i1+1];
    Array1<T>::opt unxm1=u[i2-1], unx=u[i2], unx1=u[i2+1], unx2=u[i2+2];
    int nybco=nybc+oy;
    for(int j=oy; j < nybco; j++) {
      um1[j]=unxm1[j];
      u0[j]=unx[j];
      unx1[j]=u1[j];
      unx2[j]=u2[j];
    }
  }
	
  void YDirichlet(const Array2<T>&) {}
	
  void YDirichlet2(const Array2<T>&) {}
	
  void YDirichlet(const Array2<T>& u, T b0, T b1) {
    if(homogeneous) return;
    int nxbco=nxbc+ox;
    for(int i=ox; i < nxbco; i++) {
      Array1<T>::opt ui=u[i];
      ui[j1-1]=b0;
      ui[j2+1]=b1;
    }
  }
	
  void YDirichlet2(const Array2<T>& u, T b0, T b1) {
    if(homogeneous) return;
    int nxbco=nxbc+ox;
    for(int i=ox; i < nxbco; i++) {
      Array1<T>::opt ui=u[i];
      ui[j1-2]=b0;
      ui[j1-1]=b0;
      ui[j2+1]=b1;
      ui[j2+2]=b1;
    }
  }
	
  void YDirichlet(const Array2<T>& u, const Array2<T>& b, int contract=0) {
    if(homogeneous) return;
    int nxbco,J2;
    if(contract) {nxbco=nx1bc+ox; J2=j2p;} 
    else {nxbco=nxbc+ox; J2=j2;} 
    for(int i=ox; i < nxbco; i++) {
      Array1<T>::opt ui=u[i];
      Array1<T>::opt bi=b[i];
      ui[j1-1]=bi[j1-1];
      ui[J2+1]=bi[j2+1];
    }
  }
	
  void YDirichlet2(const Array2<T>& u, const Array2<T>& b, int contract=0) {
    if(homogeneous) return;
    int nxbco,J2;
    if(contract) {nxbco=nx1bc+ox; J2=j2p;} 
    else {nxbco=nxbc+ox; J2=j2;} 
    for(int i=ox; i < nxbco; i++) {
      Array1<T>::opt ui=u[i];
      Array1<T>::opt bi=b[i];
      ui[j1-2]=bi[j1-2];
      ui[j1-1]=bi[j1-1];
      ui[J2+1]=bi[j2+1];
      ui[J2+2]=bi[j2+2];
    }
  }
	
  void YNeumann(const Array2<T>& u) {
    int nxbco=nxbc+ox;
    for(int i=ox; i < nxbco; i++) {
      Array1<T>::opt ui=u[i];
      ui[j1-1]=ui[j1+1];
      ui[j2+1]=ui[j2-1];
    }
  }
	
  void YNeumann2(const Array2<T>& u) {
    int nxbco=nxbc+ox;
    for(int i=ox; i < nxbco; i++) {
      Array1<T>::opt ui=u[i];
      ui[j1-2]=ui[j1+2];
      ui[j1-1]=ui[j1+1];
      ui[j2+1]=ui[j2-1];
      ui[j2+2]=ui[j2-2];
    }
  }
	
  void YConstant(const Array2<T>& u) {
    int nxbco=nxbc+ox;
    for(int i=ox; i < nxbco; i++) {
      Array1<T>::opt ui=u[i];
      ui[j1-1]=ui[j1];
      ui[j2+1]=ui[j2];
    }
  }
	
  void YConstant2(const Array2<T>& u) {
    int nxbco=nxbc+ox;
    for(int i=ox; i < nxbco; i++) {
      Array1<T>::opt ui=u[i];
      ui[j1-1]=ui[j1-2]=ui[j1];
      ui[j2+2]=ui[j2+1]=ui[j2];
    }
  }
	
  void YPeriodic(const Array2<T>& u) {
    int nxbco=nxbc+ox;
    for(int i=ox; i < nxbco; i++) {
      Array1<T>::opt ui=u[i];
      ui[j1-1]=ui[j2];
      ui[j2+1]=ui[j1];
    }
  }
	
  void YPeriodic2(const Array2<T>& u) {
    int nxbco=nxbc+ox;
    for(int i=ox; i < nxbco; i++) {
      Array1<T>::opt ui=u[i];
      ui[j1-2]=ui[j2-1];
      ui[j1-1]=ui[j2];
      ui[j2+1]=ui[j1];
      ui[j2+2]=ui[j1+1];
    }
  }
};
  
}

#endif
