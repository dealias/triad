#ifndef __Grid3_h__
#define __Grid3_h__ 1

#include "MultiGrid.h"

namespace Array {
  
template<class T>
class Grid3 : public Grid<Array3<T>,T> {
 protected:	
  // number of points in previous and current levels and incl. boundaries
  int nx, nx1bc, nxbc, rx, offx, ox, i1, i1p, i2, i2p;
  int ny, ny1bc, nybc, ry, offy, oy, j1, j1p, j2, j2p;
  int nz, nz1bc, nzbc, rz, offz, oz, k1, k1p, k2, k2p;
  typename Array1<Real>::opt x,y,z;
  Real hx, hxinv, hx2, hx2inv;
  Real hy, hyinv, hy2, hy2inv;
  Real hz, hzinv, hz2, hz2inv;
  Real hxyinv, hxzinv, hyzinv;
 public:
  Grid3() {radix=2; dimension=3;}
  virtual ~Grid3() {};

  virtual Limits XMeshRange()=0;
  virtual Limits YMeshRange()=0;
  virtual Limits ZMeshRange()=0;
	
  Real X(int i) {return x[i];}
  Real Y(int i) {return y[i];}
  Real Z(int i) {return z[i];}
  Real *X() {return x;}
  Real *Y() {return y;}
  Real *Z() {return z;}
  int Nx() {return nx;}
  int Ny() {return ny;}
  int Nz() {return nz;}
  int Nxbc() {return nxbc;}
  int Nybc() {return nybc;}
  int Nzbc() {return nzbc;}
  int Nx1bc() {return nx1bc;}
  int Ny1bc() {return ny1bc;}
  int Nz1bc() {return nz1bc;}
  int I1() {return i1;}
  int I2() {return i2;}
  int J1() {return j1;}
  int J2() {return j2;}
  int K1() {return k1;}
  int K2() {return k2;}
  int Ox() {return ox;}
  int Oy() {return oy;}
  int Oz() {return oz;}
  Real Hx() {return hx;}
  Real Hy() {return hy;}
  Real Hz() {return hz;}

  void Allocate(int allocate=1) {
    Mesh(x,XMeshRange(),nx,nx1bc,nxbc,hx,hxinv,hx2,hx2inv,rx,
	 offx,ox,i1,i1p,i2,i2p);
    Mesh(y,YMeshRange(),ny,ny1bc,nybc,hy,hyinv,hy2,hy2inv,ry,
	 offy,oy,j1,j1p,j2,j2p);
    Mesh(z,ZMeshRange(),nz,nz1bc,nzbc,hz,hzinv,hz2,hz2inv,rz,
	 offz,oz,k1,k1p,k2,k2p);
    hxyinv=1.0/(hx*hy); hxzinv=1.0/(hx*hz); hyzinv=1.0/(hy*hz);
    if(!allocate) return;
    d.Allocate(nxbc,nybc,nzbc,ox,oy,oz);
    if(level > 0) {
      v.Allocate(nx1bc,ny1bc,nz1bc,ox,oy,oz);
      if(nonlinear) v2.Allocate(nx1bc,ny1bc,nz1bc,ox,oy,oz);
    }
  }

  virtual void Defect(const Array3<T>& d0, const Array3<T>& u,
		      const Array3<T>& f)=0;
  virtual void Smooth(const Array3<T>& u, const Array3<T>& f)=0;
  virtual void GaussSeidel(const Array3<T>&, const Array3<T>&, int, int, int,
			   int, int, int) {};
	
  virtual void XGaussSeidel(const Array3<T>&, const Array3<T>&, int, int,
			    int, int) {}; 
  virtual void YGaussSeidel(const Array3<T>&, const Array3<T>&, int, int,
			    int, int) {} ; 
  virtual void ZGaussSeidel(const Array3<T>&, const Array3<T>&, int, int,
			    int, int) {};
  virtual void XYGaussSeidel(const Array3<T>&, const Array3<T>&,
			     int, int) {};
  virtual void XZGaussSeidel(const Array3<T>&, const Array3<T>&,
			     int, int) {};
  virtual void YZGaussSeidel(const Array3<T>&, const Array3<T>&,
			     int, int) {};

  virtual void Restrict(const Array3<T>& r, const Array3<T>& u) {
    if(&r != &u) {
      XDirichlet(r,u,1);
      YDirichlet(r,u,1);
      ZDirichlet(r,u,1);
    }
    for(int i=i1; i <= i2p; i++) {
      int i0=rx*i+offx;
      Array2<T> ri=r[i], um=u[i0-1], uz=u[i0], up=u[i0+1];
      for(int j=j1; j <= j2p; j++) {
	int j0=ry*j+offy;
	typename Array1<T>::opt rij=ri[j];
	typename Array1<T>::opt umm=um[j0-1]+offz, umz=um[j0]+offz;
	typename Array1<T>::opt ump=um[j0+1]+offz, uzm=uz[j0-1]+offz;
	typename Array1<T>::opt uzz=uz[j0]+offz, uzp=uz[j0+1]+offz;
	typename Array1<T>::opt upm=up[j0-1]+offz, upz=up[j0]+offz;
	typename Array1<T>::opt upp=up[j0+1]+offz;
	for(int k=k1; k <= k2p; k++) {
	  rij[k]=0.125*(0.5*(0.5*(0.5*(
				       umm[rz*k-1]+umm[rz*k+1]+ump[rz*k-1]+
				       ump[rz*k+1]+um[rz*j-1][rz*k-1]+
				       umm[rz*k+1]+upp[rz*k-1]+upp[rz*k+1])+
				  umm[rz*k]+umz[rz*k-1]+umz[rz*k+1]+
				  ump[rz*k]+uzm[rz*k-1]+uzm[rz*k+1]+
				  uzp[rz*k-1]+uzp[rz*k+1]+upm[rz*k]+
				  upz[rz*k-1]+upz[rz*k+1]+upp[rz*k])+
			     umz[rz*k]+uzm[rz*k]+uzz[rz*k-1]+
			     uzz[rz*k+1]+uzp[rz*k]+upz[rz*k])
			+uzz[rz*k]);
	}
      }
    }
  }
	
  virtual void SubtractProlongation(const Array3<T>& u,
				    const Array3<T>& v0) {
    for(int i=i1p; i <= i2p; i++) {
      int i0=rx*i+offx;
      Array2<T> uz=u[i0], up=u[i0+1], vz=v0[i], vp=v0[i+1];
      for(int j=j1p; j <= j2p; j++) {
	int j0=ry*j+offy;
	typename Array1<T>::opt uzz=uz[j0]+offz, uzp=uz[j0+1]+offz;
	typename Array1<T>::opt upz=up[j0]+offz, upp=up[j0+1]+offz;
	typename Array1<T>::opt vzz=vz[j], vzp=vz[j+1], vpz=vp[j], vpp=vp[j+1];
	for(int k=k1p; k <= k2p; k++) {
	  uzz[rz*k] -= vzz[k];
	  upz[rz*k] -= 0.5*(vzz[k]+vpz[k]);
	  uzp[rz*k] -= 0.5*(vzz[k]+vzp[k]);
	  upp[rz*k] -= 0.25*(vzz[k]+vpz[k]+vzp[k]+vpp[k]);
	  uzz[rz*k+1] -= 0.5*(vzz[k]+vzz[k+1]);
	  upz[rz*k+1] -= 0.25*(vzz[k]+vpz[k]+vzz[k+1]+vpz[k+1]);
	  uzp[rz*k+1] -= 0.25*(vzz[k]+vzp[k]+vzz[k+1]+vzp[k+1]);
	  upp[rz*k+1] -= 0.125*(vzz[k]+vpz[k]+vzp[k]+vpp[k]+
				vzz[k+1]+vpz[k+1]+vzp[k+1]+vpp[k+1]);
	}
      }
    }
  }
	
  virtual inline void L0inv(const Array3<T>& u, const Array3<T>& f) {};
	
  void Jacobi(const Array3<T>& u, const Array3<T>& f, Real omegah2) {
    Defect(d,u,f);
    for(int i=i1; i <= i2; i++) {
      Array2<T> di=d[i], ui=u[i];
      for(int j=j1; j <= j2; j++) {
	typename Array1<T>::opt dij=di[j], uij=ui[j];
	for(int k=k1; k <= k2; k++) {
	  uij[k] -= omegah2*dij[k];
	}
      }
    }
  }
	
  void Lexicographical(const Array3<T>& u, const Array3<T>& f) {
    GaussSeidel(u,f,0,0,0,1,1,1);
  }
	
  void XLexicographical(const Array3<T>& u, const Array3<T>& f) {
    XGaussSeidel(u,f,0,0,1,1);
  }
	
  void YLexicographical(const Array3<T>& u, const Array3<T>& f) {
    YGaussSeidel(u,f,0,0,1,1);
  }
	
  void ZLexicographical(const Array3<T>& u, const Array3<T>& f) {
    ZGaussSeidel(u,f,0,0,1,1);
  }
	
  void XYLexicographical(const Array3<T>& u, const Array3<T>& f) {
    XYGaussSeidel(u,f,0,1);
  }
	
  void XZLexicographical(const Array3<T>& u, const Array3<T>& f) {
    XZGaussSeidel(u,f,0,1);
  }
	
  void YZLexicographical(const Array3<T>& u, const Array3<T>& f) {
    YZGaussSeidel(u,f,0,1);
  }
	
  void RedBlack(const Array3<T>& u, const Array3<T>& f) {
    GaussSeidel(u,f,0,0,0,2,2,2);
    GaussSeidel(u,f,0,1,1,2,2,2);
    GaussSeidel(u,f,1,0,1,2,2,2);
    GaussSeidel(u,f,1,1,0,2,2,2);
		
    GaussSeidel(u,f,0,0,1,2,2,2);
    GaussSeidel(u,f,0,1,0,2,2,2);
    GaussSeidel(u,f,1,0,0,2,2,2);
    GaussSeidel(u,f,1,1,1,2,2,2);
  }
	
  void XZebra(const Array3<T>& u, const Array3<T>& f) {
    XGaussSeidel(u,f,0,0,2,2);
    XGaussSeidel(u,f,1,1,2,2);
    XGaussSeidel(u,f,0,1,2,2);
    XGaussSeidel(u,f,1,0,2,2);
  }
	
  void YZebra(const Array3<T>& u, const Array3<T>& f) {
    YGaussSeidel(u,f,0,0,2,2);
    YGaussSeidel(u,f,1,1,2,2);
    YGaussSeidel(u,f,0,1,2,2);
    YGaussSeidel(u,f,1,0,2,2);
  }
	
  void ZZebra(const Array3<T>& u, const Array3<T>& f) {
    ZGaussSeidel(u,f,0,0,2,2);
    ZGaussSeidel(u,f,1,1,2,2);
    ZGaussSeidel(u,f,0,1,2,2);
    ZGaussSeidel(u,f,1,0,2,2);
  }
	
  void XYZebra(const Array3<T>& u, const Array3<T>& f) {
    XYGaussSeidel(u,f,0,2);
    XYGaussSeidel(u,f,1,2);
  }
	
  void XZZebra(const Array3<T>& u, const Array3<T>& f) {
    XZGaussSeidel(u,f,0,2);
    XZGaussSeidel(u,f,1,2);
  }
	
  void YZZebra(const Array3<T>& u, const Array3<T>& f) {
    YZGaussSeidel(u,f,0,2);
    YZGaussSeidel(u,f,1,2);
  }
	
  void Sum2(const Array3<T>& u, T& s) {
    for(int i=i1; i <= i2; i++) {
      Array2<T> ui=u[i];
      for(int j=j1; j <= j2; j++) {
	typename Array1<T>::opt uij=ui[j];
	for(int k=k1; k <= k2; k++) s += abs2(uij[k]);
      }
    }
  }
	
  virtual inline void BoundaryConditions(const Array3<T>& u)=0;
	
  void XDirichlet(const Array3<T>&) {}
	
  void XDirichlet(const Array3<T>& u, T b0, T b1) {
    if(homogeneous) return;
    Array2<T> u0=u[i1-1], unx1=u[i2+1];
    for(int j=0; j < nybc; j++) {
      typename Array1<T>::opt u0j=u0[j], unx1j=unx1[j];
      for(int k=0; k < nzbc; k++) {
	u0j[k]=b0;
	unx1j[k]=b1;
      }
    }
  }
	
  void XDirichlet(const Array3<T>& u, const Array3<T>& b, int contract=0) {
    if(homogeneous) return;
    int I2,ny0bc,nz0bc;
    if(contract) {I2=i2p; ny0bc=ny1bc; nz0bc=nz1bc;}
    else {I2=i2; ny0bc=nybc; nz0bc=nzbc;} 
    Array2<T> u0=u[i1-1], unx1=u[I2+1];
    Array2<T> b0=b[i1-1], bnx1=b[i2+1];
    for(int j=0; j < ny0bc; j++) {
      typename Array1<T>::opt u0j=u0[j], unx1j=unx1[j];
      typename Array1<T>::opt b0j=b0[j], bnx1j=bnx1[j];
      for(int k=0; k < nz0bc; k++) {
	u0j[k]=b0j[k];
	unx1j[k]=bnx1j[k];
      }
    }
  }
	
  void XNeumann(const Array3<T>& u) {
    Array2<T> u0=u[i1-1], u2=u[i1+1], unxm1=u[i2-1], unx1=u[i2+1];
    for(int j=0; j < nybc; j++) {
      typename Array1<T>::opt u0j=u0[j], u2j=u2[j];
      typename Array1<T>::opt unxm1j=unxm1[j], unx1j=unx1[j];
      for(int k=0; k < nzbc; k++) {
	u0j[k]=u2j[k];
	unx1j[k]=unxm1j[k];
      }
    }
  }
	
  void XConstant(const Array3<T>& u) {
    Array2<T> u0=u[i1-1], u1=u[i1], unx=u[i2], unx1=u[i2+1];
    for(int j=0; j < nybc; j++) {
      typename Array1<T>::opt u0j=u0[j], u1j=u1[j], unxj=unx[j], unx1j=unx1[j];
      for(int k=0; k < nzbc; k++) {
	u0j[k]=u1j[k];
	unx1j[k]=unxj[k];
      }
    }
  }
	
  void XMixedA(const Array3<T>& u) {
    Array2<T> u0=u[i1-1], u2=u[i1+1];
    for(int j=0; j < nybc; j++) {
      typename Array1<T>::opt u0j=u0[j], u2j=u2[j];
      for(int k=0; k < nzbc; k++) {
	u0j[k]=u2j[k];
      }
    }
  }
	
  void XMixedB(const Array3<T>& u) {
    Array2<T> unxm1=u[i2-1], unx1=u[i2+1];
    for(int j=0; j < nybc; j++) {
      typename Array1<T>::opt unxm1j=unxm1[j], unx1j=unx1[j];
      for(int k=0; k < nzbc; k++) {
	unx1j[k]=unxm1j[k];
      }
    }
  }
	
  void XPeriodic(const Array3<T>& u) {
    Array2<T> u0=u[i1-1], u1=u[i1], unx=u[i2], unx1=u[i2+1];
    for(int j=0; j < nybc; j++) {
      typename Array1<T>::opt u0j=u0[j], u1j=u1[j], unxj=unx[j], unx1j=unx1[j];
      for(int k=0; k < nzbc; k++) {
	u0j[k]=unxj[k];
	unx1j[k]=u1j[k];
      }
    }
  }
	
  void YDirichlet(const Array3<T>&) {}
	
  void YDirichlet(const Array3<T>& u, T b0, T b1) {
    if(homogeneous) return;
    for(int i=0; i < nxbc; i++) {
      Array2<T> ui=u[i];
      typename Array1<T>::opt ui0=ui[j1-1], uiny1=ui[j2+1];
      for(int k=0; k < nzbc; k++) {
	ui0[k]=b0;
	uiny1[k]=b1;
      }
    }
  }
	
  void YDirichlet(const Array3<T>& u, const Array3<T>& b, int contract=0) {
    if(homogeneous) return;
    int nx0bc,J2,nz0bc;
    if(contract) {nx0bc=nx1bc; J2=j2p; nz0bc=nz1bc;} 
    else {nx0bc=nxbc; J2=j2; nz0bc=nzbc;} 
    for(int i=0; i < nx0bc; i++) {
      Array2<T> ui=u[i], bi=b[i];
      typename Array1<T>::opt ui0=ui[j1-1], uiny1=ui[J2+1];
      typename Array1<T>::opt bi0=bi[j1-1], biny1=bi[j2+1];
      for(int k=0; k < nz0bc; k++) {
	ui0[k]=bi0[k];
	uiny1[k]=biny1[k];
      }
    }
  }
	
  void YNeumann(const Array3<T>& u) {
    for(int i=0; i < nxbc; i++) {
      Array2<T> ui=u[i];
      typename Array1<T>::opt ui0=ui[j1-1], ui2=ui[j1+1];
      typename Array1<T>::opt uinym1=ui[j2-1], uiny1=ui[j2+1];
      for(int k=0; k < nzbc; k++) {
	ui0[k]=ui2[k];
	uiny1[k]=uinym1[k];
      }
    }
  }
	
  void YConstant(const Array3<T>& u) {
    for(int i=0; i < nxbc; i++) {
      Array2<T> ui=u[i];
      typename Array1<T>::opt ui0=ui[j1-1], ui1=ui[j1], uiny=ui[j2], uiny1=ui[j2+1];
      for(int k=0; k < nzbc; k++) {
	ui0[k]=ui1[k];
	uiny1[k]=uiny[k];
      }
    }
  }
	
  void YPeriodic(const Array3<T>& u) {
    for(int i=0; i < nxbc; i++) {
      Array2<T> ui=u[i];
      typename Array1<T>::opt ui0=ui[j1-1], ui1=ui[j1], uiny=ui[j2], uiny1=ui[j2+1];
      for(int k=0; k < nzbc; k++) {
	ui0[k]=uiny[k];
	uiny1[k]=ui1[k];
      }
    }
  }
	
  void ZDirichlet(const Array3<T>&) {}
		
  void ZDirichlet(const Array3<T>& u, T b0, T b1) {
    if(homogeneous) return;
    for(int i=0; i < nxbc; i++) {
      Array2<T> ui=u[i];
      for(int j=0; j < nybc; j++) {
	typename Array1<T>::opt uij=ui[j];
	uij[k1-1]=b0;
	uij[k2+1]=b1;
      }
    }
  }
	
  void ZDirichlet(const Array3<T>& u, const Array3<T>& b, int contract=0) {
    if(homogeneous) return;
    int nx0bc,ny0bc,K2;
    if(contract) {nx0bc=nx1bc; ny0bc=ny1bc; K2=k2p;}
    else {nx0bc=nxbc; ny0bc=nybc; K2=k2;} 
    for(int i=0; i < nx0bc; i++) {
      Array2<T> ui=u[i], bi=b[i];
      for(int j=0; j < ny0bc; j++) {
	typename Array1<T>::opt uij=ui[j];
	typename Array1<T>::opt bij=bi[j];
	uij[k1-1]=bij[k1-1];
	uij[K2+1]=bij[k2+1];
      }
    }
  }
	
  void ZNeumann(const Array3<T>& u) {
    for(int i=0; i < nxbc; i++) {
      Array2<T> ui=u[i];
      for(int j=0; j < nybc; j++) {
	typename Array1<T>::opt uij=ui[j];
	uij[k1-1]=uij[k1+1];
	uij[k2+1]=uij[k2-1];
      }
    }
  }
	
  void ZConstant(const Array3<T>& u) {
    for(int i=0; i < nxbc; i++) {
      Array2<T> ui=u[i];
      for(int j=0; j < nybc; j++) {
	typename Array1<T>::opt uij=ui[j];
	uij[k1-1]=uij[k1];
	uij[k2+1]=uij[k2];
      }
    }
  }
	
  void ZPeriodic(const Array3<T>& u) {
    for(int i=0; i < nxbc; i++) {
      Array2<T> ui=u[i];
      for(int j=0; j < nybc; j++) {
	typename Array1<T>::opt uij=ui[j];
	uij[k1-1]=uij[k2];
	uij[k2+1]=uij[k1];
      }
    }
  }
};

}
  
#endif
