#ifndef __Grid3_h__
#define __Grid3_h__ 1

#include "MultiGrid.h"

template<class T>
class Grid3 : public Grid<Array3<T>,T> {
protected:	
	// number of points in previous and current levels and incl. boundaries
	int nx1, nx, nx1bc, nxbc, sx, rx, offx, ox;
	int ny1, ny, ny1bc, nybc, sy, ry, offy, oy;
	int nz1, nz, nz1bc, nzbc, sz, rz, offz, oz;
	Array1(Real) x;
	Array1(Real) y;
	Array1(Real) z;
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
	Real Hx() {return hx;}
	Real Hy() {return hy;}
	Real Hz() {return hz;}

	void Allocate(int allocate=1) {
		Mesh(x,XMeshRange(),nx1,nx,nx1bc,nxbc,hx,hxinv,hx2,hx2inv,sx,rx,
			 offx,ox);
		Mesh(y,YMeshRange(),ny1,ny,ny1bc,nybc,hy,hyinv,hy2,hy2inv,sy,ry,
			 offy,oy);
		Mesh(z,ZMeshRange(),nz1,nz,nz1bc,nzbc,hz,hzinv,hz2,hz2inv,sz,rz,
			 offz,oz);
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
		for(int i=1; i <= nx1; i++) {
			int i2=rx*i+offx;
			Array2<T> ri=r[i], um=u[i2-1], uz=u[i2], up=u[i2+1];
			for(int j=1; j <= ny1; j++) {
				int j2=ry*j+offy;
				Array1(T) rij=ri[j];
				Array1(T) umm=um[j2-1]+offz;
				Array1(T) umz=um[j2]+offz;
				Array1(T) ump=um[j2+1]+offz;
				Array1(T) uzm=uz[j2-1]+offz;
				Array1(T) uzz=uz[j2]+offz;
				Array1(T) uzp=uz[j2+1]+offz;
				Array1(T) upm=up[j2-1]+offz;
				Array1(T) upz=up[j2]+offz;
				Array1(T) upp=up[j2+1]+offz;
				for(int k=1; k <= nz1; k++) {
					rij[k]=0.125*(0.5*(0.5*(0.5*(
						umm[rz*k-1]+umm[rz*k+1]+ump[rz*k-1]+
						ump[rz*k+1]+um[rz*j-1][rz*k-1]+umm[rz*k+1]+
						upp[rz*k-1]+upp[rz*k+1])+
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
		for(int i=sx; i <= nx1; i++) {
			int i2=rx*i+offx;
			Array2<T> uz=u[i2], up=u[i2+1], vz=v0[i], vp=v0[i+1];
			for(int j=sy; j <= ny1; j++) {
				int j2=ry*j+offy;
				Array1(T) uzz=uz[j2]+offz;
				Array1(T) uzp=uz[j2+1]+offz;
				Array1(T) upz=up[j2]+offz;
				Array1(T) upp=up[j2+1]+offz;
				Array1(T) vzz=vz[j];
				Array1(T) vzp=vz[j+1];
				Array1(T) vpz=vp[j];
				Array1(T) vpp=vp[j+1];
				for(int k=sz; k <= nz1; k++) {
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
		for(int i=1; i <= nx; i++) {
			Array2<T> di=d[i], ui=u[i];
			for(int j=1; j <= ny; j++) {
				Array1(T) dij=di[j];
				Array1(T) uij=ui[j];
				for(int k=1; k <= nz; k++) {
					uij[k] -= omegah2*dij[k];
				}
			}
		}
	}
	
	void Lexicographical(const Array3<T>& u, const Array3<T>& f) {
		GaussSeidel(u,f,1,1,1,1,1,1);
	}
	
	void XLexicographical(const Array3<T>& u, const Array3<T>& f) {
		XGaussSeidel(u,f,1,1,1,1);
	}
	
	void YLexicographical(const Array3<T>& u, const Array3<T>& f) {
		YGaussSeidel(u,f,1,1,1,1);
	}
	
	void ZLexicographical(const Array3<T>& u, const Array3<T>& f) {
		ZGaussSeidel(u,f,1,1,1,1);
	}
	
	void XYLexicographical(const Array3<T>& u, const Array3<T>& f) {
		XYGaussSeidel(u,f,1,1);
	}
	
	void XZLexicographical(const Array3<T>& u, const Array3<T>& f) {
		XZGaussSeidel(u,f,1,1);
	}
	
	void YZLexicographical(const Array3<T>& u, const Array3<T>& f) {
		YZGaussSeidel(u,f,1,1);
	}
	
	void RedBlack(const Array3<T>& u, const Array3<T>& f) {
		GaussSeidel(u,f,1,1,1,2,2,2);
		GaussSeidel(u,f,1,2,2,2,2,2);
		GaussSeidel(u,f,2,1,2,2,2,2);
		GaussSeidel(u,f,2,2,1,2,2,2);
		
		GaussSeidel(u,f,1,1,2,2,2,2);
		GaussSeidel(u,f,1,2,1,2,2,2);
		GaussSeidel(u,f,2,1,1,2,2,2);
		GaussSeidel(u,f,2,2,2,2,2,2);
	}
	
	void XZebra(const Array3<T>& u, const Array3<T>& f) {
		XGaussSeidel(u,f,1,1,2,2);
		XGaussSeidel(u,f,2,2,2,2);
		XGaussSeidel(u,f,1,2,2,2);
		XGaussSeidel(u,f,2,1,2,2);
	}
	
	void YZebra(const Array3<T>& u, const Array3<T>& f) {
		YGaussSeidel(u,f,1,1,2,2);
		YGaussSeidel(u,f,2,2,2,2);
		YGaussSeidel(u,f,1,2,2,2);
		YGaussSeidel(u,f,2,1,2,2);
	}
	
	void ZZebra(const Array3<T>& u, const Array3<T>& f) {
		ZGaussSeidel(u,f,1,1,2,2);
		ZGaussSeidel(u,f,2,2,2,2);
		ZGaussSeidel(u,f,1,2,2,2);
		ZGaussSeidel(u,f,2,1,2,2);
	}
	
	void XYZebra(const Array3<T>& u, const Array3<T>& f) {
		XYGaussSeidel(u,f,1,2);
		XYGaussSeidel(u,f,2,2);
	}
	
	void XZZebra(const Array3<T>& u, const Array3<T>& f) {
		XZGaussSeidel(u,f,1,2);
		XZGaussSeidel(u,f,2,2);
	}
	
	void YZZebra(const Array3<T>& u, const Array3<T>& f) {
		YZGaussSeidel(u,f,1,2);
		YZGaussSeidel(u,f,2,2);
	}
	
	void Sum2(const Array3<T>& u, T& s) {
		for(int i=1; i <= nx; i++) {
			Array2<T> ui=u[i];
			for(int j=1; j <= ny; j++) {
				Array1(T) uij=ui[j];
				for(int k=1; k <= nz; k++) s += abs2(uij[k]);
			}
		}
	}
	
	virtual inline void BoundaryConditions(const Array3<T>& u)=0;
	
	void XDirichlet(const Array3<T>&) {}
	
	void XDirichlet(const Array3<T>& u, T b0, T b1) {
		if(homogeneous) return;
		Array2<T> u0=u[0], unx1=u[nx+1];
		for(int j=0; j < nybc; j++) {
			Array1(T) u0j=u0[j];
			Array1(T) unx1j=unx1[j];
			for(int k=0; k < nzbc; k++) {
				u0j[k]=b0;
				unx1j[k]=b1;
			}
		}
	}
	
	void XDirichlet(const Array3<T>& u, const Array3<T>& b, int contract=0) {
		int nx0,ny0bc,nz0bc;
		if(contract) {nx0=nx1; ny0bc=ny1bc; nz0bc=nz1bc;}
		else {nx0=nx; ny0bc=nybc; nz0bc=nzbc;} 
		Array2<T> u0=u[0], unx1=u[nx0+1];
		Array2<T> b0=b[0], bnx1=b[nx+1];
		for(int j=0; j < ny0bc; j++) {
			Array1(T) u0j=u0[j];
			Array1(T) unx1j=unx1[j];
			Array1(T) b0j=b0[j];
			Array1(T) bnx1j=bnx1[j];
			for(int k=0; k < nz0bc; k++) {
				u0j[k]=b0j[k];
				unx1j[k]=bnx1j[k];
			}
		}
	}
	
	void XNeumann(const Array3<T>& u) {
		Array2<T> u0=u[0], u2=u[2], unxm1=u[nx-1], unx1=u[nx+1];
		for(int j=0; j < nybc; j++) {
			Array1(T) u0j=u0[j];
			Array1(T) u2j=u2[j];
			Array1(T) unxm1j=unxm1[j];
			Array1(T) unx1j=unx1[j];
			for(int k=0; k < nzbc; k++) {
				u0j[k]=u2j[k];
				unx1j[k]=unxm1j[k];
			}
		}
	}
	
	void XDirichletInterpolate(const Array3<T>& u, T b0, T b1) {
		if(homogeneous) {b0=b1=0.0;}
		else {b0 *= 2.0; b1 *= 2.0;}
		Array2<T> u0=u[0], u2=u[2], unxm1=u[nx-1], unx1=u[nx+1];
		for(int j=0; j < nybc; j++) {
			Array1(T) u0j=u0[j];
			Array1(T) u2j=u2[j];
			Array1(T) unxm1j=unxm1[j];
			Array1(T) unx1j=unx1[j];
			for(int k=0; k < nzbc; k++) {
				u0j[k]=b0-u2j[k];
				unx1j[k]=b1-unxm1j[k];
			}
		}
	}
	
	void XConstant(const Array3<T>& u) {
		Array2<T> u0=u[0], u1=u[1], unx=u[nx], unx1=u[nx+1];
		for(int j=0; j < nybc; j++) {
			Array1(T) u0j=u0[j];
			Array1(T) u1j=u1[j];
			Array1(T) unxj=unx[j];
			Array1(T) unx1j=unx1[j];
			for(int k=0; k < nzbc; k++) {
				u0j[k]=u1j[k];
				unx1j[k]=unxj[k];
			}
		}
	}
	
	void XMixedA(const Array3<T>& u) {
		Array2<T> u0=u[0], u2=u[2];
		for(int j=0; j < nybc; j++) {
			Array1(T) u0j=u0[j];
			Array1(T) u2j=u2[j];
			for(int k=0; k < nzbc; k++) {
				u0j[k]=u2j[k];
			}
		}
	}
	
	void XMixedB(const Array3<T>& u) {
		Array2<T> unxm1=u[nx-1], unx1=u[nx+1];
		for(int j=0; j < nybc; j++) {
			Array1(T) unxm1j=unxm1[j];
			Array1(T) unx1j=unx1[j];
			for(int k=0; k < nzbc; k++) {
				unx1j[k]=unxm1j[k];
			}
		}
	}
	
	void XPeriodic(const Array3<T>& u) {
		Array2<T> u0=u[0], u1=u[1], unx=u[nx], unx1=u[nx+1];
		for(int j=0; j < nybc; j++) {
			Array1(T) u0j=u0[j];
			Array1(T) u1j=u1[j];
			Array1(T) unxj=unx[j];
			Array1(T) unx1j=unx1[j];
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
			Array1(T) ui0=ui[0];
			Array1(T) uiny1=ui[ny+1];
			for(int k=0; k < nzbc; k++) {
				ui0[k]=b0;
				uiny1[k]=b1;
			}
		}
	}
	
	void YDirichlet(const Array3<T>& u, const Array3<T>& b, int contract=0) {
		int nx0bc,ny0,nz0bc;
		if(contract) {nx0bc=nx1bc; ny0=ny1; nz0bc=nz1bc;} 
		else {nx0bc=nxbc; ny0=ny; nz0bc=nzbc;} 
		for(int i=0; i < nx0bc; i++) {
			Array2<T> ui=u[i], bi=b[i];
			Array1(T) ui0=ui[0];
			Array1(T) uiny1=ui[ny0+1];
			Array1(T) bi0=bi[0];
			Array1(T) biny1=bi[ny+1];
			for(int k=0; k < nz0bc; k++) {
				ui0[k]=bi0[k];
				uiny1[k]=biny1[k];
			}
		}
	}
	
	void YNeumann(const Array3<T>& u) {
		for(int i=0; i < nxbc; i++) {
			Array2<T> ui=u[i];
			Array1(T) ui0=ui[0];
			Array1(T) ui2=ui[2];
			Array1(T) uinym1=ui[ny-1];
			Array1(T) uiny1=ui[ny+1];
			for(int k=0; k < nzbc; k++) {
				ui0[k]=ui2[k];
				uiny1[k]=uinym1[k];
			}
		}
	}
	
	void YDirichletInterpolate(const Array3<T>& u, T b0, T b1) {
		if(homogeneous) {b0=b1=0.0;}
		else {b0 *= 2.0; b1 *= 2.0;}
		for(int i=0; i < nxbc; i++) {
			Array2<T> ui=u[i];
			Array1(T) ui0=ui[0];
			Array1(T) ui2=ui[2];
			Array1(T) uinym1=ui[ny-1];
			Array1(T) uiny1=ui[ny+1];
			for(int k=0; k < nzbc; k++) {
				ui0[k]=b0-ui2[k];
				uiny1[k]=b1-uinym1[k];
			}
		}
	}
	
	void YConstant(const Array3<T>& u) {
		for(int i=0; i < nxbc; i++) {
			Array2<T> ui=u[i];
			Array1(T) ui0=ui[0];
			Array1(T) ui1=ui[1];
			Array1(T) uiny=ui[ny];
			Array1(T) uiny1=ui[ny+1];
			for(int k=0; k < nzbc; k++) {
				ui0[k]=ui1[k];
				uiny1[k]=uiny[k];
			}
		}
	}
	
	void YPeriodic(const Array3<T>& u) {
		for(int i=0; i < nxbc; i++) {
			Array2<T> ui=u[i];
			Array1(T) ui0=ui[0];
			Array1(T) ui1=ui[1];
			Array1(T) uiny=ui[ny];
			Array1(T) uiny1=ui[ny+1];
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
				Array1(T) uij=ui[j];
				uij[0]=b0;
				uij[nz+1]=b1;
			}
		}
	}
	
	void ZDirichlet(const Array3<T>& u, const Array3<T>& b, int contract=0) {
		int nx0bc,ny0bc,nz0;
		if(contract) {nx0bc=nx1bc; ny0bc=ny1bc; nz0=nz1;} 
		else {nx0bc=nxbc; ny0bc=nybc; nz0=nz;} 
		for(int i=0; i < nx0bc; i++) {
			Array2<T> ui=u[i], bi=b[i];
			for(int j=0; j < ny0bc; j++) {
				Array1(T) uij=ui[j];
				Array1(T) bij=bi[j];
				uij[0]=bij[0];
				uij[nz0+1]=bij[nz+1];
			}
		}
	}
	
	void ZNeumann(const Array3<T>& u) {
		for(int i=0; i < nxbc; i++) {
			Array2<T> ui=u[i];
			for(int j=0; j < nybc; j++) {
				Array1(T) uij=ui[j];
				uij[0]=uij[2];
				uij[nz+1]=uij[nz-1];
			}
		}
	}
	
	void ZDirichletInterpolate(const Array3<T>& u, T b0, T b1) {
		if(homogeneous) {b0=b1=0.0;}
		else {b0 *= 2.0; b1 *= 2.0;}
		for(int i=0; i < nxbc; i++) {
			Array2<T> ui=u[i];
			for(int j=0; j < nybc; j++) {
				Array1(T) uij=ui[j];
				uij[0]=b0-uij[2];
				uij[nz+1]=b1-uij[nz-1];
			}
		}
	}
	
	void ZConstant(const Array3<T>& u) {
		for(int i=0; i < nxbc; i++) {
			Array2<T> ui=u[i];
			for(int j=0; j < nybc; j++) {
				Array1(T) uij=ui[j];
				uij[0]=uij[1];
				uij[nz+1]=uij[nz];
			}
		}
	}
	
	void ZPeriodic(const Array3<T>& u) {
		for(int i=0; i < nxbc; i++) {
			Array2<T> ui=u[i];
			for(int j=0; j < nybc; j++) {
				Array1(T) uij=ui[j];
				uij[0]=uij[nz];
				uij[nz+1]=uij[1];
			}
		}
	}
};

#endif
