#ifndef __Grid2_h__
#define __Grid2_h__ 1

#include "MultiGrid.h"

template<class T>
class Grid2 : public Grid<Array2<T>,T> {
protected:	
	// number of points in previous and current levels and incl. boundaries
	int nx1, nx, nx1bc, nxbc, sx, rx, offx, ox;
	int ny1, ny, ny1bc, nybc, sy, ry, offy, oy;
	Array1(Real) x;
	Array1(Real) y;
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
	Real Hx() {return hx;}
	Real Hy() {return hy;}

	void Allocate(int allocate=1) {
		Mesh(x,XMeshRange(),nx1,nx,nx1bc,nxbc,hx,hxinv,hx2,hx2inv,sx,rx,
			 offx,ox);
		Mesh(y,YMeshRange(),ny1,ny,ny1bc,nybc,hy,hyinv,hy2,hy2inv,sy,ry,
			 offy,oy);
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
		for(int i=1; i <= nx1; i++) {
			int i2=rx*i+offx;
			Array1(T) ri=r[i];
			Array1(T) um=u[i2-1]+offy;
			Array1(T) uz=u[i2]+offy;
			Array1(T) up=u[i2+1]+offy;
			for(int j=1; j <= ny1; j++) {
				ri[j]=0.25*(0.5*(0.5*(um[ry*j-1]+um[ry*j+1]+up[ry*j-1]+
									  up[ry*j+1])+
								 um[ry*j]+uz[ry*j-1]+uz[ry*j+1]+up[ry*j])+
							uz[ry*j]);
			}
		}
	}
	
	virtual void SubtractProlongation(const Array2<T>& u,
									  const Array2<T>& v0) {
		for(int i=sx; i <= nx1; i++) {
			int i2=rx*i+offx;
			Array1(T) vz=v0[i];
			Array1(T) vp=v0[i+1];
			Array1(T) uz=u[i2]+offy;
			Array1(T) up=u[i2+1]+offy;
			for(int j=sy; j <= ny1; j++) {
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
		for(int i=1; i <= nx; i++) {
			Array1(T) di=d[i];
			Array1(T) ui=u[i];
			for(int j=1; j <= ny; j++) {
					ui[j] -= omegah2*di[j];
			}
		}
	}
	
	void Lexicographical(const Array2<T>& u, const Array2<T>& f) {
		GaussSeidel(u,f,1,1,1,1);
	}
	
	void XLexicographical(const Array2<T>& u, const Array2<T>& f) {
		XGaussSeidel(u,f,1,1);
	}
	
	void YLexicographical(const Array2<T>& u, const Array2<T>& f) {
		YGaussSeidel(u,f,1,1);
	}
	
	void RedBlack(const Array2<T>& u, const Array2<T>& f) {
		GaussSeidel(u,f,1,1,2,2);
		GaussSeidel(u,f,2,2,2,2);
		
		GaussSeidel(u,f,1,2,2,2);
		GaussSeidel(u,f,2,1,2,2);
	}
	
	void XZebra(const Array2<T>& u, const Array2<T>& f) {
		XGaussSeidel(u,f,1,2);
		XGaussSeidel(u,f,2,2);
	}
	
	void YZebra(const Array2<T>& u, const Array2<T>& f) {
		YGaussSeidel(u,f,1,2);
		YGaussSeidel(u,f,2,2);
	}
	
	void Sum2(const Array2<T>& u, T& s) {
		for(int i=1; i <= nx; i++) {
			Array1(T) ui=u[i];
			for(int j=1; j <= ny; j++) s += abs2(ui[j]);
		}
	}
	
	virtual inline void BoundaryConditions(const Array2<T>& u)=0;
	
	void XDirichlet(const Array2<T>&) {}
	
	void XDirichlet2(const Array2<T>&) {}
		
	void XDirichlet(const Array2<T>& u, T b0, T b1) {
		if(homogeneous) return;
		Array1(T) u0=u[0];
		Array1(T) unx1=u[nx+1];
		for(int j=0; j < nybc; j++) {
			u0[j]=b0;
			unx1[j]=b1;
		}
	}
	
	void XDirichlet2(const Array2<T>& u, T b0, T b1) {
		if(homogeneous) return;
		Array1(T) um1=u[-1];
		Array1(T) u0=u[0];
		Array1(T) unx1=u[nx+1];
		Array1(T) unx2=u[nx+2];
		int nybco=nybc+oy;
		for(int j=oy; j < nybco; j++) {
			um1[j]=b0;
			u0[j]=b0;
			unx1[j]=b1;
			unx2[j]=b1;
		}
	}
	
	void XDirichlet(const Array2<T>& u, const Array2<T>& b, int contract=0) {
		int nx0,ny0bc;
		if(contract) {nx0=nx1; ny0bc=ny1bc;} 
		else {nx0=nx; ny0bc=nybc;} 
		Array1(T) u0=u[0];
		Array1(T) unx1=u[nx0+1];
		Array1(T) b0=b[0];
		Array1(T) bnx1=b[nx+1];
		for(int j=0; j < ny0bc; j++) {
			u0[j]=b0[j];
			unx1[j]=bnx1[j];
		}
	}
	
	void XDirichlet2(const Array2<T>& u, const Array2<T>& b, int contract=0) {
		int nx0,ny0bc;
		if(contract) {nx0=nx1; ny0bc=ny1bc;} 
		else {nx0=nx; ny0bc=nybc;} 
		Array1(T) um1=u[-1];
		Array1(T) u0=u[0];
		Array1(T) unx1=u[nx0+1];
		Array1(T) unx2=u[nx0+2];
		Array1(T) bm1=b[-1];
		Array1(T) b0=b[0];
		Array1(T) bnx1=b[nx+1];
		Array1(T) bnx2=b[nx+2];
		int ny0bco=ny0bc+oy;
		for(int j=oy; j < ny0bco; j++) {
			um1[j]=bm1[j];
			u0[j]=b0[j];
			unx1[j]=bnx1[j];
			unx2[j]=bnx2[j];
		}
	}
	
	void XNeumann(const Array2<T>& u) {
		Array1(T) u0=u[0];
		Array1(T) u2=u[2];
		Array1(T) unxm1=u[nx-1];
		Array1(T) unx1=u[nx+1];
		for(int j=0; j < nybc; j++) {
			u0[j]=u2[j];
			unx1[j]=unxm1[j];
		}
	}
	
	void XNeumann2(const Array2<T>& u) {
		Array1(T) um1=u[-1];
		Array1(T) u0=u[0];
		Array1(T) u2=u[2];
		Array1(T) u3=u[3];
		Array1(T) unxm2=u[nx-2];
		Array1(T) unxm1=u[nx-1];
		Array1(T) unx1=u[nx+1];
		Array1(T) unx2=u[nx+2];
		int nybco=nybc+oy;
		for(int j=oy; j < nybco; j++) {
			um1[j]=u3[j];
			u0[j]=u2[j];
			unx1[j]=unxm1[j];
			unx2[j]=unxm2[j];
		}
	}
	
	void XDirichletInterpolate(const Array2<T>& u, T b0, T b1) {
		if(homogeneous) {b0=b1=0.0;}
		else {b0 *= 2.0; b1 *= 2.0;}
		Array1(T) u0=u[0];
		Array1(T) u2=u[2];
		Array1(T) unxm1=u[nx-1];
		Array1(T) unx1=u[nx+1];
		for(int j=0; j < nybc; j++) {
			u0[j]=b0-u2[j];
			unx1[j]=b1-unxm1[j];
		}
	}
	
	void XDirichletInterpolate2(const Array2<T>& u, T b0, T b1) {
		if(homogeneous) {b0=b1=0.0;}
		else {b0 *= 2.0; b1 *= 2.0;}
		Array1(T) um1=u[-1];
		Array1(T) u0=u[0];
		Array1(T) u2=u[2];
		Array1(T) u3=u[3];
		Array1(T) unxm2=u[nx-2];
		Array1(T) unxm1=u[nx-1];
		Array1(T) unx1=u[nx+1];
		Array1(T) unx2=u[nx+2];
		int nybco=nybc+oy;
		for(int j=oy; j < nybco; j++) {
			um1[j]=b0-u3[j];
			u0[j]=b0-u2[j];
			unx1[j]=b1-unxm1[j];
			unx2[j]=b1-unxm2[j];
		}
	}
	
	void XConstant(const Array2<T>& u) {
		Array1(T) u0=u[0];
		Array1(T) u1=u[1];
		Array1(T) unx=u[nx];
		Array1(T) unx1=u[nx+1];
		for(int j=0; j < nybc; j++) {
			u0[j]=u1[j];
			unx1[j]=unx[j];
		}
	}
	
	void XConstant2(const Array2<T>& u) {
		Array1(T) um1=u[-1];
		Array1(T) u0=u[0];
		Array1(T) u1=u[1];
		Array1(T) unx=u[nx];
		Array1(T) unx1=u[nx+1];
		Array1(T) unx2=u[nx+2];
		int nybco=nybc+oy;
		for(int j=oy; j < nybco; j++) {
			u0[j]=um1[j]=u1[j];
			unx2[j]=unx1[j]=unx[j];
		}
	}
	
	void XMixedA(const Array2<T>& u) {
		Array1(T) u0=u[0];
		Array1(T) u2=u[2];
		for(int j=0; j < nybc; j++) {
			u0[j]=u2[j];
		}
	}
	
	void XMixedA2(const Array2<T>& u) {
		Array1(T) um1=u[-1];
		Array1(T) u0=u[0];
		Array1(T) u2=u[2];
		Array1(T) u3=u[3];
		int nybco=nybc+oy;
		for(int j=oy; j < nybco; j++) {
			um1[j]=u3[j];
			u0[j]=u2[j];
		}
	}
	
	void XMixedB(const Array2<T>& u) {
		Array1(T) unxm1=u[nx-1];
		Array1(T) unx1=u[nx+1];
		for(int j=0; j < nybc; j++) {
			unx1[j]=unxm1[j];
		}
	}
	
	void XMixedB2(const Array2<T>& u) {
		Array1(T) unxm2=u[nx-2];
		Array1(T) unxm1=u[nx-1];
		Array1(T) unx1=u[nx+1];
		Array1(T) unx2=u[nx+2];
		int nybco=nybc+oy;
		for(int j=oy; j < nybco; j++) {
			unx1[j]=unxm1[j];
			unx2[j]=unxm2[j];
		}
	}
	
	void XPeriodic(const Array2<T>& u) {
		Array1(T) u0=u[0];
		Array1(T) u1=u[1];
		Array1(T) unx=u[nx];
		Array1(T) unx1=u[nx+1];
		for(int j=0; j < nybc; j++) {
			u0[j]=unx[j];
			unx1[j]=u1[j];
		}
	}
	
	void XPeriodic2(const Array2<T>& u) {
		Array1(T) um1=u[-1];
		Array1(T) u0=u[0];
		Array1(T) u1=u[1];
		Array1(T) u2=u[2];
		Array1(T) unxm1=u[nx-1];
		Array1(T) unx=u[nx];
		Array1(T) unx1=u[nx+1];
		Array1(T) unx2=u[nx+2];
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
		for(int i=0; i < nxbc; i++) {
			Array1(T) ui=u[i];
			ui[0]=b0;
			ui[ny+1]=b1;
		}
	}
	
	void YDirichlet2(const Array2<T>& u, T b0, T b1) {
		if(homogeneous) return;
		int nxbco=nxbc+ox;
		for(int i=ox; i < nxbco; i++) {
			Array1(T) ui=u[i];
			ui[-1]=b0;
			ui[0]=b0;
			ui[ny+1]=b1;
			ui[ny+2]=b1;
		}
	}
	
	void YDirichlet(const Array2<T>& u, const Array2<T>& b, int contract=0) {
		int nx0bc,ny0;
		if(contract) {nx0bc=nx1bc; ny0=ny1;} 
		else {nx0bc=nxbc; ny0=ny;} 
		for(int i=0; i < nx0bc; i++) {
			Array1(T) ui=u[i];
			Array1(T) bi=b[i];
			ui[0]=bi[0];
			ui[ny0+1]=bi[ny+1];
		}
	}
	
	void YDirichlet2(const Array2<T>& u, const Array2<T>& b, int contract=0) {
		int nx0bc,ny0;
		if(contract) {nx0bc=nx1bc; ny0=ny1;} 
		else {nx0bc=nxbc; ny0=ny;} 
		int nx0bco=nx0bc+ox;
		for(int i=ox; i < nx0bco; i++) {
			Array1(T) ui=u[i];
			Array1(T) bi=b[i];
			ui[-1]=bi[-1];
			ui[0]=bi[0];
			ui[ny0+1]=bi[ny+1];
			ui[ny0+2]=bi[ny+2];
		}
	}
	
	void YNeumann(const Array2<T>& u) {
		for(int i=0; i < nxbc; i++) {
			Array1(T) ui=u[i];
			ui[0]=ui[2];
			ui[ny+1]=ui[ny-1];
		}
	}
	
	void YNeumann2(const Array2<T>& u) {
		int nxbco=nxbc+ox;
		for(int i=ox; i < nxbco; i++) {
			Array1(T) ui=u[i];
			ui[-1]=ui[3];
			ui[0]=ui[2];
			ui[ny+1]=ui[ny-1];
			ui[ny+2]=ui[ny-2];
		}
	}
	
	void YDirichletInterpolate(const Array2<T>& u, T b0, T b1) {
		if(homogeneous) {b0=b1=0.0;}
		else {b0 *= 2.0; b1 *= 2.0;}
		for(int i=0; i < nxbc; i++) {
			Array1(T) ui=u[i];
			ui[0]=b0-ui[2];
			ui[ny+1]=b1-ui[ny-1];
		}
	}
	
	void YDirichletInterpolate2(const Array2<T>& u, T b0, T b1) {
		if(homogeneous) {b0=b1=0.0;}
		else {b0 *= 2.0; b1 *= 2.0;}
		int nxbco=nxbc+ox;
		for(int i=ox; i < nxbco; i++) {
			Array1(T) ui=u[i];
			ui[-1]=b0-ui[3];
			ui[0]=b0-ui[2];
			ui[ny+1]=b1-ui[ny-1];
			ui[ny+2]=b1-ui[ny-2];
		}
	}
	
	void YConstant(const Array2<T>& u) {
		for(int i=0; i < nxbc; i++) {
			Array1(T) ui=u[i];
			ui[0]=ui[1];
			ui[ny+1]=ui[ny];
		}
	}
	
	void YConstant2(const Array2<T>& u) {
		int nxbco=nxbc+ox;
		for(int i=ox; i < nxbco; i++) {
			Array1(T) ui=u[i];
			ui[0]=ui[-1]=ui[1];
			ui[ny+2]=ui[ny+1]=ui[ny];
		}
	}
	
	void YPeriodic(const Array2<T>& u) {
		for(int i=0; i < nxbc; i++) {
			Array1(T) ui=u[i];
			ui[0]=ui[ny];
			ui[ny+1]=ui[1];
		}
	}
	
	void YPeriodic2(const Array2<T>& u) {
		int nxbco=nxbc+ox;
		for(int i=ox; i < nxbco; i++) {
			Array1(T) ui=u[i];
			ui[-1]=ui[ny-1];
			ui[0]=ui[ny];
			ui[ny+1]=ui[1];
			ui[ny+2]=ui[2];
		}
	}
};

#endif
