#ifndef __Grid1_h__
#define __Grid1_h__ 1

#include "MultiGrid.h"

template<class T>
class Grid1 : public Grid<Array1<T>,T> {
protected:	
	// number of points in previous and current levels and incl. boundaries
	int nx1, nx, nx1bc, nxbc, sx, rx, offx;
	Array1(Real) x;
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
	Real Hx() {return hx;}
	
	void Allocate(int allocate=1) {
		Mesh(x,XMeshRange(),nx1,nx,nx1bc,nxbc,hx,hxinv,hx2,hx2inv,sx,rx,offx);
		if(!allocate) return;
		d.Allocate(nxbc);
		if(level > 0) {
			v.Allocate(nx1bc);
			if(nonlinear) v2.Allocate(nx1bc);
		}
	}

	virtual void Defect(const Array1<T>& d0, const Array1<T>& u, const
						Array1<T>& f)=0; 
	virtual void Smooth(const Array1<T>& u, const Array1<T>& f)=0;
	
	virtual void GaussSeidel(const Array1<T>&, const Array1<T>&, int, int) {};
	
	void Restrict(const Array1<T>& r, const Array1<T>& u) {
		if(&r != &u) XDirichlet(r,u,1);
		Array1(T) u0=u+offx;
		for(int i=1; i <= nx1; i++)
			r[i]=0.5*(0.5*(u0[rx*i-1]+u0[rx*i+1])+u0[rx*i]);
	}
	
	void SubtractProlongation(const Array1<T>& u, const Array1<T>& v0) {
		Array1(T) u0=u+offx;
		for(int i=sx; i <= nx1; i++) {
			u0[rx*i] -= v0[i];
			u0[rx*i+1] -= 0.5*(v0[i]+v0[i+1]);
		}
	}
	
	void Jacobi(const Array1<T>& u, const Array1<T>& f, Real omegah2) {
		Defect(d,u,f);
		for(int i=1; i <= nx; i++) u[i] -= omegah2*d[i];
	}
	
	void Lexicographical(const Array1<T>& u, const Array1<T>& f) {
		GaussSeidel(u,f,1,1);
	}
	
	void RedBlack(const Array1<T>& u, const Array1<T>& f) {
		GaussSeidel(u,f,1,2);
		GaussSeidel(u,f,2,2);
	}
	
	void Sum2(const Array1<T>& u, T& s) {
		for(int i=1; i <= nx; i++) s += abs2(u[i]);
	}

	inline virtual void BoundaryConditions(const Array1<T>& u)=0;
	
	void XDirichlet(const Array1<T>&) {}
	
	void XDirichlet(const Array1<T>& u, T b0, T b1) {
		if(homogenous) return;
		u[0]=b0;
		u[nx+1]=b1;
	}
	
	void XDirichlet(const Array1<T>& u, const Array1<T>& b, int contract=0) {
		int nx0;
		if(contract) nx0=nx1;
		else nx0=nx;
		u[0]=b[0];
		u[nx0+1]=b[nx+1];
	}
	
	void XNeumann(const Array1<T>& u) {
		u[0]=u[2];
		u[nx+1]=u[nx-1];
	}
	
	void XDirichletInterpolate(const Array1<T>& u, T b0, T b1) {
		if(homogeneous) {b0=b1=0.0;}
		else {b0 *= 2.0; b1 *= 2.0;}
		u[0]=b0-u[2];
		u[nx+1]=b1-u[nx-1];
	}
	
	void XConstant(const Array1<T>& u) {
		u[0]=u[1];
		u[nx+1]=u[nx];
	}
	
	void XMixed0(const Array1<T>& u) {
		u[0]=u[2];
	}
	
	void XMixed1(const Array1<T>& u) {
		u[nx+1]=u[nx-1];
	}
	
	void XPeriodic(const Array1<T>& u) {
		u[0]=u[nx];
		u[nx+1]=u[1];
	}
};

#endif
