#ifndef __MultiGrid_h__
#define __MultiGrid_h__ 1

#include "Array.h"
#include "Tridiagonal.h"

namespace Array {

class BC {
 protected:
  int internal,external;
  int offset; // Offset in origin
  int ioff;   // Index of first ghost point
  int active; // Index of first active point
  int poffset;// Offset to first prolongation point
 public:	
  BC() : offset(0), ioff(0), active(1), poffset(-1) {}
  int Internal() const {return internal;}
  int External() const {return external;}
  int Offset() const {return offset;}
  int Ioff() const {return ioff;}
  int Active() const {return active;}
  int Poffset() const {return poffset;}
  virtual int Resolution(int radix, int lvl) const =0;
};

class DirichletBC : public BC {
 public:	
  DirichletBC() {internal=2; external=0;}
  int Resolution(int radix, int lvl) const {return pow(radix,lvl+1)-1;}
};

class ExtendedDirichletBC : public BC {
 public:	
  ExtendedDirichletBC() {internal=2; external=2; active=2; offset=-1;}
  int Resolution(int radix, int lvl) const {return pow(radix,lvl)-1;}
};

class NeumannBC : public BC {
 public:	
  NeumannBC() {internal=0; external=2; offset=-1; poffset=0;}
  int Resolution(int radix, int lvl) const {return pow(radix,lvl)+1;}
};

class PeriodicBC : public BC {
 public:	
  PeriodicBC() {internal=1; external=1;}
  int Resolution(int radix, int lvl) const {return pow(radix,lvl+1);}
};

class DirichletBC2 : public BC {
 public:	
  DirichletBC2() {internal=2; external=2; ioff=-1;}
  int Resolution(int radix, int lvl) const {return pow(radix,lvl+1)-1;}
};

class ExtendedDirichletBC2 : public BC {
 public:	
  ExtendedDirichletBC2() {
    internal=2; external=4; ioff=-1; active=2; offset=-1;
  }
  int Resolution(int radix, int lvl) const {return pow(radix,lvl)-1;}
};

class NeumannBC2 : public BC {
 public:	
  NeumannBC2() {internal=0; external=4; ioff=-1; offset=-1; poffset=0;}
  int Resolution(int radix, int lvl) const {return pow(radix,lvl)+1;}
};

class PeriodicBC2 : public BC {
 public:	
  PeriodicBC2() {internal=1; external=3; ioff=-1;}
  int Resolution(int radix, int lvl) const {return pow(radix,lvl+1);}
};

const DirichletBC Dirichlet[1];
const DirichletBC2 Dirichlet2[1];

const ExtendedDirichletBC ExtendedDirichlet[1];
const ExtendedDirichletBC2 ExtendedDirichlet2[1];

const NeumannBC Neumann[1];
const NeumannBC2 Neumann2[1];

const PeriodicBC Periodic[1];
const PeriodicBC2 Periodic2[1];

const PeriodicBC MixedB[1];
const PeriodicBC2 MixedB2[1];

class Limits {
 public:
  Real min, max;
  const BC *bc;
  int skiplevels; // number of levels to skip in this direction
  int n0; // number of points in lowest level
  Limits() {};
  Limits(Real min_, Real max_, const BC *bc_=Dirichlet, int skiplevels_=0,
	 int n0_=1) :
    min(min_), max(max_), bc(bc_), skiplevels(skiplevels_), n0(n0_) {
    if(skiplevels < 0) skiplevels=0;
  }
};

template<class T, class V>
class Grid {
 protected:
  int level;
  int homogeneous;
  Grid<T,V> *parent;
  int nonlinear;
  int radix;
  int dimension;
  T v, v2, d;
  int allpoints; // total number of points in grid (incl. boundary points)
 public:
  virtual ~Grid() {};
  int AllPoints() {return allpoints;}
  virtual void Allocate(int allocate=1)=0;
	
  virtual void Initialize(int level0, int homogeneous0=0,
			  Grid<T,V> *parent0=NULL,
			  int nonlinear0=0, int allocate=1) {
    level=level0; homogeneous=homogeneous0;
    nonlinear=nonlinear0; parent=parent0;
    allpoints=1;
    Allocate(allocate);
    if(!allocate) return;
    d=(V) 0.0;
    if(level > 0) v=(V) 0.0;
  }
	
  virtual void Mesh(Array1<Real>::opt &x, Limits limits, int& n,
		    int& n1bc, int& nbc, Real& h, Real& hinv,
		    Real& h2, Real& h2inv, int& r, int& offset, int &ioff,
		    int& start, int& startp, int& stop, int& stopp) {
    // number of points in one direction
    int lvl=max(level-limits.skiplevels,0);
    n=limits.n0*limits.bc->Resolution(radix,lvl);
    int n1=(lvl > 0) ? limits.n0*limits.bc->Resolution(radix,lvl-1) : n;
    int bcpts=limits.bc->Internal()+limits.bc->External();
    nbc=n+bcpts;
    n1bc=n1+bcpts;
    allpoints *= nbc;
    h=(limits.max-limits.min)/(n+limits.bc->Internal()-1);
    hinv=1.0/h;
    h2=h*h; h2inv=hinv*hinv;
    ioff=limits.bc->Ioff();
#ifdef NDEBUG		
    x=new Real[nbc]-ioff;
#else		
    x.Allocate(nbc,ioff);
#endif		
    offset=limits.bc->Offset();
    for(int i=ioff; i < nbc+ioff; i++)
      x[i]=limits.min+(i+offset)*h;
    start=limits.bc->Active();
    if(level <= limits.skiplevels) {startp=start; r=1; offset=0;}
    else {startp=start+limits.bc->Poffset(); r=radix;}
    stop=limits.bc->Active()+n-1;
    stopp=limits.bc->Active()+n1-1;
  }
	
  virtual void Defect(const T& d0, const T& u, const T& f)=0;
  virtual void Smooth(const T& u, const T& f)=0;
  virtual void Restrict(const T& r, const T& u)=0;
  virtual void SubtractProlongation(const T& u, const T& v0)=0;
  virtual void BoundaryConditions(const T&)=0;
  virtual void L0inv(const T&, const T&)=0;
  virtual inline void SubtractKernel(const T&, const T&) {};
	
  int Solve(const T& u, const T& f, int nu1=0, int gamma=1, int nu2=1,
	    int singular=0, int niter=1, Real rate=0.0, int statistics=0,
	    V *pdefect0=NULL, V *pdefect=NULL) {
    // u and f must be distinct
    int i,it=0;
    if(level == 0) {L0inv(u,f); return 0;}
    if(level == 1) gamma=1;
		
    V defect=0.0, defect0=0.0;
		
    if(rate) statistics=1;
		
    if(statistics) {
      Defect(d,u,f);
      defect0=Norm(d);
    }
		
    for(;;) {
      for(i=0; i < nu1; i++) {Smooth(u,f); BoundaryConditions(u);}
			
      if(!statistics || nu1) Defect(d,u,f);

      int homogeneous0=homogeneous;
      homogeneous=1;
      BoundaryConditions(d);
      homogeneous=homogeneous0;

      Restrict(d,d);
		
      if(nonlinear) {
	Restrict(v,u);
	parent->BoundaryConditions(v);
	parent->Defect(d,v,d);
	v2=v;
	for(i=0; i < gamma; i++) parent->Solve(v2,d,nu1,gamma,nu2);
	v -= v2;
      } else {
	v=(V) 0.0;
	for(i=0; i < gamma; i++) parent->Solve(v,d,nu1,gamma,nu2);
      }		
		
      SubtractProlongation(u,v);
      BoundaryConditions(u);
      for(i=0; i < nu2; i++) {Smooth(u,f); BoundaryConditions(u);}
			
      it++;
			
      if(statistics) {
	Defect(d,u,f);
	defect=Norm(d);
	if(rate && rate*defect <= defect0) break;
      }
			
      if(it == niter) {
	if(rate) it=0;
	break;
      }
    }
		
		
    if(singular) {
      if(!statistics) Defect(d,u,f);
      SubtractKernel(u,d);
      BoundaryConditions(u);
      if(statistics) defect=Norm(d);
    }
		
    if(statistics) {
      if(pdefect0) *pdefect0=defect0;
      if(pdefect) *pdefect=defect;
    }
		
    return it;
  }
	
  virtual void Lu(const T& u, const T& f) { // u and f must be distinct
    d=(V) 0.0;
    Defect(f,u,d);
  }
	
  void ComputeForce(const T& u, const T& f) {
    BoundaryConditions(u);
    Lu(u,f);
  }
	
  virtual void ReportValue(V error) {
    cout << error;
  }
	
  virtual void ReportValue(V error, V lasterror) {
    cout << error << "\t" << divide0(lasterror,error);
  }
	
  virtual void ReportHeader() {
    cout << "iter\tdefect\tratio" << endl;
  }
	
  virtual void ReportHeaderError() {
    cout << "iter\tdefect\tratio\terror\tratio" << endl;
  }
	
  virtual void Report(V defect, int it) {
    cout << it << "\t";
    ReportValue(defect);
    cout << endl;
  }
	
  virtual void Report(V defect, V defect0, int it) {
    cout << it << "\t";
    ReportValue(defect,defect0);
    cout << endl;
  }
	
  virtual void Report(V defect, int it, V error) {
    cout << it << "\t";
    ReportValue(defect);
    cout << "\t\t";
    ReportValue(error);
    cout << endl;
  }
	
  virtual void Report(V defect, V defect0, int it, V error, V error0) {
    cout << it << "\t";
    ReportValue(defect,defect0);
    cout << "\t\t";
    ReportValue(error,error0);
    cout << endl;
  }
	
  virtual void Sum2(const T& u, V& s)=0;
  V Deviation(V sum) {return sqrt(sum/(allpoints-1));}
  V Norm(const T& u) {V s=0.0; Sum2(u,s); return Deviation(s);}
  V DefectNorm(const T& u, const T& f) {Defect(d,u,f); return Norm(d);}
	
  void Compare(const T& u, const T& uexact, const T& f,
	       V& defect, V& error, V offset=0.0) {
    V s=0.0;
    defect=DefectNorm(u,f);
    for (int i=0; i < allpoints; i++) s += abs2(u(i)-uexact(i)-offset);
    error=Deviation(s);
  }
};

template<class G>
class MultiGrid {
 protected:
  G *grid;
  int nlevel;
 public:
  MultiGrid(int nlevel0, int nonlinear=0) : nlevel(nlevel0) {
    grid=new G[nlevel];
    for(int i=0; i < nlevel; i++) 
      grid[i].Initialize(i,i < nlevel-1,
			 (i == 0 ? NULL : grid+i-1),nonlinear);
  }
	
  G& Grid(int i) {return grid[i];}
  G& Fine() {return grid[nlevel-1];}
};

}

#endif
