#ifndef __Tridiagonal_h__
#define __Tridiagonal_h__ 1

#include <new>
#include "Array.h"

namespace Array {
  
template<class T>
inline void Allocate(T *&work, int &nwork, int n) 
{
  if(n > nwork) work=new(work,nwork=n) T;
}

// Solve the problem Lu=f for u given f, subject to the Dirichlet boundary
// condtions u[0]=u_0 and u[n+1]=u_{n+1},
// where L is the n x (n+2) matrix
//
// [b a b           ]
// [  b a b         ]
// [    b a b       ]
// [       ...      ]
// [         b a b  ]
// [           b a b]
//
// Note: u and f need not be distinct.

template<class T>
inline void Poisson1(int n, T *u, const T *f, Real b, Real a)
{
  if(n < 1) msg(ERROR,"Invalid matrix size (%d)",n);
	
  static int nwork=0;
  static Real *work=NULL;
  Allocate(work,nwork,n+1);
  int i;
	
  Real temp=b/a;
  work[1]=-temp;
  Real binv=1.0/b;
  u[1]=(f[1]*binv-u[0])*temp;
	
  for(i=2; i <= n; i++)	{
    Real temp=b/(a+b*work[i-1]);
    work[i]=-temp;
    u[i]=(f[i]*binv-u[i-1])*temp;
  }

  for(i=n; i >= 1; i--) u[i] += work[i]*u[i+1];
}

// Solve the problem Lu=f for u given f, subject to the Dirichlet boundary
// condtions u[0]=u_0 and u[n+1]=u_{n+1},
// where L is the n x (n+2) matrix
//
// [c1 a1 b1              ]
// [   c2 a2 b2           ]
// [      c3 a3 b3        ]
// [          ...         ]
// [              cn an bn]
//
// [work] is an optional work area of size Real [n+1] that may be set to a
// or b.
//
// Note: u and f need not be distinct.

template<class T>
inline void tridiagonal(int n, T *u, const T *f, const Real *c, const Real *a,
			const Real *b, Real *work)
{
  if(n < 1) msg(ERROR,"Invalid matrix size (%d)",n);
	
  static int nwork=0;
  static Real *work0=NULL;
  if(work == NULL) {Allocate(work0,nwork,n+1); work=work0;}
  int i;
	
  Real temp=1.0/a[1];
  work[1]=-b[1]*temp;
  u[1]=(f[1]-c[1]*u[0])*temp;
	
  for(i=2; i <= n; i++) {
    Real temp=1.0/(a[i]+c[i]*work[i-1]);
    work[i]=-b[i]*temp;
    u[i]=(f[i]-c[i]*u[i-1])*temp;
  }

  for(i=n; i >= 1; i--) u[i] += work[i]*u[i+1];
}

template<class T>
inline void tridiagonal(int n, T *u, const T *f, const Real *c, const Real *a,
			const Real *b)
{
  tridiagonal(n,u,f,c,a,b,(Real *) NULL);
}

// Solve multiple problems Lu=f for u given f and u[0] and u[n+1],
// where L is the n x (n+2) matrix
//
// [c1 a1 b1              ]
// [   c2 a2 b2           ]
// [      c3 a3 b3        ]
// [          ...         ]
// [              cn an bn]
//
// m is the number of tridiagonal systems to solve,
// inc1 is the stride between the elements of each vector,
// inc2 is the stride between the first elements of the vectors,
//
// [work] is an optional work area of size Real [n*inc1+m*inc2]
// that may be set to a or b.
//
// Note: u and f need not be distinct.

template<class T>
void mtridiagonal(int n, T *u, const T *f, const Real *c, const Real *a,
		  const Real *b, int m, int inc1, int inc2, Real *work)
{
  int jstop=m*inc2;
		
#ifndef _CRAYMVP
  if(inc1 == 1) {
    for(int j=0; j < jstop; j += inc2) 
      tridiagonal(n,u+j,f+j,c+j,a+j,b+j,work);
    return;
  }
#endif	
	
  if(n < 1) msg(ERROR,"Invalid matrix size (%d)",n);

  static int nwork=0;
  static Real *work0=NULL;
  if(work == NULL) {Allocate(work0,nwork,n*inc1+m*inc2); work=work0;}
  int i;
	
  const Real *ai=a+inc1, *bi=b+inc1, *ci=c+inc1;
  Real *worki=work+inc1;
  const T *fi=f+inc1;
  T *ui=u+inc1;
  for(int j=0; j < jstop; j += inc2) {
    Real temp=1.0/ai[j];
    worki[j]=-bi[j]*temp;
    ui[j]=(fi[j]-ci[j]*u[j])*temp;
  }
	
  for(i=2; i <= n; i++) {
    ai += inc1; bi += inc1; ci += inc1;
    Real *workim1=worki;
    T *uim1=ui;
    worki += inc1; ui += inc1; fi += inc1;
    for(int j=0; j < jstop; j += inc2) {
      Real temp=1.0/(ai[j]+ci[j]*workim1[j]);
      worki[j]=-bi[j]*temp;
      ui[j]=(fi[j]-ci[j]*uim1[j])*temp;
    }
  }

  T *uip1=ui+inc1;
  for(i=n; i >= 1; i--) {
    for(int j=0; j < jstop; j += inc2) {
      ui[j] += worki[j]*uip1[j];
    }
    worki -= inc1; uip1=ui; ui -= inc1;
  }
}

template<class T>
inline void mtridiagonal(int n, T *u, const T *f, const Real *c, const Real *a,
			 const Real *b, int m, int inc1, int inc2)
{
  mtridiagonal(n,u,f,c,a,b,m,inc1,inc2,(Real *) NULL);
}

// Solve the problem Lu=f for u given f, where L is the n x n matrix
//
// [ a b        b]
// [ b a b       ]
// [   b a b     ]
// [     ...     ]
// [        b a b]
// [ b        b a]
//
// On exit: u[0]=u[n], u[n+1]=u[1].
//
// Note: u and f need not be distinct.

template<class T>
void Poisson1p(int n, T *u, T *f, Real b, Real a)
{
  if(n < 1) msg(ERROR,"Invalid matrix size (%d)",n);
  if(n == 1) {u[2]=u[1]=u[0]=f[1]/a; return;}
  if(n == 2) {
    Real factor=1.0/(a*a-b*b);
    T temp=(a*f[1]-b*f[2])*factor;
    u[2]=u[0]=(-b*f[1]+a*f[2])*factor;
    u[3]=u[1]=temp;
    return;
  }
	
  int i;
  static int nsize=0;
  static Real *gamma=NULL, *delta=NULL;
	
  int size=n-1;
  if(size > nsize) {
    nsize=size;
    gamma=new(gamma,nsize) Real;
    delta=new(delta,nsize) Real;
  }
	
  Real binv=1.0/b;
  Real ainv=b/a;
  delta[1]=gamma[1]=ainv;
  Real alpha=a-b*ainv;
  u[1]=f[1]*binv*ainv;
  T fn=f[n]-b*u[1];
  Real beta=b;

  for(i=2; i <= n-2; i++)	{
    Real alphainv=b/(a-b*gamma[i-1]);
    beta *= -gamma[i-1];
    gamma[i]=alphainv;
    u[i]=(f[i]*binv-u[i-1])*alphainv;
    fn -= beta*u[i];
    delta[i]=-delta[i-1]*alphainv;
    alpha -= beta*delta[i];
  }
	
  Real alphainv=b/(a-b*gamma[n-2]);
  u[n-1]=(f[n-1]*binv-u[n-2])*alphainv;
  beta=b-beta*gamma[n-2];
  Real dnm1=(1.0-delta[n-2])*alphainv;
  T temp=u[0]=u[n]=(fn-beta*u[n-1])/(alpha-beta*dnm1);
  u[n-1] -= dnm1*temp;
	
  for(i=n-2; i >= 1; i--) u[i] -= gamma[i]*u[i+1]+delta[i]*temp;

  u[n+1]=u[1];
}

// Solve the problem Lu=f for u given f, where L is the n x n matrix
//
// [ -2b   b                   b ]
// [   b -2b   b                 ]
// [       b -2b   b             ]
// [              ...            ]
// [                   b -2b   b ]
// [   b                   b -2b ]
//
// On exit: u[0]=u[n], u[n+1]=u[1].
//
// For n > 2 this matrix is singular; the solution can only be determined
// to within an arbitrary constant.
//
// Note: u and f need not be distinct.

template<class T>
void Poisson1P(int n, T *u, T *f, Real b)
{
  if(n < 1) msg(ERROR,"Invalid matrix size (%d)",n);
  if(n == 1) {u[2]=u[1]=u[0]=-f[1]/(2.0*b); return;}
  if(n == 2) {
    Real factor=-1.0/(3.0*b);
    T temp=(2.0*f[1]+f[2])*factor;
    u[2]=u[0]=(f[1]+2.0*f[2])*factor;
    u[3]=u[1]=temp;
    return;
  }
	
  int i;
  static int nsize=0;
  static Real *gamma=NULL;
	
  int size=n-1;
  if(size > nsize) {
    nsize=size;
    gamma=new(gamma,nsize) Real;
  }
	
  Real alpha=-1.5*b;
  Real beta=b;
  Real binv=1.0/b;
  Real alphainv=0.5;
  Real delta=-0.5;
	
  u[1]=-0.5*f[1]*binv;

  for(i=2; i <= n-2; i++)	{
    beta *= alphainv;
    gamma[i]=alphainv=1.0/(2.0-alphainv);
    u[i]=(u[i-1]-f[i]*binv)*alphainv;
    delta *= alphainv;
    alpha -= beta*delta;
  }
	
  u[n]=u[0]=0.0;
  u[n-1]=(u[n-2]-f[n-1]*binv)/(2.0-alphainv);
	
  for(i=n-2; i >= 2; i--) u[i] += gamma[i]*u[i+1];
  u[1] += 0.5*u[2];
  u[n+1]=u[1];
}

// Solve the problem Lu=f for u given f, where L is the n x n matrix
//
// [ a1 b1          c1 ]
// [ c2 a2 b2          ]
// [    c3 a3 b3       ]
// [         ...       ]
// [ bn          cn an ]
//
// [gamma] and [delta] are optional work areas of size Real [n-1]
// that may be set to b and c, respectively. 
//
// On exit: u[0]=u[n], u[n+1]=u[1].
//
// Note: u and f need not be distinct.

template<class T>
void tridiagonalp(int n, T *u, const T *f, const Real *c, const Real *a,
		  const Real *b, Real *gamma, Real *delta)
{
  if(n < 1) msg(ERROR,"Invalid matrix size (%d)",n);
  if(n == 1) {u[2]=u[1]=u[0]=f[1]/a[1]; return;}
  if(n == 2) {
    Real factor=1.0/(a[1]*a[2]-c[1]*b[2]);
    T temp=(a[2]*f[1]-c[1]*f[2])*factor;
    u[2]=u[0]=(a[1]*f[2]-b[2]*f[1])*factor;
    u[3]=u[1]=temp;
    return;
  }
	
  int i;
  static int ngamma=0, ndelta=0;
  static Real *gamma0=NULL, *delta0=NULL;
	
  int size=n-1;
  if(gamma == NULL) {Allocate(gamma0,ngamma,size); gamma=gamma0;}
  if(delta == NULL) {Allocate(delta0,ndelta,size); delta=delta0;}
	
  Real ainv=1.0/a[1];
  gamma[1]=b[1]*ainv;
  delta[1]=c[1]*ainv;
  u[1]=f[1]*ainv;
  Real beta=b[n];
  T fn=f[n]-beta*u[1];
  Real alpha=a[n]-beta*delta[1];

  for(i=2; i <= n-2; i++)	{
    Real alphainv=1/(a[i]-c[i]*gamma[i-1]);
    beta *= -gamma[i-1];
    gamma[i]=b[i]*alphainv;
    u[i]=(f[i]-c[i]*u[i-1])*alphainv;
    fn -= beta*u[i];
    delta[i]=-c[i]*delta[i-1]*alphainv;
    alpha -= beta*delta[i];
  }
	
  Real alphainv=1.0/(a[n-1]-c[n-1]*gamma[n-2]);
  u[n-1]=(f[n-1]-c[n-1]*u[n-2])*alphainv;
  beta=c[n]-beta*gamma[n-2];
  Real dnm1=(b[n-1]-c[n-1]*delta[n-2])*alphainv;
  T temp=u[0]=u[n]=(fn-beta*u[n-1])/(alpha-beta*dnm1);
  u[n-1] -= dnm1*temp;
	
  for(i=n-2; i >= 1; i--) u[i] -= gamma[i]*u[i+1]+delta[i]*temp;
  u[n+1]=u[1];
}

template<class T>
inline void tridiagonalp(int n, T *u, const T *f,
			 const Real *c, const Real *a, const Real *b)
{
  tridiagonalp(n,u,f,c,a,b,(Real *) NULL,(Real *) NULL);
}

// Solve the problem Lu=f for u given f, where L is the n x n matrix
//
// [ a1 b1          c1 ]
// [ c2 a2 b2          ]
// [    c3 a3 b3       ]
// [         ...       ]
// [ bn          cn an ]
//
// where a_i+b_i+c_i=0 for i=1,...,n.
//
// For n > 2 this matrix is singular; the solution can only be determined
// to within an arbitrary constant.
//
// [gamma] is an optional work area of size Real [n-1], which may be set to b.
//
// On exit: u[0]=u[n], u[n+1]=u[1].
//
// Note: u and f need not be distinct.

template<class T>
void tridiagonalP(int n, T *u, const T *f, const Real *c, const Real *a,
		  const Real *b, Real *gamma)
{
  if(n < 1) msg(ERROR,"Invalid matrix size (%d)",n);
  if(n == 1) {u[2]=u[1]=u[0]=f[1]/a[1]; return;}
  if(n == 2) {
    Real factor=1.0/(a[1]*a[2]-c[1]*b[2]);
    T temp=(a[2]*f[1]-c[1]*f[2])*factor;
    u[2]=u[0]=(a[1]*f[2]-b[2]*f[1])*factor;
    u[3]=u[1]=temp;
    return;
  }
	
  int i;
  static int ngamma=0;
  static Real *gamma0=NULL;
	
  int size=n-1;
  if(gamma == NULL) {Allocate(gamma0,ngamma,size); gamma=gamma0;}
	
  Real ainv=1.0/a[1];
  gamma[1]=b[1]*ainv;
  Real delta=c[1]*ainv;
  u[1]=f[1]*ainv;
  Real beta=b[n];
  Real alpha=a[n]-beta*delta;

  for(i=2; i <= n-2; i++)	{
    Real alphainv=1/(a[i]-c[i]*gamma[i-1]);
    beta *= -gamma[i-1];
    gamma[i]=b[i]*alphainv;
    u[i]=(f[i]-c[i]*u[i-1])*alphainv;
    delta *= -c[i]*alphainv;
    alpha -= beta*delta;
  }
	
  u[n]=u[0]=0.0;
  u[n-1]=(f[n-1]-c[n-1]*u[n-2])/(a[n-1]-c[n-1]*gamma[n-2]);
	
  for(i=n-2; i >= 1; i--) u[i] -= gamma[i]*u[i+1];
  u[n+1]=u[1];
}

template<class T>
inline void tridiagonalP(int n, T *u, const T *f,
			 const Real *c, const Real *a, const Real *b)
{
  tridiagonalP(n,u,f,c,a,b,(Real *) NULL);
}

// Solve multiple problems Lu=f for u given f, where L is the n x n matrix
//
// [ a1 b1        c1 ]
// [ c2 a2 b2        ]
// [   c3 a3 b3      ]
// [       ...       ]
// [ bn        cn an ]
//
// m is the number of tridiagonal systems to solve,
// inc1 is the stride between the elements of each vector,
// inc2 is the stride between the first elements of the vectors,
//
// [gamma] and [delta] are optional work areas of size Real [(n-2)*inc1+m*inc2]
// that may be set to b and c, respectively.
//
// On entry: the ith element of the jth vector is stored in u[i*inc1+j*inc2],
// where i=0,...,n+1 and j=0,...,m-1.
//
// On exit: u[j*inc2]=u[n*inc1+j*inc2], u[(n+1)*inc1+j*inc2]=u[inc1+j*inc2],
// for each j=0,...,m-1.
//
// Note: u and f need not be distinct.

template<class T>
void mtridiagonalp(int n, T *u, const T *f,
		   const Real *c, const Real *a, const Real *b,
		   int m, int inc1, int inc2, Real *gamma, Real *delta)
{
  int jstop=m*inc2;
		
#ifndef _CRAYMVP
  if(inc1 == 1) {
    for(int j=0; j < jstop; j += inc2)
      tridiagonalp(n,u+j,f+j,c+j,a+j,b+j,gamma,delta);
    return;
  }
#endif	
	
  if(n < 1) msg(ERROR,"Invalid matrix size (%d)",n);

  if(n == 1) {
    T *u1=u+inc1, *u2=u1+inc1;
    const T *f1=f+inc1, *a1=a+inc1;
    for(int j=0; j < jstop; j += inc2) u2[j]=u1[j]=u[j]=f1[j]/a1[j];
    return;
  }
	
  if(n == 2) {
    T *u1=u+inc1, *u2=u1+inc1, *u3=u2+inc1;
    const T *f1=f+inc1, *f2=f1+inc1;
    const T *a1=a+inc1, *a2=a1+inc1;
    const T *c1=c+inc1, *b2=b+2*inc1;
    for(int j=0; j < jstop; j += inc2) {
      Real factor=1.0/(a1[j]*a2[j]-c1[j]*b2[j]);
      T temp=(a2[j]*f1[j]-c1[j]*f2[j])*factor;
      u2[j]=u[j]=(a1[j]*f2[j]-b2[j]*f1[j])*factor;
      u3[j]=u1[j]=temp;
    }
    return;
  }
	
  int i;
  static int msize=0, ngamma=0, ndelta=0;
  static Real *alpha=NULL, *beta=NULL, *gamma0=NULL, *delta0=NULL;
  static T *Fn=NULL;
	
  if(m > msize) {
    msize=m;
    alpha=new(alpha,msize) Real;
    beta=new(beta,msize) Real;
    Fn=new(Fn,msize) T;
  }
	
  int size=(n-2)*inc1+m*inc2;
  if(gamma == NULL) {Allocate(gamma0,ngamma,size); gamma=gamma0;}
  if(delta == NULL) {Allocate(delta0,ndelta,size); delta=delta0;}
	
  int j,ninc1=n*inc1;
  const Real *an=a+ninc1, *bn=b+ninc1, *cn=c+ninc1, *fn=f+ninc1;
	
  const Real *ai=a+inc1, *bi=b+inc1, *ci=c+inc1;
  Real *gammai=gamma+inc1, *deltai=delta+inc1;
  const T *fi=f+inc1;
  T *ui=u+inc1, *pFn=Fn;
  Real *palpha=alpha, *pbeta=beta; 
	
  for(j=0; j < jstop; j += inc2) {
    Real ainv=1.0/ai[j];
    gammai[j]=bi[j]*ainv;
    deltai[j]=ci[j]*ainv;
    ui[j]=fi[j]*ainv;
    Real beta0=*(pbeta++)=bn[j];
    *(pFn++)=fn[j]-beta0*ui[j];
    *(palpha++)=an[j]-beta0*deltai[j];
  }

  for(i=2; i <= n-2; i++) {
    ai += inc1; bi += inc1; ci += inc1;
    Real *gammaim1=gammai, *deltaim1=deltai;
    T *uim1=ui, *pFn=Fn;
    gammai += inc1; deltai += inc1; ui += inc1; fi += inc1;
    Real *palpha=alpha, *pbeta=beta; 
    for(j=0; j < jstop; j += inc2) {
      Real alphainv=1/(ai[j]-ci[j]*gammaim1[j]);
      (*pbeta) *= -gammaim1[j];
      gammai[j]=bi[j]*alphainv;
      ui[j]=(fi[j]-ci[j]*uim1[j])*alphainv;
      *(pFn++) -= (*pbeta)*ui[j];
      deltai[j]=-ci[j]*deltaim1[j]*alphainv;
      *(palpha++) -= (*pbeta)*deltai[j];
      pbeta++;
    }
  }
	
  pFn=Fn;
  palpha=alpha, pbeta=beta; 
  T *unm1=ui+inc1, *un=unm1+inc1;
  const Real *anm1=ai+inc1, *bnm1=bi+inc1, *cnm1=ci+inc1;
  const T *fnm1=fi+inc1;
  for(j=0; j < jstop; j += inc2) {
    Real alphainv=1.0/(anm1[j]-cnm1[j]*gammai[j]);
    unm1[j]=(fnm1[j]-cnm1[j]*ui[j])*alphainv;
    Real beta0=cn[j]-*(pbeta++)*gammai[j];
    Real dnm1=(bnm1[j]-cnm1[j]*deltai[j])*alphainv;
    u[j]=un[j]=((*pFn++)-beta0*unm1[j])/(*(palpha++)-beta0*dnm1);
    unm1[j] -= dnm1*u[j];
  }
	
  T *uip1=unm1;
  for(i=n-2; i >= 1; i--) {
    for(j=0; j < jstop; j += inc2) {
      ui[j] -= gammai[j]*uip1[j]+deltai[j]*u[j];
    }
    gammai -= inc1; deltai -= inc1; uip1=ui; ui -= inc1;
  }
  T *u1=u+inc1, *unp1=un+inc1;
  for(j=0; j < jstop; j += inc2) unp1[j]=u1[j];
}

template<class T>
inline void mtridiagonalp(int n, T *u, const T *f,
			  const Real *c, const Real *a, const Real *b,
			  int m, int inc1, int inc2)
{
  mtridiagonalp(n,u,f,c,a,b,m,inc1,inc2,(Real *) NULL,(Real *) NULL);
}

class CubicSpline() {
protected:
  vector D;
public:
  Calculate(Real t) {
    unsigned int n=C.Size();
    
    static DynVector<V> temp(n);
    if(n > temp.Alloc()) temp.Realloc(n);
    static Array1<Real>::opt work0;
      if(!Active(work0)) work0.Allocate(n-1);
    
      work[0]=-0.5;
      u[0]=f[0]*0.5;
	
      for(i=1; i < n-1; i++) {
	Real temp=1.0/(4.0+work[i-1]);
	work[i]=-temp;
	u[i]=(f[i]-u[i-1])*temp;
      }
      u[n-1]=(f[n-1]-u[n-2])/(2.0+work[n-2]);
      
      for(i=n-2; i >= 0; i--) u[i] += work[i]*u[i+1];

      return ;
  }
  
  Interpolate()
}


}

#endif
