#ifndef __Tridiagonal_h__
#define __Tridiagonal_h__ 1

#include <new>
#include "Array.h"

namespace Array {
  
template<class T, class S>
inline void CheckReallocate(T& A, S& B, unsigned int n, unsigned int& old)
{
  if(n > old) {Reallocate(A,n); Reallocate(B,n); old=n;}
}


template<class T, class S, class U>
inline void CheckReallocate(T& A, S& B, U& C,
			    unsigned int n, unsigned int& old)
{
  if(n > old) {Reallocate(A,n); Reallocate(B,n); Reallocate(C,n); old=n;}
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
// work is an optional work array of size n+1, which may be set to a or b.
//
// Note: u and f need not be distinct.

template<class T, class C>
inline void tridiagonal(unsigned int n,
			const typename array1<T>::opt& u,
			const typename array1<T>::opt& f,
			const typename array1<C>::opt& c,
			const typename array1<C>::opt& a,
			const typename array1<C>::opt& b,
			const typename array1<C>::opt& work)
{
  if(n < 1) msg(ERROR,"Invalid matrix size (%d)",n);
	
  C temp=1.0/a[1];
  work[1]=-b[1]*temp;
  u[1]=(f[1]-c[1]*u[0])*temp;
	
  for(unsigned int i=2; i <= n; i++) {
    C temp=1.0/(a[i]+c[i]*work[i-1]);
    work[i]=-b[i]*temp;
    u[i]=(f[i]-c[i]*u[i-1])*temp;
  }

  for(unsigned int i=n; i >= 1; i--) u[i] += work[i]*u[i+1];
}

template<class T, class C>
inline void tridiagonal(unsigned int n,
			const typename array1<T>::opt& u,
			const typename array1<T>::opt& f,
			const typename array1<C>::opt& c,
			const typename array1<C>::opt& a,
			const typename array1<C>::opt& b)
{
  static typename array1<C>::opt work;
  static unsigned int size=0;
  CheckReallocate(work,n+1,size);
  tridiagonal<T,C>(n,u,f,c,a,b,work);
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
// work is an optional work array of size n*inc1+m*inc2,
// which may be set to a or b.
//
// Note: u and f need not be distinct.

template<class T, class C>
inline void mtridiagonal(unsigned int n,
			 const typename array1<T>::opt& u,
			 const typename array1<T>::opt& f,
			 const typename array1<C>::opt& c,
			 const typename array1<C>::opt& a,
			 const typename array1<C>::opt& b,
			 unsigned int m, unsigned int inc1, unsigned int inc2,
			 const typename array1<C>::opt& work)
{
  unsigned int jstop=m*inc2;
		
#ifndef _CRAYMVP
  if(inc1 == 1) {
    for(unsigned int j=0; j < jstop; j += inc2) 
      tridiagonal<T,C>(n,u+j,f+j,c+j,a+j,b+j,work);
    return;
  }
#endif	
	
  if(n < 1) msg(ERROR,"Invalid matrix size (%d)",n);

  typename array1<C>::opt ai=a+inc1, bi=b+inc1, ci=c+inc1, worki=work+inc1;
  typename array1<T>::opt fi=f+inc1, ui=u+inc1;
  for(unsigned int j=0; j < jstop; j += inc2) {
    C temp=1.0/ai[j];
    worki[j]=-bi[j]*temp;
    ui[j]=(fi[j]-ci[j]*u[j])*temp;
  }
	
  for(unsigned int i=2; i <= n; i++) {
    ai += inc1; bi += inc1; ci += inc1;
    typename array1<C>::opt workim1=worki;
    typename array1<T>::opt uim1=ui;
    worki += inc1; ui += inc1; fi += inc1;
    for(unsigned int j=0; j < jstop; j += inc2) {
      C temp=1.0/(ai[j]+ci[j]*workim1[j]);
      worki[j]=-bi[j]*temp;
      ui[j]=(fi[j]-ci[j]*uim1[j])*temp;
    }
  }

  typename array1<T>::opt uip1=ui+inc1;
  for(unsigned int i=n; i >= 1; i--) {
    for(unsigned int j=0; j < jstop; j += inc2) {
      ui[j] += worki[j]*uip1[j];
    }
    worki -= inc1; uip1=ui; ui -= inc1;
  }
}

template<class T, class C>
inline void mtridiagonal(unsigned int n,
			 const typename array1<T>::opt& u,
			 const typename array1<T>::opt& f,
			 const typename array1<C>::opt& c,
			 const typename array1<C>::opt& a,
			 const typename array1<C>::opt& b,
			 unsigned int m, unsigned int inc1, unsigned int inc2)
{
  static typename array1<C>::opt work;
  static unsigned int size=0;
  CheckReallocate(work,n*inc1+m*inc2,size);
  mtridiagonal<T,C>(n,u,f,c,a,b,m,inc1,inc2,work);
}

// Solve the problem Lu=f for u given f, where L is the n x n matrix
//
// [ a1 b1          c1 ]
// [ c2 a2 b2          ]
// [    c3 a3 b3       ]
// [         ...       ]
// [ bn          cn an ]
//
// gamma and delta are optional work areas of size n-1, which may be set to 
// b and c, respectively.
//
// On exit: u[0]=u[n], u[n+1]=u[1].
//
// Note: u and f need not be distinct.

template<class T, class C>
inline void tridiagonalp(unsigned int n,
			 const typename array1<T>::opt& u,
			 const typename array1<T>::opt& f,
			 const typename array1<C>::opt& c,
			 const typename array1<C>::opt& a,
			 const typename array1<C>::opt& b,
			 const typename array1<C>::opt& gamma,
			 const typename array1<C>::opt& delta)
{
  if(n < 1) msg(ERROR,"Invalid matrix size (%d)",n);
  if(n == 1) {u[2]=u[1]=u[0]=f[1]/a[1]; return;}
  if(n == 2) {
    C factor=1.0/(a[1]*a[2]-c[1]*b[2]);
    T temp=(a[2]*f[1]-c[1]*f[2])*factor;
    u[2]=u[0]=(a[1]*f[2]-b[2]*f[1])*factor;
    u[3]=u[1]=temp;
    return;
  }
	
  C ainv=1.0/a[1];
  gamma[1]=b[1]*ainv;
  delta[1]=c[1]*ainv;
  u[1]=f[1]*ainv;
  C beta=b[n];
  T fn=f[n]-beta*u[1];
  C alpha=a[n]-beta*delta[1];

  for(unsigned int i=2; i <= n-2; i++)	{
    C alphainv=1/(a[i]-c[i]*gamma[i-1]);
    beta *= -gamma[i-1];
    gamma[i]=b[i]*alphainv;
    u[i]=(f[i]-c[i]*u[i-1])*alphainv;
    fn -= beta*u[i];
    delta[i]=-c[i]*delta[i-1]*alphainv;
    alpha -= beta*delta[i];
  }
	
  C alphainv=1.0/(a[n-1]-c[n-1]*gamma[n-2]);
  u[n-1]=(f[n-1]-c[n-1]*u[n-2])*alphainv;
  beta=c[n]-beta*gamma[n-2];
  C dnm1=(b[n-1]-c[n-1]*delta[n-2])*alphainv;
  T temp=u[0]=u[n]=(fn-beta*u[n-1])/(alpha-beta*dnm1);
  u[n-1] -= dnm1*temp;
	
  for(unsigned int i=n-2; i >= 1; i--) u[i] -= gamma[i]*u[i+1]+delta[i]*temp;
  u[n+1]=u[1];
}

template<class T, class C>
inline void tridiagonalp(unsigned int n,
			 const typename array1<T>::opt& u,
			 const typename array1<T>::opt& f,
			 const typename array1<C>::opt& c,
			 const typename array1<C>::opt& a,
			 const typename array1<C>::opt& b)
{
  if(n < 1) msg(ERROR,"Invalid matrix size (%d)",n);
  static typename array1<C>::opt gamma,delta;
  static unsigned int size=0;
  CheckReallocate(gamma,delta,n-1,size);
  tridiagonalp<T,C>(n,u,f,c,a,b,gamma,delta);
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
// gamma is an optional work area of size n-1, which may be set to b.
//
// On exit: u[0]=u[n], u[n+1]=u[1].
//
// Note: u and f need not be distinct.

template<class T, class C>
inline void tridiagonalP(unsigned int n,
			 const typename array1<T>::opt& u,
			 const typename array1<T>::opt& f,
			 const typename array1<C>::opt& c,
			 const typename array1<C>::opt& a,
			 const typename array1<C>::opt& b,
			 const typename array1<C>::opt& gamma)
{
  if(n < 1) msg(ERROR,"Invalid matrix size (%d)",n);
  if(n == 1) {u[2]=u[1]=u[0]=f[1]/a[1]; return;}
  if(n == 2) {
    C factor=1.0/(a[1]*a[2]-c[1]*b[2]);
    T temp=(a[2]*f[1]-c[1]*f[2])*factor;
    u[2]=u[0]=(a[1]*f[2]-b[2]*f[1])*factor;
    u[3]=u[1]=temp;
    return;
  }
	
  C ainv=1.0/a[1];
  gamma[1]=b[1]*ainv;
  C delta=c[1]*ainv;
  u[1]=f[1]*ainv;
  C beta=b[n];
  C alpha=a[n]-beta*delta;

  for(unsigned int i=2; i <= n-2; i++)	{
    C alphainv=1/(a[i]-c[i]*gamma[i-1]);
    beta *= -gamma[i-1];
    gamma[i]=b[i]*alphainv;
    u[i]=(f[i]-c[i]*u[i-1])*alphainv;
    delta *= -c[i]*alphainv;
    alpha -= beta*delta;
  }
	
  u[n]=u[0]=0.0;
  u[n-1]=(f[n-1]-c[n-1]*u[n-2])/(a[n-1]-c[n-1]*gamma[n-2]);
	
  for(unsigned int i=n-2; i >= 1; i--) u[i] -= gamma[i]*u[i+1];
  u[n+1]=u[1];
}

template<class T, class C>
inline void tridiagonalP(unsigned int n,
			 const typename array1<T>::opt& u,
			 const typename array1<T>::opt& f,
			 const typename array1<C>::opt& c,
			 const typename array1<C>::opt& a,
			 const typename array1<C>::opt& b)
{
  if(n < 1) msg(ERROR,"Invalid matrix size (%d)",n);
  static typename array1<C>::opt gamma;
  static unsigned int size=0;
  CheckReallocate(gamma,n-1,size);
  tridiagonalP<T,C>(n,u,f,c,a,b,gamma);
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
// Fn is an optional work array of size m.
// alpha and beta are optional work array of size m.
// gamma and delta are optional work arrays of size (n-2)*inc1+m*inc2,
// which may be set to b and c, respectively.
//
// On entry: the ith element of the jth vector is stored in u[i*inc1+j*inc2],
// where i=0,...,n+1 and j=0,...,m-1.
//
// On exit: u[j*inc2]=u[n*inc1+j*inc2], u[(n+1)*inc1+j*inc2]=u[inc1+j*inc2],
// for each j=0,...,m-1.
//
// Note: u and f need not be distinct.

template<class T, class C>
inline void mtridiagonalp(unsigned int n,
			  const typename array1<T>::opt& u,
			  const typename array1<T>::opt& f,
			  const typename array1<C>::opt& c,
			  const typename array1<C>::opt& a,
			  const typename array1<C>::opt& b,
			  unsigned int m, unsigned int inc1, unsigned int inc2,
			  const typename array1<T>::opt& Fn,
			  const typename array1<C>::opt& alpha,
			  const typename array1<C>::opt& beta,
			  const typename array1<C>::opt& gamma,
			  const typename array1<C>::opt& delta)
{
  unsigned int jstop=m*inc2;
		
#ifndef _CRAYMVP
  if(inc1 == 1) {
    for(unsigned int j=0; j < jstop; j += inc2)
      tridiagonalp<T,C>(n,u+j,f+j,c+j,a+j,b+j,gamma,delta);
    return;
  }
#endif	
	
  if(n < 1) msg(ERROR,"Invalid matrix size (%d)",n);

  if(n == 1) {
    typename array1<C>::opt a1=a+inc1;
    typename array1<T>::opt u1=u+inc1, u2=u1+inc1, f1=f+inc1;
    for(unsigned int j=0; j < jstop; j += inc2) u2[j]=u1[j]=u[j]=f1[j]/a1[j];
    return;
  }
	
  if(n == 2) {
    typename array1<T>::opt u1=u+inc1, u2=u1+inc1, u3=u2+inc1;
    typename array1<T>::opt f1=f+inc1, f2=f1+inc1;
    typename array1<C>::opt a1=a+inc1, a2=a1+inc1, c1=c+inc1, b2=b+2*inc1;
    for(unsigned int j=0; j < jstop; j += inc2) {
      C factor=1.0/(a1[j]*a2[j]-c1[j]*b2[j]);
      T temp=(a2[j]*f1[j]-c1[j]*f2[j])*factor;
      u2[j]=u[j]=(a1[j]*f2[j]-b2[j]*f1[j])*factor;
      u3[j]=u1[j]=temp;
    }
    return;
  }
	
  unsigned int ninc1=n*inc1;
  typename array1<C>::opt an=a+ninc1, bn=b+ninc1, cn=c+ninc1, fn=f+ninc1;
  typename array1<C>::opt ai=a+inc1, bi=b+inc1, ci=c+inc1;
  typename array1<C>::opt gammai=gamma+inc1, deltai=delta+inc1;
  typename array1<T>::opt ui=u+inc1, fi=f+inc1;
	
  for(unsigned int j=0, k=0; j < jstop; j += inc2, k++) {
    C ainv=1.0/ai[j];
    gammai[j]=bi[j]*ainv;
    deltai[j]=ci[j]*ainv;
    ui[j]=fi[j]*ainv;
    C beta0=beta[k]=bn[j];
    Fn[k]=fn[j]-beta0*ui[j];
    alpha[k]=an[j]-beta0*deltai[j];
  }

  for(unsigned int i=2; i <= n-2; i++) {
    ai += inc1; bi += inc1; ci += inc1;
    typename array1<C>::opt gammaim1=gammai, deltaim1=deltai;
    typename array1<T>::opt uim1=ui;
    gammai += inc1; deltai += inc1; ui += inc1; fi += inc1;
    for(unsigned int j=0, k=0; j < jstop; j += inc2, k++) {
      C alphainv=1/(ai[j]-ci[j]*gammaim1[j]);
      beta[k] *= -gammaim1[j];
      gammai[j]=bi[j]*alphainv;
      ui[j]=(fi[j]-ci[j]*uim1[j])*alphainv;
      Fn[k] -= beta[k]*ui[j];
      deltai[j]=-ci[j]*deltaim1[j]*alphainv;
      alpha[k] -= beta[k]*deltai[j];
    }
  }
	
  typename array1<T>::opt unm1=ui+inc1, un=unm1+inc1, fnm1=fi+inc1;
  typename array1<C>::opt anm1=ai+inc1, bnm1=bi+inc1, cnm1=ci+inc1;
  for(unsigned int j=0, k=0; j < jstop; j += inc2, k++) {
    C alphainv=1.0/(anm1[j]-cnm1[j]*gammai[j]);
    unm1[j]=(fnm1[j]-cnm1[j]*ui[j])*alphainv;
    C beta0=cn[j]-beta[k]*gammai[j];
    C dnm1=(bnm1[j]-cnm1[j]*deltai[j])*alphainv;
    u[j]=un[j]=(Fn[k]-beta0*unm1[j])/(alpha[k]-beta0*dnm1);
    unm1[j] -= dnm1*u[j];
  }
	
  typename array1<T>::opt uip1=unm1;
  for(unsigned int i=n-2; i >= 1; i--) {
    for(unsigned int j=0; j < jstop; j += inc2) {
      ui[j] -= gammai[j]*uip1[j]+deltai[j]*u[j];
    }
    gammai -= inc1; deltai -= inc1; uip1=ui; ui -= inc1;
  }
  typename array1<T>::opt u1=u+inc1, unp1=un+inc1;
  for(unsigned int j=0; j < jstop; j += inc2) unp1[j]=u1[j];
}

template<class T, class C>
inline void mtridiagonalp(unsigned int n,
			  const typename array1<T>::opt& u,
			  const typename array1<T>::opt& f,
			  const typename array1<C>::opt& c,
			  const typename array1<C>::opt& a,
			  const typename array1<C>::opt& b,
			  unsigned int m, unsigned int inc1, unsigned int inc2)
{
  if(n < 1) msg(ERROR,"Invalid matrix size (%d)",n);
  static typename array1<T>::opt Fn;
  static typename array1<C>::opt alpha,beta;
  static unsigned int msize=0;
  static typename array1<C>::opt gamma,delta;
  static unsigned int size=0;
  if(n > 2) {
    CheckReallocate(Fn,alpha,beta,m,msize);
    CheckReallocate(gamma,delta,(n-2)*inc1+m*inc2,size);
  }
  mtridiagonalp<T,C>(n,u,f,c,a,b,m,inc1,inc2,Fn,alpha,beta,gamma,delta);
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

template<class T, class C>
inline void Poisson1(unsigned int n,
		     const typename array1<T>::opt& u,
		     const typename array1<T>::opt& f,
		     const C& b, const C& a,
		     const typename array1<C>::opt& work)
{
  if(n < 1) msg(ERROR,"Invalid matrix size (%d)",n);
	
  C temp=b/a;
  work[1]=-temp;
  C binv=1.0/b;
  u[1]=(f[1]*binv-u[0])*temp;
	
  for(unsigned int i=2; i <= n; i++)	{
    C temp=b/(a+b*work[i-1]);
    work[i]=-temp;
    u[i]=(f[i]*binv-u[i-1])*temp;
  }

  for(unsigned int i=n; i >= 1; i--) u[i] += work[i]*u[i+1];
}

template<class T, class C>
inline void Poisson1(unsigned int n,
		     const typename array1<T>::opt& u,
		     const typename array1<T>::opt& f,
		     const C& b, const C& a)
{
  static typename array1<C>::opt work;
  static unsigned int size=0;
  CheckReallocate(work,n+1,size);
  Poisson1<T,C>(n,u,f,b,a,work);
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
// gamma and delta are optional work arrays of size n-1.
//
// On exit: u[0]=u[n], u[n+1]=u[1].
//
// Note: u and f need not be distinct.

template<class T, class C>
inline void Poisson1p(unsigned int n,
		      const typename array1<T>::opt& u,
		      const typename array1<T>::opt& f,
		      const C& b, const C& a,
		      const typename array1<C>::opt& gamma,
		      const typename array1<C>::opt& delta)
{
  if(n < 1) msg(ERROR,"Invalid matrix size (%d)",n);
  if(n == 1) {u[2]=u[1]=u[0]=f[1]/a; return;}
  if(n == 2) {
    C factor=1.0/(a*a-b*b);
    T temp=(a*f[1]-b*f[2])*factor;
    u[2]=u[0]=(-b*f[1]+a*f[2])*factor;
    u[3]=u[1]=temp;
    return;
  }
	
  C binv=1.0/b;
  C ainv=b/a;
  delta[1]=gamma[1]=ainv;
  C alpha=a-b*ainv;
  u[1]=f[1]*binv*ainv;
  T fn=f[n]-b*u[1];
  C beta=b;

  for(unsigned int i=2; i <= n-2; i++) {
    C alphainv=b/(a-b*gamma[i-1]);
    beta *= -gamma[i-1];
    gamma[i]=alphainv;
    u[i]=(f[i]*binv-u[i-1])*alphainv;
    fn -= beta*u[i];
    delta[i]=-delta[i-1]*alphainv;
    alpha -= beta*delta[i];
  }
	
  C alphainv=b/(a-b*gamma[n-2]);
  u[n-1]=(f[n-1]*binv-u[n-2])*alphainv;
  beta=b-beta*gamma[n-2];
  C dnm1=(1.0-delta[n-2])*alphainv;
  T temp=u[0]=u[n]=(fn-beta*u[n-1])/(alpha-beta*dnm1);
  u[n-1] -= dnm1*temp;
	
  for(unsigned int i=n-2; i >= 1; i--) 
    u[i] -= gamma[i]*u[i+1]+delta[i]*temp;

  u[n+1]=u[1];
}

template<class T, class C>
inline void Poisson1p(unsigned int n,
		      const typename array1<T>::opt& u,
		      const typename array1<T>::opt& f,
		      const C& b, const C& a)
{
  if(n < 1) msg(ERROR,"Invalid matrix size (%d)",n);
  static typename array1<C>::opt gamma,delta;
  static unsigned int size=0;
  CheckReallocate(gamma,delta,n-1,size);
  Poisson1p<T,C>(n,u,f,b,a,gamma,delta);
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
// gamma is an optional work array of size n-1.
//
// On exit: u[0]=u[n], u[n+1]=u[1].
//
// For n > 2 this matrix is singular; the solution can only be determined
// to within an arbitrary constant.
//
// Note: u and f need not be distinct.

template<class T, class C>
inline void Poisson1P(unsigned int n,
		      const typename array1<T>::opt& u,
		      const typename array1<T>::opt& f,
		      C b, const typename array1<C>::opt& gamma)
{
  if(n < 1) msg(ERROR,"Invalid matrix size (%d)",n);
  if(n == 1) {u[2]=u[1]=u[0]=-f[1]/(2.0*b); return;}
  if(n == 2) {
    C factor=-1.0/(3.0*b);
    T temp=(2.0*f[1]+f[2])*factor;
    u[2]=u[0]=(f[1]+2.0*f[2])*factor;
    u[3]=u[1]=temp;
    return;
  }
	
  C alpha=-1.5*b;
  C beta=b;
  C binv=1.0/b;
  C alphainv=0.5;
  C delta=-0.5;
	
  u[1]=-0.5*f[1]*binv;

  for(unsigned int i=2; i <= n-2; i++)	{
    beta *= alphainv;
    gamma[i]=alphainv=1.0/(2.0-alphainv);
    u[i]=(u[i-1]-f[i]*binv)*alphainv;
    delta *= alphainv;
    alpha -= beta*delta;
  }
	
  u[n]=u[0]=0.0;
  u[n-1]=(u[n-2]-f[n-1]*binv)/(2.0-alphainv);
	
  for(unsigned int i=n-2; i >= 2; i--) u[i] += gamma[i]*u[i+1];
  u[1] += 0.5*u[2];
  u[n+1]=u[1];
}

template<class T, class C>
inline void Poisson1P(unsigned int n,
		      const typename array1<T>::opt& u,
		      const typename array1<T>::opt& f, const C& b)
{
  if(n < 1) msg(ERROR,"Invalid matrix size (%d)",n);
  static typename array1<C>::opt gamma;
  static unsigned int size=0;
  CheckReallocate(gamma,n-1,size);
  Poisson1P<T,C>(n,u,f,b,gamma);
}

// Solve the problem Lu=f for u given f, where L is the n x n matrix
//
// [ a0 b0            ]
// [ c1 a1 b1         ]
// [    c2 a2 b2      ]
// [           ...    ]
// [            cm am ]
//
// where m=n-1. Note: u and f need not be distinct.
// work is an optional work area of size n-1, which may be set to a or b.

template<class T, class C>
inline void Tridiagonal(unsigned int n, 
			const typename array1<T>::opt& u,
			const typename array1<T>::opt& f,
			const typename array1<C>::opt& c,
			const typename array1<C>::opt& a,
			const typename array1<C>::opt& b,
			const typename array1<C>::opt& work)
{
  if(n < 1) exit(-1);
//  if(n < 1) msg(ERROR,"Invalid matrix size (%d)",n);
    
  C temp=1.0/a[0];
  u[0]=f[0]*temp;
  
  if(n == 1) return;
  
  work[0]=-b[0]*temp;
	
  for(unsigned int i=1; i < n-1; i++) {
    C temp=1.0/(a[i]+c[i]*work[i-1]);
    work[i]=-b[i]*temp;
    u[i]=(f[i]-c[i]*u[i-1])*temp;
  }
  
  temp=1.0/(a[n-1]+c[n-1]*work[n-2]);
  u[n-1]=(f[n-1]-c[n-1]*u[n-2])*temp;

  for(int i=(int)n-2; i >= 0; i--) u[i] += work[i]*u[i+1];
}

template<class T, class C>
inline void Tridiagonal(unsigned int n, 
			const typename array1<T>::opt& u,
			const typename array1<T>::opt& f,
			const typename array1<C>::opt& c,
			const typename array1<C>::opt& a,
			const typename array1<C>::opt& b)
{
  if(n < 1) exit(-1);
  static typename array1<C>::opt work;
  static unsigned int worksize=0;
  CheckReallocate(work,n-1,worksize);
  Tridiagonal<T,C>(n,u,f,c,a,b,work);
}

// Do a binary search of an ordered array to find an interval containing key.
// Return the index corresponding to the left-hand endpoint of the matching
// interval, or -1 if key is less than the first element.
template<class X>
unsigned int bintsearch(X key, unsigned int n, const typename array1<X>::opt x)
{
  if(key < x[0]) return -1;
  if(key >= x[n-1]) return n-1;
  
  unsigned int l=0;
  unsigned int u=n-1;
	
  while (l < u) {
    unsigned int i=(l + u)/2;
    if(x[i] <= key && key < x[i+1]) return i;
    if(key < x[i]) u=i;
    else l=i + 1;
  }
//  msg(ERROR,"Statement not reachable");
  exit(1);
}

template<class Y, class X>
class CubicSpline {
protected:
  static typename array1<X>::opt a,b,c;
  static typename array1<Y>::opt y2,f;
  static unsigned int size;
public:
  CubicSpline() {size=0;}
  CubicSpline(unsigned int n, const typename array1<X>::opt x,
	      const typename array1<Y>::opt y) {
    if(n > size) {
      Reallocate(y2,n);
      Reallocate(a,n);
      Reallocate(b,n);
      Reallocate(c,n);
      Reallocate(f,n);
      size=n;
    }
    for(unsigned int i=1; i < n-1; i++) {
      c[i]=(x[i]-x[i-1])/6.0;
      a[i]=(x[i+1]-x[i-1])/3.0;
      b[i]=(x[i+1]-x[i])/6.0;
      f[i]=(y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]);
    }
      
    y2[0]=y2[n-1]=0.0;
    if(n > 2) Tridiagonal<Y,X>(n-2,y2+1,f+1,c+1,a+1,b+1);
    cout << y2 << endl;
    return;
  }
  
  Y Interpolate(unsigned int n, const typename array1<X>::opt x,
	      const typename array1<Y>::opt y, X x0) {
    return bintsearch(x0,n,x);
  }
};

template<class Y, class X>
unsigned int CubicSpline<Y,X>::size=0;

template<class Y, class X>
typename array1<X>::opt CubicSpline<Y,X>::a;

template<class Y, class X>
typename array1<X>::opt CubicSpline<Y,X>::b;

template<class Y, class X>
typename array1<X>::opt CubicSpline<Y,X>::c;

template<class Y, class X>
typename array1<Y>::opt CubicSpline<Y,X>::y2;

template<class Y, class X>
typename array1<Y>::opt CubicSpline<Y,X>::f;


}
#endif
