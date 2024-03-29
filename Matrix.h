/* Matrix.h:  A matrix class build upon Array.h
Copyright (C) 1999-2004 John C. Bowman (bowman@math.ualberta.ca)

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. */

#ifndef __Matrix_h__
#define __Matrix_h__ 1

#define __MATRIX_H_VERSION__ 1.14

#include "Array.h"
#include <cassert>
#include <cmath>

namespace Array {
  
template<class T>
inline void MatrixSwap(T& p, T& q)
{
  T temp; temp=p; p=q; q=temp;
}

template<class T>
inline array2<T> Identity(size_t n, size_t m)
{
  array2<T> A(n,m);
  A.Identity();
  A.Hold();
  return A;
}

template<class T>
inline void Trans(const array2<T>& A)
{
  size_t nx=A.Nx();
  size_t ny=A.Ny();
  assert(nx == ny);
	
  for(size_t i=1; i < nx; i++) {
    typename array1<T>::opt ATi=A()+i;
    typename array1<T>::opt Ai=A[i];
    for(int j=0; j < i; j++) {
      T temp=Ai[j];
      Ai[j]=ATi[j*ny];
      ATi[j*ny]=temp;
    }
  }
}

template<class T>
inline void Conj(const array2<T>& A)
{
  size_t nx=A.Nx();
  size_t ny=A.Ny();
  assert(nx == ny);
	
  A(0)=conj(A(0));
  for(size_t i=1; i < nx; i++) {
    typename array1<T>::opt ATi=A()+i;
    typename array1<T>::opt Ai=A[i];
    Ai[i]=conj(Ai[i]);
    for(size_t j=0; j < i; j++) {
      T temp=conj(Ai[j]);
      Ai[j]=conj(ATi[j*ny]);
      ATi[j*ny]=temp;
    }
  }
}

template<class T>
inline void Trans(const array2<T>& A, const array2<T>& B)
{
  size_t nx=A.Nx();
  size_t ny=A.Ny();
  assert(&A != &B && nx == B.Ny() && ny == B.Nx());
	
  for(size_t i=0; i < nx; i++) {
    typename array1<T>::opt Ai=A[i];
    typename array1<T>::opt BTi=B()+i;
    for(size_t j=0; j < ny; j++) {
      Ai[j]=BTi[j*nx];
    }
  }
  B.Purge();
}

template<class T>
inline void Conj(const array2<T>& A, const array2<T>& B)
{
  size_t nx=A.Nx();
  size_t ny=A.Ny();
  assert(&A != &B && nx == B.Ny() && ny == B.Nx());
	
  for(size_t i=0; i < nx; i++) {
    typename array1<T>::opt Ai=A[i];
    typename array1<T>::opt BTi=B()+i;
    for(size_t j=0; j < ny; j++) {
      Ai[j]=conj(BTi[j*nx]);
    }
  }
  B.Purge();
}

template<class T>
inline void Conj(const array1<T>& A, const array1<T>& B)
{
  size_t size=A.Size();
  assert(size == B.Size());
  for(size_t i=0; i < size; i++) A(i)=conj(B(i));
  B.Purge();
}

template<class T>
inline array1<T> conj(const array1<T>& B)
{
  array1<T> A(B.Nx());
  Conj(A,B);
  A.Hold();
  return A;
}

template<class T>
inline array2<T> conj(const array2<T>& B)
{
  array1<T> A(B.Nx(),B.Ny());
  Conj(A,B);
  A.Hold();
  return A;
}

template<class T>
inline array2<T> trans(const array2<T>& B)
{
  array2<T> A(B.Nx(),B.Ny());
  Trans(A,B);
  A.Hold();
  return A;
}

template<class T>
inline array2<T> operator ~ (const array2<T>& B)
{
  array2<T> A(B.Ny(),B.Nx());
  Trans(A,B);
  A.Hold();
  return A;
}

template<class T>
inline array1<T> operator * (const array1<T>& B)
{
  array1<T> A(B.Nx());
  Conj(A,B);
  A.Hold();
  return A;
}

template<class T>
inline array2<T> operator * (const array2<T>& B)
{
  array2<T> A(B.Ny(),B.Nx());
  Conj(A,B);
  A.Hold();
  return A;
}

template<class T>
inline void real(const array1<double>& A, const array1<T>& B)
{
  size_t size=A.Size();
  assert(size == B.Size());
  for(size_t i=0; i < size; i++) A(i)=real(B(i));
  B.Purge();
  return;
}

template<class T>
inline array1<double> real(const array1<T>& B)
{
  size_t n=B.Nx();
  array1<double> A(n);
  real(A,B);
  A.Hold();
  return A;
}

template<class T>
inline void imag(const array1<double>& A, const array1<T>& B)
{
  size_t size=A.Size();
  assert(size == B.Size());
  for(size_t i=0; i < size; i++) A(i)=imag(B(i));
  B.Purge();
  return;
}

template<class T>
inline array1<double> imag(const array1<T>& B)
{
  size_t n=B.Nx();
  array1<double> A(n);
  imag(A,B);
  A.Hold();
  return A;
}

template<class T, class U, class V>
inline void Mult(const array2<T>& A, const array2<U>& B, const array2<V>& C)
{
  size_t n=C.Size();
  
  static typename array1<V>::opt temp;
  static size_t tempsize=0;
  CheckReallocate(temp,n,tempsize);
	
  size_t bny=B.Ny();
  size_t cny=C.Ny();
  assert(bny == C.Nx() && cny == A.Ny());
	
  array2<V> CT(cny,C.Nx(),temp);
	
  Trans(CT,C);
  C.Purge();
	
  size_t bnx=B.Nx();
  assert(bnx == A.Nx());
	
  if(&A != &B) {
    for(size_t i=0; i < bnx; i++) {
      typename array1<T>::opt Ai=A[i];
      typename array1<U>::opt Bi=B[i];
      for(size_t k=0; k < cny; k++) {
	T sum=0.0;
	typename array1<V>::opt CTk=CT[k];
	for(size_t j=0; j < bny; j++) {
	  sum += Bi[j]*CTk[j];
	}
	Ai[k]=sum;
      }
    }
    B.Purge();
  } else {
    static typename array1<T>::opt work;
    static size_t worksize=0;
    CheckReallocate(work,cny,worksize);
		
    for(size_t i=0; i < bnx; i++) {
      typename array1<T>::opt Ai=A[i];
      size_t k;
      for(k=0; k < cny; k++) {
	T sum=0.0;
	typename array1<V>::opt CTk=CT[k];
	for(size_t j=0; j < bny; j++) {
	  sum += Ai[j]*CTk[j];
	}
	work[k]=sum;
      }
      for(k=0; k < cny; k++) Ai[k]=work[k];
    }
  }		
}

template<class T, class U, class V>
inline void Mult(const array1<T>& A, const array1<U>& B, const array2<V>& C)
{
  size_t n=C.Size();
  
  static typename array1<V>::opt temp;
  static size_t tempsize=0;
  CheckReallocate(temp,n,tempsize);
  
  size_t bny=B.Nx();
  size_t cny=C.Ny();
  assert(bny == C.Nx() && cny == A.Nx());
	
  array2<V> CT(cny,C.Nx(),temp);
	
  Trans(CT,C);
  C.Purge();
	
  if(&A != &B) {
    for(size_t k=0; k < cny; k++) {
      T sum=0.0;
      typename array1<V>::opt CTk=CT[k];
      for(size_t j=0; j < bny; j++) {
	sum += B[j]*CTk[j];
      }
      A[k]=sum;
    }
    B.Purge();
  } else {
    static typename array1<T>::opt work;
    static size_t worksize=0;
    CheckReallocate(work,cny,worksize);
    
    size_t k;
    for(k=0; k < cny; k++) {
      T sum=0.0;
      typename array1<V>::opt CTk=CT[k];
      for(size_t j=0; j < bny; j++) {
	sum += A[j]*CTk[j];
      }
      work[k]=sum;
    }
    for(k=0; k < cny; k++) A[k]=work[k];
  }		
}

template<class T, class U, class V, class W>
inline void MultAdd(const array1<T>& A, const array1<U>& B, const array2<V>& C,
		    const array1<W>& D)
{
  size_t n=C.Size();
  
  static typename array1<V>::opt temp;
  static size_t tempsize=0;
  CheckReallocate(temp,n,tempsize);
  
  size_t bny=B.Nx();
  size_t cny=C.Ny();
  assert(bny == C.Nx() && cny == A.Nx());
	
  array2<V> CT(cny,C.Nx(),temp);
	
  Trans(CT,C);
  C.Purge();
	
  if(&A != &B) {
    for(size_t k=0; k < cny; k++) {
      T sum=0.0;
      typename array1<V>::opt CTk=CT[k];
      for(size_t j=0; j < bny; j++) {
	sum += B[j]*CTk[j];
      }
      A[k]=sum+D[k];
    }
    B.Purge();
  } else {
    static typename array1<T>::opt work;
    static size_t worksize=0;
    CheckReallocate(work,cny,worksize);
		
    size_t k;
    for(k=0; k < cny; k++) {
      T sum=0.0;
      typename array1<V>::opt CTk=CT[k];
      for(size_t j=0; j < bny; j++) {
	sum += A[j]*CTk[j];
      }
      work[k]=sum;
    }
    for(k=0; k < cny; k++) A[k]=work[k]+D[k];
  }		
}

template<class T, class U, class V>
inline void Mult(const array1<T>& A, const array2<U>& B, const array1<V>& C)
{
  size_t bny=B.Ny();
  size_t bnx=B.Nx();
  assert(bnx == A.Nx() && bny == C.Nx());
	
  for(size_t i=0; i < bnx; i++) {
    typename array1<U>::opt Bi=B[i];
    T sum=0.0;
    for(size_t j=0; j < bny; j++) {
      sum += Bi[j]*C[j];
    }
    A[i]=sum;
  }
  B.Purge();
  C.Purge();
}

template<class T, class U, class V, class W>
inline void MultAdd(const array1<T>& A, const array2<U>& B, const array1<V>& C,
		    const array1<W>& D)
{
  size_t bny=B.Ny();
  size_t bnx=B.Nx();
  assert(bny == C.Nx() && bnx == A.Nx());
	
  for(size_t i=0; i < bnx; i++) {
    typename array1<U>::opt Bi=B[i];
    T sum=0.0;
    for(size_t j=0; j < bny; j++) {
      sum += Bi[j]*C[j];
    }
    A[i]=sum+D[i];
  }
  B.Purge();
  C.Purge();
}

template<class T, class U, class V, class W>
inline void MultAdd(const array2<T>& A, const array2<U>& B, const array2<V>& C,
		    const array2<W>& D)
{
  size_t bny=B.Ny();
  assert(bny == C.Nx());
  size_t n=C.Size();
  
  static typename array1<V>::opt temp;
  static size_t tempsize=0;
  CheckReallocate(temp,n,tempsize);
	
  size_t cny=C.Ny();
  assert(cny == A.Ny() && cny == D.Ny());
  array2<V> CT(cny,C.Nx(),temp);
	
  Trans(CT,C);
  C.Purge();
	
  size_t bnx=B.Nx();
  assert(bnx == A.Nx() && bnx == D.Nx());
	
  if(&A != &B) {
    for(size_t i=0; i < bnx; i++) {
      typename array1<T>::opt Ai=A[i];
      typename array1<U>::opt Bi=B[i];
      typename array1<W>::opt Di=D[i];
      for(size_t k=0; k < cny; k++) {
	T sum=0.0;
	typename array1<T>::opt CTk=CT[k];
	for(size_t j=0; j < bny; j++) {
	  sum += Bi[j]*CTk[j];
	}
	Ai[k]=sum+Di[k];
      }
    }
    B.Purge();
  } else {
    static typename array1<T>::opt work;
    static size_t worksize=0;
    CheckReallocate(work,cny,worksize);
		
    for(size_t i=0; i < bnx; i++) {
      typename array1<T>::opt Ai=A[i];
      size_t k;
      for(k=0; k < cny; k++) {
	T sum=0.0;
	typename array1<V>::opt CTk=CT[k];
	for(size_t j=0; j < bny; j++) {
	  sum += Ai[j]*CTk[j];
	}
	work[k]=sum;
      }
      typename array1<W>::opt Di=D[i];
      for(k=0; k < cny; k++) Ai[k]=work[k]+Di[k];
    }
  }		
}

template<class T, class U>
inline array2<T> const operator * (const array2<T>& B, const array2<U>& C)
{
  array2<T> A(B.Nx(),C.Ny());
  Mult(A,B,C);
  A.Hold();
  return A;
}

template<class T, class U>
inline array1<T> const operator * (const array2<T>& B, const array1<U>& C)
{
  array1<T> A(B.Nx());
  Mult(A,B,C);
  A.Hold();
  return A;
}

template<class T, class U>
inline array1<T> const operator * (const array1<T>& B, const array2<U>& C)
{
  array1<T> A(B.Nx());
  Mult(A,B,C);
  A.Hold();
  return A;
}

template<class T, class U, class V>
inline void Mult(const array1<T>& A, const array1<U>& B, V C)
{
  size_t size=A.Size(); 
  for(size_t i=0; i < size; i++) A(i)=B(i)*C;
  B.Purge();
}

template<class T, class U>
inline array2<T> operator * (const array2<T>& B, U C)
{
  array2<T> A(B.Nx(),B.Ny());
  Mult(A,B,C);
  A.Hold();
  return A;
}

template<class T, class U>
inline array1<T> operator * (const array1<T>& B, U C)
{
  array1<T> A(B.Nx());
  Mult(A,B,C);
  A.Hold();
  return A;
}

template<class T, class U>
inline array2<T> operator * (U C, const array2<T>& B)
{
  array2<T> A(B.Nx(),B.Ny());
  Mult(A,B,C);
  A.Hold();
  return A;
}

template<class T, class U>
inline array1<T> operator * (U C, const array1<T>& B)
{
  array1<T> A(B.Nx());
  Mult(A,B,C);
  A.Hold();
  return A;
}

template<class T>
inline array2<T>& array2<T>::operator *= (const array2<T>& A)
{
  Mult(*this,*this,A);
  return *this;
}

template<class T>
void GaussJordan(const array2<T>& a, const array2<T>& b)
{
  // Linear equation solution by Gauss-Jordan elimination.
  // On input a is an n x n matrix
  // and b is an n x m matrix containing the m right-hand side vectors.
  // On output, a is replaced by its matrix inverse, and b is replaced by  
  // the corresponding set of solution vectors.

  // The integer arrays pivot, indx_r, and indx_c are used for 
  // bookkeeping on the pivoting.
  size_t n=a.Nx();
  size_t m=b.Ny();
  assert(&a != &b && n == a.Ny() && n == b.Nx());
	
  static typename array1<int>::opt indx_c;
  static typename array1<int>::opt indx_r;
  static typename array1<int>::opt pivot;
  static size_t tempsize=0;
  if(n > tempsize) {
    Reallocate(indx_c,n);
    Reallocate(indx_r,n);
    Reallocate(pivot,n);
    tempsize=n;
  }

  for(size_t j=0; j < n; j++) {
    pivot[j]=0;
  }
  size_t col=0, row=0;
  for(size_t i=0; i < n; i++) {
    // This is the main loop over the columns to be reduced.
    double big=0.0;
    for(size_t j=0; j < n; j++) {
      typename array1<T>::opt aj=a[j];
      // This is the outer loop of the search for a pivot element.
      if(pivot[j] != 1) {
	for(size_t k=0; k < n; k++) {
	  if(pivot[k] == 0) {
	    double temp=abs2(aj[k]);
	    if(temp >= big) {
	      big=temp;
	      row=j;
	      col=k;
	    }
	  } else if(pivot[k] > 1) ArrayExit("Singular matrix");
	}
      }
    }
    ++(pivot[col]);
    // We now have the pivot element, so we interchange rows, if
    // needed, to put  the pivot element on the diagonal. The columns
    // are not physically interchanged, only relabeled: 
    // indx_c[i], the column of the ith pivot element, is the ith
    // column that is reduced, while indx_r[i] is the row in which that
    // pivot element was  originally located. If indx_r[i]=indx_c[i]
    // there is an implied column interchange. With this form of
    // bookkeeping, the solution b's will end up in the correct order,
    // and the inverse matrix will be scrambled by columns.
    typename array1<T>::opt acol=a[col];
    typename array1<T>::opt bcol=b[col];
    if(row != col) {
      typename array1<T>::opt arow=a[row];
      typename array1<T>::opt brow=b[row];
      for(size_t l=0; l < n; l++) MatrixSwap(arow[l], acol[l]);
      for(size_t l=0; l < m; l++) MatrixSwap(brow[l], bcol[l]);
    }
    indx_r[i]=row; 
    // We are now ready to divide the pivot row by the pivot element,
    // located at row and col. 
    indx_c[i]=col;
    if(acol[col] == 0.0) ArrayExit("Singular matrix");
    T pivinv=1.0/acol[col];
    acol[col]=1.0;
    for(size_t l=0; l < n; l++) acol[l] *= pivinv;
    for(size_t l=0; l < m; l++) bcol[l] *= pivinv;
    for(size_t ll=0; ll < n; ll++) {
      // Next, we reduce the rows...
      if(ll != col) { 
				//...except for the pivot one, of course.
	typename array1<T>::opt all=a[ll], bll=b[ll];
	T dum=all[col];
	all[col]=0.0;
	for(size_t l=0; l < n; l++) all[l] -= acol[l]*dum;
	for(size_t l=0; l < m; l++) bll[l] -= bcol[l]*dum;
      }
    }
  }
	
  // This is the end of the main loop over columns of the reduction. 
  // It only remains to unscramble the solution in view of the column
  // interchanges. We do this by interchanging pairs of
  // columns in the reverse order that the permutation was built up.
  for(int l=n-1; l >= 0; l--) {
    if(indx_r[l] != indx_c[l]) {
      for(size_t k=0; k < n; k++) 
	MatrixSwap(a[k][indx_r[l]], a[k][indx_c[l]]);
    }
  }
}
	
// Compute A=C^{-1} B

template<class T>
inline void Divide(const array2<T>& A, const array2<T>& C, const array2<T>& B)
{
  if(&A != &B) {
    assert(A.Nx() == B.Nx() && A.Ny() == B.Ny());
    A.Load(B);
    B.Purge();
  }
	
  size_t n=C.Size();
  
  static typename array1<T>::opt temp;
  static size_t tempsize=0;
  CheckReallocate(temp,n,tempsize);
	
  array2<T> D(C.Nx(),C.Ny(),temp);
  D.Load(C);
  C.Purge();
  GaussJordan(D,A);
}
	
// Compute A=C^{-1} B
  
template<class T>
inline void Divide(const array1<T>& A, const array2<T>& C, const array1<T>& B)
{
  if(&A != &B) {
    assert(A.Nx() == B.Nx());
    A.Load(B);
    B.Purge();
  }
	
  size_t n=C.Size();
  
  static typename array1<T>::opt temp;
  static size_t tempsize=0;
  CheckReallocate(temp,n,tempsize);
	
  array2<T> D(C.Nx(),C.Ny(),temp);
  array2<T> A2(A.Nx(),1,A);
	
  D.Load(C);
  C.Purge();
  GaussJordan(D,A2);
}
	
template<class T>
inline void Divide(const array2<T>& A, T C, const array2<T>& B)
{
  T Cinv=1.0/C;
  size_t size=A.Size(); 
  for(size_t i=0; i < size; i++) A(i)=B(i)*Cinv;
  B.Purge();
}

template<class T>
inline array2<T> operator / (const array2<T>& B, T C)
{
  array2<T> A(B.Nx(),B.Ny());
  Divide(A,C,B);
  A.Hold();
  return A;
}

template<class T>
inline void Add(const array1<T>& A, const array1<T>& B, const array1<T>& C)
{
  size_t size=A.Size();
  assert(size == B.Size() && size == C.Size());
  for(size_t i=0; i < size; i++) A(i)=B(i)+C(i);
  B.Purge();
  C.Purge();
}

template<class T>
inline void Add(const array2<T>& A, const array2<T>& B, const array2<T>& C)
{
  assert(A.Nx() == B.Nx() && A.Ny() == B.Ny() && 
	 B.Nx() == C.Nx() && B.Ny() == C.Ny());
	
  size_t size=A.Size();
  for(size_t i=0; i < size; i++) A(i)=B(i)+C(i);
  B.Purge();
  C.Purge();
}

template<class T>
inline void Sub(const array1<T>& A, const array1<T>& B, const array1<T>& C)
{
  size_t size=A.Size();
  assert(size == B.Size() && size == C.Size());
  for(size_t i=0; i < size; i++) A(i)=B(i)-C(i);
  B.Purge();
  C.Purge();
}

template<class T>
inline void Sub(const array2<T>& A, const array2<T>& B, const array2<T>& C)
{
  assert(A.Nx() == B.Nx() && A.Ny() == B.Ny() && 
	 B.Nx() == C.Nx() && B.Ny() == C.Ny());
	
  size_t size=A.Size();
  for(size_t i=0; i < size; i++) A(i)=B(i)-C(i);
  B.Purge();
  C.Purge();
}

template<class T>
inline void Unary(const array1<T>& A, const array1<T>& B)
{
  size_t size=A.Size();
  assert(size == B.Size());
  for(size_t i=0; i < size; i++) A(i)=-B(i);
  B.Purge();
}

template<class T>
inline void Unary(const array2<T>& A, const array2<T>& B)
{
  assert(A.Nx() == B.Nx() && A.Ny() == B.Ny());
  size_t size=A.Size();
  for(size_t i=0; i < size; i++) A(i)=-B(i);
  B.Purge();
}

template<class T>
inline void Add(array1<T>& A, const array1<T>& B, T C)
{
  size_t size=A.Size();
  assert(size == B.Size());
  for(size_t i=0; i < size; i++) A(i)=B(i)+C;
  B.Purge();
}

template<class T>
inline void Add(array2<T>& A, const array2<T>& B, T C)
{
  assert(A.Nx() == B.Nx() && A.Ny() == B.Ny());
  size_t size=A.Size();
  for(size_t i=0; i < size; i++) A(i)=B(i)+C;
  B.Purge();
}

template<class T>
inline void Sub(array1<T>& A, const array1<T>& B, T C)
{
  size_t size=A.Size();
  assert(size == B.Size());
  for(size_t i=0; i < size; i++) A(i)=B(i)-C;
  B.Purge();
}

template<class T>
inline void Sub(array2<T>& A, const array2<T>& B, T C)
{
  assert(A.Nx() == B.Nx() && A.Ny() == B.Ny());
  size_t size=A.Size();
  for(size_t i=0; i < size; i++) A(i)=B(i)-C;
  B.Purge();
}

template<class T>
inline void Add(const array1<T>& A, T B, const array1<T>& C)
{
  size_t size=A.Size();
  assert(size == C.Size());
  for(size_t i=0; i < size; i++) A(i)=B+C(i);
  C.Purge();
}

template<class T>
inline void Add(const array2<T>& A, T B, const array2<T>& C)
{
  assert(A.Nx() == C.Nx() && A.Ny() == C.Ny());
  size_t size=A.Size();
  for(size_t i=0; i < size; i++) A(i)=B+C(i);
  C.Purge();
}

template<class T>
inline void Sub(const array1<T>& A, T B, const array1<T>& C)
{
  size_t size=A.Size();
  assert(size == C.Size());
  for(size_t i=0; i < size; i++) A(i)=B-C(i);
  C.Purge();
}

template<class T>
inline void Sub(const array2<T>& A, T B, const array2<T>& C)
{
  assert(A.Nx() == C.Nx() && A.Ny() == C.Ny());
  size_t size=A.Size();
  for(size_t i=0; i < size; i++) A(i)=B-C(i);
  C.Purge();
}

template<class T>
inline array1<T> operator + (const array1<T>& B, const array1<T>& C)
{
  array1<T> A(B.Nx());
  Add(A,B,C);
  A.Hold();
  return A;
}

template<class T>
inline array2<T> operator + (const array2<T>& B, const array2<T>& C)
{
  array2<T> A(B.Nx(),B.Ny());
  Add(A,B,C);
  A.Hold();
  return A;
}

template<class T>
inline array1<T> operator - (const array1<T>& B, const array1<T>& C)
{
  array1<T> A(B.Nx());
  Sub(A,B,C);
  A.Hold();
  return A;
}

template<class T>
inline array2<T> operator - (const array2<T>& B, const array2<T>& C)
{
  array2<T> A(B.Nx(),B.Ny());
  Sub(A,B,C);
  A.Hold();
  return A;
}

template<class T>
inline array2<T> operator + (const array2<T>& B, T C)
{
  array2<T> A(B.Nx(),B.Ny());
  Add(A,B,C);
  A.Hold();
  return A;
}

template<class T>
inline array2<T> operator - (const array2<T>& B, T C)
{
  array2<T> A(B.Nx(),B.Ny());
  Sub(A,B,C);
  A.Hold();
  return A;
}

template<class T>
inline array2<T> operator + (T B, const array2<T>& C)
{
  array2<T> A(C.Nx(),C.Ny());
  Add(A,B,C);
  A.Hold();
  return A;
}

template<class T>
inline array2<T> operator - (T B, const array2<T>& C)
{
  array2<T> A(C.Nx(),C.Ny());
  Sub(A,B,C);
  A.Hold();
  return A;
}

template<class T>
inline array2<T> operator - (const array2<T>& B)
{
  array2<T> A(B.Nx(),B.Ny());
  Unary(A,B);
  A.Hold();
  return A;
}

// Compute the matrix exponential B of a square matrix A
template<class T, class U>
inline void Exp(array2<T>& B, const array2<U>& A)
{
  size_t n=A.Nx();
  assert(n == A.Ny());
  size_t n2=n*n;
	
  double Amax=0.0;
  for(size_t j=0; j < n2; j++) {
    double value=abs2(A(j));
    if(value > Amax) Amax=value;
  }
	
  size_t j=1+((int) (floor(0.5*log2(Amax))+0.5));
  if(j < 0) j=0;
  double scale=1.0/(1 << j);
  B=A*scale;
  A.Purge();
	
  size_t q=1;
  double delta=DBL_EPSILON;
  double epsilon=8.0;
  while (epsilon > delta) {
    epsilon /= 8*(4*q*q-1);
    q++;
  }
  q--;
	
  static typename array1<T>::opt Dtemp;
  static typename array1<T>::opt Ntemp;
  static typename array1<T>::opt Xtemp;
  static size_t tempsize=0;
  if(n2 > tempsize) {
    Reallocate(Dtemp,n2);
    Reallocate(Ntemp,n2);
    Reallocate(Xtemp,n2); 
    tempsize=n2;
  }
  
  array2<T> D(n,n,Dtemp);
  array2<T> N(n,n,Ntemp);
  array2<T> X(n,n,Xtemp);
	
  D.Identity();
  X=N=D;
	
  double c=1.0;
  double sign=-1.0;
	
  for(size_t k=0; k < q; k++) {
    X=B*X;
    c *= (q-k)/((double) ((q+q-k)*(k+1)));
    N += c*X;
    D += (sign*c)*X;
    sign = -sign;
  }
	
  Divide(B,D,N); // B=D^{-1} N
	
  for(size_t k=0; k < j; k++) B *= B;
}	

template<class T>
const array2<T> exp(const array2<T>& A)
{
  size_t n=A.Nx();
  array2<T> B(n,n);
  Exp(B,A);
  B.Hold();
  return B;
}

}

#endif
