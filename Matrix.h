/* Matrix.h:  A matrix class build upon Array.h
Copyright (C) 1999 John C. Bowman (bowman@math.ualberta.ca)

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

#define __MATRIX_H_VERSION__ 1.10J

#include "DynVector.h"
#include "Array.h"
#include <assert.h>

template<class T>
inline void MatrixSwap(T& p, T& q)
{
	T temp; temp=p; p=q; q=temp;
}

template<class T>
inline void SetIdentity(Array2<T>& A)
{
	A=0.0;
	unsigned int ny1=A.Ny()+1, size=A.Size();
	for(unsigned int i=0; i < size; i += ny1) A(i)=1.0;
}

template<class T>
inline Array2<T> Identity(unsigned int n, unsigned int m, T var=(T) 0.0)
{
	Array2<T> A(n,m);
	SetIdentity(A);
	A.Hold();
	return A;
}

template<class T>
inline void Trans(const Array2<T>& A)
{
	unsigned int nx=A.Nx();
	unsigned int ny=A.Ny();
	assert(nx == ny);
	
	for(unsigned int i=1; i < nx; i++) {
		T *ATi=A()+i;
		T *Ai=A[i];
		for(int j=0; j < i; j++) {
			T temp=Ai[j];
			Ai[j]=ATi[j*ny];
			ATi[j*ny]=temp;
		}
	}
}

template<class T>
inline void Conj(const Array2<T>& A)
{
	unsigned int nx=A.Nx();
	unsigned int ny=A.Ny();
	assert(nx == ny);
	
	A(0)=conj(A(0));
	for(unsigned int i=1; i < nx; i++) {
		T *ATi=A()+i;
		T *Ai=A[i];
		Ai[i]=conj(Ai[i]);
		for(unsigned int j=0; j < i; j++) {
			T temp=conj(Ai[j]);
			Ai[j]=conj(ATi[j*ny]);
			ATi[j*ny]=temp;
		}
	}
}

template<class T>
inline void Trans(const Array2<T>& A, const Array2<T>& B)
{
	unsigned int nx=A.Nx();
	unsigned int ny=A.Ny();
	assert(&A != &B && nx == B.Ny() && ny == B.Nx());
	
	for(unsigned int i=0; i < nx; i++) {
		T *Ai=A[i];
		T *BTi=B()+i;
		for(unsigned int j=0; j < ny; j++) {
			Ai[j]=BTi[j*nx];
		}
	}
	B.Purge();
}

template<class T>
inline void Conj(const Array2<T>& A, const Array2<T>& B)
{
	unsigned int nx=A.Nx();
	unsigned int ny=A.Ny();
	assert(&A != &B && nx == B.Ny() && ny == B.Nx());
	
	for(unsigned int i=0; i < nx; i++) {
		T *Ai=A[i];
		T *BTi=B()+i;
		for(unsigned int j=0; j < ny; j++) {
			Ai[j]=conj(BTi[j*nx]);
		}
	}
	B.Purge();
}

template<class T>
inline void Conj(const Array1<T>& A, const Array1<T>& B)
{
	unsigned int size=A.Size();
	assert(size == B.Size());
	for(unsigned int i=0; i < size; i++) A(i)=conj(B(i));
	B.Purge();
}

template<class T>
inline Array1<T> conj(const Array1<T>& B)
{
	Array1<T> A(B.Nx());
	Conj(A,B);
	A.Hold();
	return A;
}

template<class T>
inline Array2<T> operator ~ (const Array2<T>& B)
{
	Array2<T> A(B.Ny(),B.Nx());
	Trans(A,B);
	A.Hold();
	return A;
}

template<class T>
inline Array1<T> operator * (const Array1<T>& B)
{
	Array1<T> A(B.Nx());
	Conj(A,B);
	A.Hold();
	return A;
}

template<class T>
inline Array2<T> operator * (const Array2<T>& B)
{
	Array2<T> A(B.Ny(),B.Nx());
	Conj(A,B);
	A.Hold();
	return A;
}

template<class T>
inline void real(const Array1<Real>& A, const Array1<T>& B)
{
	unsigned int size=A.Size();
	assert(size == B.Size());
	for(unsigned int i=0; i < size; i++) A(i)=real(B(i));
	B.Purge();
	return;
}

template<class T>
inline Array1<Real> real(const Array1<T>& B)
{
	unsigned int n=B.Nx();
	Array1<Real> A(n);
	real(A,B);
	A.Hold();
	return A;
}

template<class T>
inline void imag(const Array1<Real>& A, const Array1<T>& B)
{
	unsigned int size=A.Size();
	assert(size == B.Size());
	for(unsigned int i=0; i < size; i++) A(i)=imag(B(i));
	B.Purge();
	return;
}

template<class T>
inline Array1<Real> imag(const Array1<T>& B)
{
	unsigned int n=B.Nx();
	Array1<Real> A(n);
	imag(A,B);
	A.Hold();
	return A;
}

template<class T, class U, class V>
inline void Mult(const Array2<T>& A, const Array2<U>& B, const Array2<V>& C)
{
	unsigned int n=C.Size();
	static DynVector<V> temp(n);
	if(n > temp.Alloc()) temp.Resize(n);
	
	unsigned int bny=B.Ny();
	unsigned int cny=C.Ny();
	assert(bny == C.Nx() && cny == A.Ny());
	
	Array2<V> CT(cny,C.Nx(),temp);
	
	Trans(CT,C);
	C.Purge();
	
	unsigned int bnx=B.Nx();
	assert(bnx == A.Nx());
	
	if(&A != &B) {
		for(unsigned int i=0; i < bnx; i++) {
			Array1(T) Ai=A[i];
			Array1(U) Bi=B[i];
			for(unsigned int k=0; k < cny; k++) {
				T sum=0.0;
				Array1(V) CTk=CT[k];
				for(unsigned int j=0; j < bny; j++) {
					sum += Bi[j]*CTk[j];
				}
				Ai[k]=sum;
			}
		}
		B.Purge();
	} else {
		static DynVector<T> temp2(cny);
		if(cny > temp2.Alloc()) temp2.Resize(cny);
		T *work=temp2;
		
		for(unsigned int i=0; i < bnx; i++) {
			Array1(T) Ai=A[i];
			unsigned int k;
			for(k=0; k < cny; k++) {
				T sum=0.0;
				Array1(V) CTk=CT[k];
				for(unsigned int j=0; j < bny; j++) {
					sum += Ai[j]*CTk[j];
				}
				work[k]=sum;
			}
			for(k=0; k < cny; k++) Ai[k]=work[k];
		}
	}		
}

template<class T, class U, class V>
inline void Mult(const Array1<T>& A, const Array1<U>& B, const Array2<V>& C)
{
	unsigned int n=C.Size();
	static DynVector<V> temp(n);
	if(n > temp.Alloc()) temp.Resize(n);
	
	unsigned int bny=B.Nx();
	unsigned int cny=C.Ny();
	assert(bny == C.Nx() && cny == A.Nx());
	
	Array2<V> CT(cny,C.Nx(),temp);
	
	Trans(CT,C);
	C.Purge();
	
	if(&A != &B) {
		for(unsigned int k=0; k < cny; k++) {
			T sum=0.0;
			Array1(V) CTk=CT[k];
			for(unsigned int j=0; j < bny; j++) {
				sum += B[j]*CTk[j];
			}
			A[k]=sum;
		}
		B.Purge();
	} else {
		static DynVector<T> temp2(cny);
		if(cny > temp2.Alloc()) temp2.Resize(cny);
		T *work=temp2;
		
		unsigned int k;
		for(k=0; k < cny; k++) {
			T sum=0.0;
			Array1(V) CTk=CT[k];
			for(unsigned int j=0; j < bny; j++) {
				sum += A[j]*CTk[j];
			}
			work[k]=sum;
		}
		for(k=0; k < cny; k++) A[k]=work[k];
	}		
}

template<class T, class U, class V, class W>
inline void MultAdd(const Array1<T>& A, const Array1<U>& B, const Array2<V>& C,
				 const Array1<W>& D)
{
	unsigned int n=C.Size();
	static DynVector<V> temp(n);
	if(n > temp.Alloc()) temp.Resize(n);
	
	unsigned int bny=B.Nx();
	unsigned int cny=C.Ny();
	assert(bny == C.Nx() && cny == A.Nx());
	
	Array2<V> CT(cny,C.Nx(),temp);
	
	Trans(CT,C);
	C.Purge();
	
	if(&A != &B) {
		for(unsigned int k=0; k < cny; k++) {
			T sum=0.0;
			Array1(V) CTk=CT[k];
			for(unsigned int j=0; j < bny; j++) {
				sum += B[j]*CTk[j];
			}
			A[k]=sum+D[k];
		}
		B.Purge();
	} else {
		static DynVector<T> temp2(cny);
		if(cny > temp2.Alloc()) temp2.Resize(cny);
		T *work=temp2;
		
		unsigned int k;
		for(k=0; k < cny; k++) {
			T sum=0.0;
			Array1(V) CTk=CT[k];
			for(unsigned int j=0; j < bny; j++) {
				sum += A[j]*CTk[j];
			}
			work[k]=sum;
		}
		for(k=0; k < cny; k++) A[k]=work[k]+D[k];
	}		
}

template<class T, class U, class V>
inline void Mult(const Array1<T>& A, const Array2<U>& B, const Array1<V>& C)
{
	unsigned int bny=B.Ny();
	unsigned int bnx=B.Nx();
	assert(bnx == A.Nx() && bny == C.Nx());
	
	for(unsigned int i=0; i < bnx; i++) {
		Array1(U) Bi=B[i];
		T sum=0.0;
		for(unsigned int j=0; j < bny; j++) {
			sum += Bi[j]*C[j];
		}
		A[i]=sum;
	}
	B.Purge();
	C.Purge();
}

template<class T, class U, class V, class W>
inline void MultAdd(const Array1<T>& A, const Array2<U>& B, const Array1<V>& C,
					const Array1<W>& D)
{
	unsigned int bny=B.Ny();
	unsigned int bnx=B.Nx();
	assert(bny == C.Nx() && bnx == A.Nx());
	
	for(unsigned int i=0; i < bnx; i++) {
		Array1(U) Bi=B[i];
		T sum=0.0;
		for(unsigned int j=0; j < bny; j++) {
			sum += Bi[j]*C[j];
		}
		A[i]=sum+D[i];
	}
	B.Purge();
	C.Purge();
}

template<class T, class U, class V, class W>
inline void MultAdd(const Array2<T>& A, const Array2<U>& B, const Array2<V>& C,
					const Array2<W>& D)
{
	unsigned int bny=B.Ny();
	assert(bny == C.Nx());

	unsigned int n=C.Size();
	static DynVector<V> temp(n);
	if(n > temp.Alloc()) temp.Resize(n);
	
	unsigned int cny=C.Ny();
	assert(cny == A.Ny() && cny == D.Ny());
	Array2<V> CT(cny,C.Nx(),temp);
	
	Trans(CT,C);
	C.Purge();
	
	unsigned int bnx=B.Nx();
	assert(bnx == A.Nx() && bnx == D.Nx());
	
	if(&A != &B) {
		for(unsigned int i=0; i < bnx; i++) {
			Array1(T) Ai=A[i];
			Array1(U) Bi=B[i];
			Array1(W) Di=D[i];
			for(unsigned int k=0; k < cny; k++) {
				T sum=0.0;
				Array1(T) CTk=CT[k];
				for(unsigned int j=0; j < bny; j++) {
					sum += Bi[j]*CTk[j];
				}
				Ai[k]=sum+Di[k];
			}
		}
		B.Purge();
	} else {
		static DynVector<T> temp2(cny);
		if(cny > temp2.Alloc()) temp2.Resize(cny);
		T *work=temp2;
		
		for(unsigned int i=0; i < bnx; i++) {
			Array1(T) Ai=A[i];
			unsigned int k;
			for(k=0; k < cny; k++) {
				T sum=0.0;
				Array1(V) CTk=CT[k];
				for(unsigned int j=0; j < bny; j++) {
					sum += Ai[j]*CTk[j];
				}
				work[k]=sum;
			}
			Array1(W) Di=D[i];
			for(k=0; k < cny; k++) Ai[k]=work[k]+Di[k];
		}
	}		
}

template<class T, class U>
inline Array2<T> const operator * (const Array2<T>& B, const Array2<U>& C)
{
	Array2<T> A(B.Nx(),C.Ny());
	Mult(A,B,C);
	A.Hold();
	return A;
}

template<class T, class U>
inline Array1<T> const operator * (const Array2<T>& B, const Array1<U>& C)
{
	Array1<T> A(B.Nx());
	Mult(A,B,C);
	A.Hold();
	return A;
}

template<class T, class U>
inline Array1<T> const operator * (const Array1<T>& B, const Array2<U>& C)
{
	Array1<T> A(B.Nx());
	Mult(A,B,C);
	A.Hold();
	return A;
}

template<class T, class U, class V>
inline void Mult(const array1<T>& A, const array1<U>& B, V C)
{
	unsigned int size=A.Size(); 
	for(unsigned int i=0; i < size; i++) A(i)=B(i)*C;
	B.Purge();
}

template<class T, class U>
inline Array2<T> operator * (const Array2<T>& B, U C)
{
	Array2<T> A(B.Nx(),B.Ny());
	Mult(A,B,C);
	A.Hold();
	return A;
}

template<class T, class U>
inline Array1<T> operator * (const Array1<T>& B, U C)
{
	Array1<T> A(B.Nx());
	Mult(A,B,C);
	A.Hold();
	return A;
}

template<class T, class U>
inline Array2<T> operator * (U C, const Array2<T>& B)
{
	Array2<T> A(B.Nx(),B.Ny());
	Mult(A,B,C);
	A.Hold();
	return A;
}

template<class T, class U>
inline Array1<T> operator * (U C, const Array1<T>& B)
{
	Array1<T> A(B.Nx());
	Mult(A,B,C);
	A.Hold();
	return A;
}

template<class T>
inline Array2<T>& Array2<T>::operator *= (const Array2<T>& A)
{
	Mult(*this,*this,A);
	return *this;
}

#if 0 // Defined in Array.h
template<class T>
inline Array2<T>& Array2<T>::operator *= (T A)
{
	Mult(*this,*this,A);
	return *this;
}
#endif

template<class T>
void GaussJordan(const Array2<T>& a, const Array2<T>& b)
{
// Linear equation solution by Gauss-Jordan elimination.
// a[1..n][1..n] is the input matrix.
// b[1..n][1..m] is input containing the m right-hand side vectors.
// On output, a is replaced by its matrix inverse, and b is replaced by  
// the corresponding set of solution vectors.

// The integer arrays pivot, indx_r, and indx_c are used for 
// bookkeeping on the pivoting.
	unsigned int n=a.Nx();
	unsigned int m=b.Ny();
	assert(&a != &b && n == a.Ny() && n == b.Nx());
	
	unsigned int indx_c[n], indx_r[n], pivot[n];

	for(unsigned int j=0; j < n; j++) {
		pivot[j]=0;
	}
	unsigned int col=0, row=0;
	for(unsigned int i=0; i < n; i++) {
		// This is the main loop over the columns to be reduced.
		double big=0.0;
		for(unsigned int j=0; j < n; j++) {
			Array1(T) aj=a[j];
			// This is the outer loop of the search for a pivot element.
			if(pivot[j] != 1) {
				for(unsigned int k=0; k < n; k++) {
					if(pivot[k] == 0) {
						double temp=abs2(aj[k]);
						if(temp >= big) {
							big=temp;
							row=j;
							col=k;
						}
					} else if(pivot[k] > 1) 
					{__ARRAY_EXIT("Singular matrix");}
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
		Array1(T) acol=a[col];
		Array1(T) bcol=b[col];
		if(row != col) {
			Array1(T) arow=a[row];
			Array1(T) brow=b[row];
			for(unsigned int l=0; l < n; l++) MatrixSwap(arow[l], acol[l]);
			for(unsigned int l=0; l < m; l++) MatrixSwap(brow[l], bcol[l]);
		}
		indx_r[i]=row; 
		// We are now ready to divide the pivot row by the pivot element,
		// located at row and col. 
		indx_c[i]=col;
		if(acol[col] == 0.0) {__ARRAY_EXIT("Singular matrix");}
		T pivinv=1.0/acol[col];
		acol[col]=1.0;
		for(unsigned int l=0; l < n; l++) acol[l] *= pivinv;
		for(unsigned int l=0; l < m; l++) bcol[l] *= pivinv;
		for(unsigned int ll=0; ll < n; ll++) {
			// Next, we reduce the rows...
			if(ll != col) { 
				//...except for the pivot one, of course.
				Array1(T) all=a[ll];
				Array1(T) bll=b[ll];
				T dum=all[col];
				all[col]=0.0;
				for(unsigned int l=0; l < n; l++) all[l] -= acol[l]*dum;
				for(unsigned int l=0; l < m; l++) bll[l] -= bcol[l]*dum;
			}
		}
	}
	
	// This is the end of the main loop over columns of the reduction. 
    // It only remains to unscramble the solution in view of the column
	// interchanges. We do this by interchanging pairs of
	// columns in the reverse order that the permutation was built up.
	for(int l=n-1; l >= 0; l--) {
		if(indx_r[l] != indx_c[l]) {
			for(unsigned int k=0; k < n; k++) 
				MatrixSwap(a[k][indx_r[l]], a[k][indx_c[l]]);
		}
	}
}
	
// Compute A=C^{-1} B

template<class T>
inline void Divide(const Array2<T>& A, const Array2<T>& C, const Array2<T>& B)
{
	if(&A != &B) {
		assert(A.Nx() == B.Nx() && A.Ny() == B.Ny());
		A.Load(B);
		B.Purge();
	}
	
	unsigned int n=C.Size();
	static DynVector<T> temp(n);
	if(n > temp.Alloc()) temp.Resize(n);
	
	Array2<T> D(C.Nx(),C.Ny(),temp);
	D.Load(C);
	C.Purge();
	GaussJordan(D,A);
};
	
template<class T>
inline void Divide(const Array1<T>& A, const Array2<T>& C, const Array1<T>& B)
{
	if(&A != &B) {
		assert(A.Nx() == B.Nx());
		A.Load(B);
		B.Purge();
	}
	
	unsigned int n=C.Size();
	static DynVector<T> temp(n);
	if(n > temp.Alloc()) temp.Resize(n);
	
	Array2<T> D(C.Nx(),C.Ny(),temp);
	Array2<T> A2(A.Nx(),1,A);
	
	D.Load(C);
	C.Purge();
	GaussJordan(D,A2);
};
	
template<class T>
inline void Divide(const Array2<T>& A, T C, const Array2<T>& B)
{
	T Cinv=1.0/C;
	unsigned int size=A.Size(); 
	for(unsigned int i=0; i < size; i++) A(i)=B(i)*Cinv;
	B.Purge();
}

template<class T>
inline Array2<T> operator / (const Array2<T>& B, T C)
{
	Array2<T> A(B.Nx(),B.Ny());
	Divide(A,C,B);
	A.Hold();
	return A;
}

template<class T>
inline void Add(const array1<T>& A, const array1<T>& B, const array1<T>& C)
{
	unsigned int size=A.Size();
	assert(size == B.Size() && size == C.Size());
	for(unsigned int i=0; i < size; i++) A(i)=B(i)+C(i);
	B.Purge();
	C.Purge();
}

template<class T>
inline void Add(const Array2<T>& A, const Array2<T>& B, const Array2<T>& C)
{
	assert(A.Nx() == B.Nx() && A.Ny() == B.Ny() && 
		   B.Nx() == C.Nx() && B.Ny() == C.Ny());
	
	unsigned int size=A.Size();
	for(unsigned int i=0; i < size; i++) A(i)=B(i)+C(i);
	B.Purge();
	C.Purge();
}

template<class T>
inline void Sub(const array1<T>& A, const array1<T>& B, const array1<T>& C)
{
	unsigned int size=A.Size();
	assert(size == B.Size() && size == C.Size());
	for(unsigned int i=0; i < size; i++) A(i)=B(i)-C(i);
	B.Purge();
	C.Purge();
}

template<class T>
inline void Sub(const Array2<T>& A, const Array2<T>& B, const Array2<T>& C)
{
	assert(A.Nx() == B.Nx() && A.Ny() == B.Ny() && 
		   B.Nx() == C.Nx() && B.Ny() == C.Ny());
	
	unsigned int size=A.Size();
	for(unsigned int i=0; i < size; i++) A(i)=B(i)-C(i);
	B.Purge();
	C.Purge();
}

template<class T>
inline void Unary(const array1<T>& A, const array1<T>& B)
{
	unsigned int size=A.Size();
	assert(size == B.Size());
	for(unsigned int i=0; i < size; i++) A(i)=-B(i);
	B.Purge();
}

template<class T>
inline void Unary(const Array2<T>& A, const Array2<T>& B)
{
	assert(A.Nx() == B.Nx() && A.Ny() == B.Ny());
	unsigned int size=A.Size();
	for(unsigned int i=0; i < size; i++) A(i)=-B(i);
	B.Purge();
}

template<class T>
inline void Add(array1<T>& A, const array1<T>& B, T C)
{
	unsigned int size=A.Size();
	assert(size == B.Size());
	for(unsigned int i=0; i < size; i++) A(i)=B(i)+C;
	B.Purge();
}

template<class T>
inline void Add(Array2<T>& A, const Array2<T>& B, T C)
{
	assert(A.Nx() == B.Nx() && A.Ny() == B.Ny());
	unsigned int size=A.Size();
	for(unsigned int i=0; i < size; i++) A(i)=B(i)+C;
	B.Purge();
}

template<class T>
inline void Sub(array1<T>& A, const array1<T>& B, T C)
{
	unsigned int size=A.Size();
	assert(size == B.Size());
	for(unsigned int i=0; i < size; i++) A(i)=B(i)-C;
	B.Purge();
}

template<class T>
inline void Sub(Array2<T>& A, const Array2<T>& B, T C)
{
	assert(A.Nx() == B.Nx() && A.Ny() == B.Ny());
	unsigned int size=A.Size();
	for(unsigned int i=0; i < size; i++) A(i)=B(i)-C;
	B.Purge();
}

template<class T>
inline void Add(const array1<T>& A, T B, const array1<T>& C)
{
	unsigned int size=A.Size();
	assert(size == C.Size());
	for(unsigned int i=0; i < size; i++) A(i)=B+C(i);
	C.Purge();
}

template<class T>
inline void Add(const Array2<T>& A, T B, const Array2<T>& C)
{
	assert(A.Nx() == C.Nx() && A.Ny() == C.Ny());
	unsigned int size=A.Size();
	for(unsigned int i=0; i < size; i++) A(i)=B+C(i);
	C.Purge();
}

template<class T>
inline void Sub(const array1<T>& A, T B, const array1<T>& C)
{
	unsigned int size=A.Size();
	assert(size == C.Size());
	for(unsigned int i=0; i < size; i++) A(i)=B-C(i);
	C.Purge();
}

template<class T>
inline void Sub(const Array2<T>& A, T B, const Array2<T>& C)
{
	assert(A.Nx() == C.Nx() && A.Ny() == C.Ny());
	unsigned int size=A.Size();
	for(unsigned int i=0; i < size; i++) A(i)=B-C(i);
	C.Purge();
}

template<class T>
inline Array1<T> operator + (const Array1<T>& B, const Array1<T>& C)
{
	Array1<T> A(B.Nx());
	Add(A,B,C);
	A.Hold();
	return A;
}

template<class T>
inline Array2<T> operator + (const Array2<T>& B, const Array2<T>& C)
{
	Array2<T> A(B.Nx(),B.Ny());
	Add(A,B,C);
	A.Hold();
	return A;
}

template<class T>
inline Array1<T> operator - (const Array1<T>& B, const Array1<T>& C)
{
	Array1<T> A(B.Nx());
	Sub(A,B,C);
	A.Hold();
	return A;
}

template<class T>
inline Array2<T> operator - (const Array2<T>& B, const Array2<T>& C)
{
	Array2<T> A(B.Nx(),B.Ny());
	Sub(A,B,C);
	A.Hold();
	return A;
}

template<class T>
inline Array2<T> operator + (const Array2<T>& B, T C)
{
	Array2<T> A(B.Nx(),B.Ny());
	Add(A,B,C);
	A.Hold();
	return A;
}

template<class T>
inline Array2<T> operator - (const Array2<T>& B, T C)
{
	Array2<T> A(B.Nx(),B.Ny());
	Sub(A,B,C);
	A.Hold();
	return A;
}

template<class T>
inline Array2<T> operator + (T B, const Array2<T>& C)
{
	Array2<T> A(C.Nx(),C.Ny());
	Add(A,B,C);
	A.Hold();
	return A;
}

template<class T>
inline Array2<T> operator - (T B, const Array2<T>& C)
{
	Array2<T> A(C.Nx(),C.Ny());
	Sub(A,B,C);
	A.Hold();
	return A;
}

template<class T>
inline Array2<T> operator - (const Array2<T>& B)
{
	Array2<T> A(B.Nx(),B.Ny());
	Unary(A,B);
	A.Hold();
	return A;
}

#define log2 __log2
#define pow2 __pow2

// Compute the matrix exponential B of a square matrix A
template<class T, class U>
inline void Exp(Array2<T>& B, const Array2<U>& A)
{
	unsigned int n=A.Nx();
    assert(n == A.Ny());
	unsigned int n2=n*n;
	
	double Amax=0.0;
	for(unsigned int j=0; j < n2; j++) Amax=max(Amax,abs2(A(j)));
	Amax=sqrt(Amax);
	
	int j=max(0,1+((int) (floor(log2(Amax))+0.5)));
	double scale=pow2(-j);
	B=A*scale;
	A.Purge();
	
	unsigned int q=1;
//	double delta=DBL_EPSILON;
	double delta=1.0e-11;
	double epsilon=8.0;
	while (epsilon > delta) {
		epsilon /= 8*(4*q*q-1);
		q++;
	}
	q--;
	
	// Improve to reuse work arrays here: JCB
	
	Array2<T> D(n,n);
	Array2<T> N(n,n);
	Array2<T> X(n,n);
	
	D.Identity();
	X=N=D;
	
	double c=1.0;
	double sign=-1.0;
	
	for(unsigned int k=0; k < q; k++) {
		X=B*X;
		c *= (q-k)/((double) ((q+q-k)*(k+1)));
		N += c*X;
		D += (sign*c)*X;
		sign = -sign;
	}
	
	Divide(B,D,N); // B=D^{-1} N
	
	for(int k=0; k < j; k++) B *= B;
}	

template<class T>
const Array2<T> exp(const Array2<T>& A)
{
	unsigned int n=A.Nx();
	Array2<T> B(n,n);
	Exp(B,A);
	B.Hold();
	return B;
}

#endif
