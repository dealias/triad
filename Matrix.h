#ifndef __Matrix_h__
#define __Matrix_h__ 1

#include "DynVector.h"
#include "Array.h"
#include <assert.h>

template<class T>
inline void SetIdentity(Array2<T>& A)
{
	A=0.0;
	int ny1=A.Ny()+1, size=A.Size();
	for(int i=0; i < size; i += ny1) A(i)=1.0;
}

template<class T>
inline Array2<T> Identity(int n, int m, T)
{
	Array2<T> A(n,m);
	SetIdentity(A);
	A.Freeze();
	return A;
}

template<class T>
inline void Trans(const Array2<T>& A)
{
	int nx=A.Nx();
	int ny=A.Ny();
	assert(nx == ny);
	
	for(int i=1; i < nx; i++) {
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
	int nx=A.Nx();
	int ny=A.Ny();
	assert(nx == ny);
	
	A(0)=conj(A(0));
	for(int i=1; i < nx; i++) {
		T *ATi=A()+i;
		T *Ai=A[i];
		Ai[i]=conj(Ai[i]);
		for(int j=0; j < i; j++) {
			T temp=conj(Ai[j]);
			Ai[j]=conj(ATi[j*ny]);
			ATi[j*ny]=temp;
		}
	}
}

template<class T>
inline void Trans(const Array2<T>& A, const Array2<T>& B)
{
	int nx=A.Nx();
	int ny=A.Ny();
	assert(&A != &B && nx == B.Ny() && ny == B.Nx());
	
	for(int i=0; i < nx; i++) {
		T *Ai=A[i];
		T *BTi=B()+i;
		for(int j=0; j < ny; j++) {
			Ai[j]=BTi[j*nx];
		}
	}
	B.Purge();
}

template<class T>
inline void Conjugate(const Array2<T>& A, const Array2<T>& B)
{
	int nx=A.Nx();
	int ny=A.Ny();
	assert(&A != &B && nx == B.Ny() && ny == B.Nx());
	
	for(int i=0; i < nx; i++) {
		T *Ai=A[i];
		T *BTi=B()+i;
		for(int j=0; j < ny; j++) {
			Ai[j]=conj(BTi[j*nx]);
		}
	}
	B.Purge();
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
inline Array2<T> operator * (const Array2<T>& B)
{
	Array2<T> A(B.Ny(),B.Nx());
	Conjugate(A,B);
	A.Hold();
	return A;
}

template<class T>
inline void Mult(const Array2<T>& A, const Array2<T>& B, const Array2<T>& C)
{
	int bny=B.Ny();
	assert(bny == C.Nx());

	int n=C.Size();
	static DynVector<T> temp(n);
	if(n > temp.Size()) temp.Resize(n);
	
	int cny=C.Ny();
	Array2<T> CT(cny,C.Nx(),temp);
	
	Trans(CT,C);
	C.Purge();
	
	int bnx=B.Nx();
	if(&A != &B) {
		for(int i=0; i < bnx; i++) {
			Array1(T) Ai=A[i];
			Array1(T) Bi=B[i];
			for(int k=0; k < cny; k++) {
				T sum=0.0;
				Array1(T) CTk=CT[k];
				for(int j=0; j < bny; j++) {
					sum += Bi[j]*CTk[j];
				}
				Ai[k]=sum;
			}
		}
		B.Purge();
	} else {
		static DynVector<T> temp2(cny);
		if(cny > temp2.Size()) temp2.Resize(cny);
		T *work=temp2;
		
		for(int i=0; i < bnx; i++) {
			Array1(T) Ai=A[i];
			int k;
			for(k=0; k < cny; k++) {
				T sum=0.0;
				Array1(T) CTk=CT[k];
				for(int j=0; j < bny; j++) {
					sum += Ai[j]*CTk[j];
				}
				work[k]=sum;
			}
			for(k=0; k < cny; k++) Ai[k]=work[k];
		}
	}		
}

template<class T>
inline void MultAdd(const Array2<T>& A, const Array2<T>& B, const Array2<T>& C,
					const Array2<T>& D)
{
	int bny=B.Ny();
	assert(bny == C.Nx());

	int n=C.Size();
	static DynVector<T> temp(n);
	if(n > temp.Size()) temp.Resize(n);
	
	int cny=C.Ny();
	Array2<T> CT(cny,C.Nx(),temp);
	
	Trans(CT,C);
	C.Purge();
	
	int bnx=B.Nx();
	if(&A != &B) {
		for(int i=0; i < bnx; i++) {
			Array1(T) Ai=A[i];
			Array1(T) Bi=B[i];
			Array1(T) Di=D[i];
			for(int k=0; k < cny; k++) {
				T sum=0.0;
				Array1(T) CTk=CT[k];
				for(int j=0; j < bny; j++) {
					sum += Bi[j]*CTk[j];
				}
				Ai[k]=sum+Di[k];
			}
		}
		B.Purge();
	} else {
		static DynVector<T> temp2(cny);
		if(cny > temp2.Size()) temp2.Resize(cny);
		T *work=temp2;
		
		for(int i=0; i < bnx; i++) {
			Array1(T) Ai=A[i];
			int k;
			for(k=0; k < cny; k++) {
				T sum=0.0;
				Array1(T) CTk=CT[k];
				for(int j=0; j < bny; j++) {
					sum += Ai[j]*CTk[j];
				}
				work[k]=sum;
			}
			Array1(T) Di=D[i];
			for(int k=0; k < cny; k++) Ai[k]=work[k]+Di[k];
		}
	}		
}

template<class T>
inline Array2<T> const operator * (const Array2<T>& B, const Array2<T>& C)
{
	Array2<T> A(B.Nx(),C.Ny());
	Mult(A,B,C);
	A.Hold();
	return A;
}

template<class T>
inline void Mult(const Array2<T>& A, const Array2<T>& B, T C)
{
	int size=A.Size(); 
	for(int i=0; i < size; i++) A(i)=B(i)*C;
	B.Purge();
}

template<class T>
inline Array2<T> operator * (const Array2<T>& B, T C)
{
	Array2<T> A(B.Nx(),B.Ny());
	Mult(A,B,C);
	A.Hold();
	return A;
}

template<class T>
inline Array2<T> operator * (T C, const Array2<T>& B)
{
	Array2<T> A(B.Nx(),B.Ny());
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

template<class T>
inline Array2<T>& Array2<T>::operator *= (T A)
{
	Mult(*this,*this,A);
	return *this;
}


template<class T>
inline void Divide(const Array2<T>& A, const Array2<T>& B, T C)
{
	T Cinv=1.0/C;
	int size=A.Size(); 
	for(int i=0; i < size; i++) A(i)=B(i)*Cinv;
	B.Purge();
}

template<class T>
inline Array2<T> operator / (const Array2<T>& B, T C)
{
	Array2<T> A(B.Nx(),B.Ny());
	Divide(A,B,C);
	A.Hold();
	return A;
}

template<class T>
inline void Add(const Array2<T>& A, const Array2<T>& B, const Array2<T>& C)
{
	assert(A.Nx() == B.Nx() && A.Ny() == B.Ny() && 
		   B.Nx() == C.Nx() && B.Ny() == C.Ny());
	
	int size=A.Size();
	for(int i=0; i < size; i++) A(i)=B(i)+C(i);
	B.Purge();
	C.Purge();
}

template<class T>
inline void Sub(const Array2<T>& A, const Array2<T>& B, const Array2<T>& C)
{
	assert(A.Nx() == B.Nx() && A.Ny() == B.Ny() && 
		   B.Nx() == C.Nx() && B.Ny() == C.Ny());
	
	int size=A.Size();
	for(int i=0; i < size; i++) A(i)=B(i)-C(i);
	B.Purge();
	C.Purge();
}

template<class T>
inline void Unary(const Array2<T>& A, const Array2<T>& B)
{
	assert(A.Nx() == B.Nx() && A.Ny() == B.Ny());
	int size=A.Size();
	for(int i=0; i < size; i++) A(i)=-B(i);
	B.Purge();
}

template<class T>
inline void Add(Array2<T>& A, const Array2<T>& B, T C)
{
	assert(A.Nx() == B.Nx() && A.Ny() == B.Ny());
	int size=A.Size();
	for(int i=0; i < size; i++) A(i)=B(i)+C;
	B.Purge();
}

template<class T>
inline void Sub(Array2<T>& A, const Array2<T>& B, T C)
{
	assert(A.Nx() == B.Nx() && A.Ny() == B.Ny());
	int size=A.Size();
	for(int i=0; i < size; i++) A(i)=B(i)-C;
	B.Purge();
}

template<class T>
inline void Add(const Array2<T>& A, T B, const Array2<T>& C)
{
	assert(A.Nx() == C.Nx() && A.Ny() == C.Ny());
	int size=A.Size();
	for(int i=0; i < size; i++) A(i)=B+C(i);
	C.Purge();
}

template<class T>
inline void Sub(const Array2<T>& A, T B, const Array2<T>& C)
{
	assert(A.Nx() == C.Nx() && A.Ny() == C.Ny());
	int size=A.Size();
	for(int i=0; i < size; i++) A(i)=B-C(i);
	C.Purge();
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
	Array2<T> A(B.Nx(),B.Ny());
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
inline Array2<T> operator + (const Array2<T>& B, const Array2<T>& C)
{
	Array2<T> A(B.Nx(),B.Ny());
	Add(A,B,C);
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
inline Array2<T> operator - (const Array2<T>& B)
{
	Array2<T> A(B.Nx(),B.Ny());
	Unary(A,B);
	A.Hold();
	return A;
}

#endif
