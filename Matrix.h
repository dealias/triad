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
	Array1(double) a=A;
	for(int i=0; i < size; i += ny1) a[i]=1.0;
}

template<class T>
inline Array2<T> Identity(int n, int m, T)
{
	Array2<T> A(n,m);
	A.Freeze();
	SetIdentity(A);
	return A;
}

template<class T>
inline void Trans(const Array2<T>& A, const Array2<T>& B)
{
	assert(&A != &B);
	for(int i=0; i < B.Nx(); i++) {
		Array1(T) Bi=B[i];
		for(int j=0; j < B.Ny(); j++) {
			A(j,i)=Bi[j];
		}
	}
}

template<class T>
inline void Conjugate(const Array2<T>& A, const Array2<T>& B)
{
	assert(&A != &B);
	for(int i=0; i < B.Nx(); i++) {
		Array1(T) Bi=B[i];
		for(int j=0; j < B.Ny(); j++) {
			A(j,i)=conj(Bi[j]);
		}
	}
}

template<class T>
inline Array2<T> operator * (const Array2<T>& B)
{
	Array2<T> A(B.Ny(),B.Nx());
	A.Freeze();
	Conjugate(A,B);
	return A;
}

template<class T>
inline void Mult(const Array2<T>& A, const Array2<T>& B, const Array2<T>& C)
{
	assert(B.Ny() == C.Nx());

	int n=C.Nx()*C.Ny();
	static DynVector<T> temp(n);
	if(n > temp.Size()) temp.Resize(n);
	
	Array2<T> CT(C.Ny(),C.Nx(),temp());
	
	Trans(CT,C);
	
	if(&A != &B) {
		for(int i=0; i < B.Nx(); i++) {
			Array1(T) Ai=A[i], Bi=B[i];
			for(int k=0; k < C.Ny(); k++) {
				T sum=0.0;
				Array1(T) CTk=CT[k];
				for(int j=0; j < B.Ny(); j++) {
					sum += Bi[j]*CTk[j];
				}
				Ai[k]=sum;
			}
		}
	} else {
		int n=C.Ny();
		static DynVector<T> temp2(n);
		if(n > temp2.Size()) temp2.Resize(n);
		T *work=temp2();
		
		for(int i=0; i < A.Nx(); i++) {
			Array1(T) Ai=A[i];
			int k;
			for(k=0; k < C.Ny(); k++) {
				T sum=0.0;
				Array1(T) CTk=CT[k];
				for(int j=0; j < A.Ny(); j++) {
					sum += Ai[j]*CTk[j];
				}
				work[k]=sum;
			}
			for(k=0; k < C.Ny(); k++) Ai[k]=work[k];
		}
	}		
}

template<class T>
inline void MultAdd(const Array2<T>& A, const Array2<T>& B, const Array2<T>& C,
					const Array2<T>& D)
{
	assert(B.Ny() == C.Nx());

	int n=C.Nx()*C.Ny();
	static DynVector<T> temp(n);
	if(n > temp.Size()) temp.Resize(n);
	
	Array2<T> CT(C.Ny(),C.Nx(),temp());
	
	Trans(CT,C);
	
	if(&A != &B) {
		for(int i=0; i < B.Nx(); i++) {
			Array1(T) Ai=A[i], Bi=B[i], Di=D[i];
			for(int k=0; k < C.Ny(); k++) {
				T sum=0.0;
				Array1(T) CTk=CT[k];
				for(int j=0; j < B.Ny(); j++) {
					sum += Bi[j]*CTk[j];
				}
				Ai[k]=sum+Di[k];
			}
		}
	} else {
		int n=C.Ny();
		static DynVector<T> temp2(n);
		if(n > temp2.Size()) temp2.Resize(n);
		T *work=temp2();
		
		for(int i=0; i < A.Nx(); i++) {
			Array1(T) Ai=A[i];
			for(int k=0; k < C.Ny(); k++) {
				T sum=0.0;
				Array1(T) CTk=CT[k];
				for(int j=0; j < A.Ny(); j++) {
					sum += Ai[j]*CTk[j];
				}
				work[k]=sum;
			}
			Array1a(T) Di=D[i]; // Work around GCC-2.8.1 typedef bug
			for(int k=0; k < C.Ny(); k++) Ai[k]=work[k]+Di[k];
		}
	}		
}

template<class T>
inline Array2<T> operator * (const Array2<T>& B, const Array2<T>& C)
{
	Array2<T> A(B.Nx(),C.Ny());
	A.Freeze();
	Mult(A,B,C);
	return A;
}

template<class T>
inline void Mult(const Array2<T>& A, const Array2<T>& B, T C)
{
	int size=A.Size(); 
	Array1(double) a=A, b=B;
	for(int i=0; i < size; i++) a[i]=b[i]*C;
}

template<class T>
inline Array2<T> operator * (const Array2<T>& B, T C)
{
	Array2<T> A(B.Nx(),B.Ny());
	A.Freeze();
	Mult(A,B,C);
	return A;
}

template<class T>
inline Array2<T> operator * (T C, const Array2<T>& B)
{
	Array2<T> A(B.Nx(),B.Ny());
	A.Freeze();
	Mult(A,B,C);
	return A;
}

template<class T>
inline void Divide(const Array2<T>& A, const Array2<T>& B, T C)
{
	int size=A.Size(); 
	T Cinv=1.0/C;
	Array1(double) a=A, b=B;
	for(int i=0; i < size; i++) a[i]=b[i]*Cinv;
}

template<class T>
inline Array2<T> operator / (const Array2<T>& B, T C)
{
	Array2<T> A(B.Nx(),B.Ny());
	A.Freeze();
	Divide(A,B,C);
	return A;
}

template<class T>
inline void Add(Array2<T>& A, const Array2<T>& B, T C)
{
	A=B;
	A += C;
}

template<class T>
inline Array2<T> operator + (const Array2<T>& B, T C)
{
	Array2<T> A(B.Nx(),B.Ny());
	A.Freeze();
	Add(A,B,C);
	return A;
}

template<class T>
inline Array2<T> operator + (T C, const Array2<T>& B)
{
	Array2<T> A(B.Nx(),B.Ny());
	A.Freeze();
	Add(A,B,C);
	return A;
}

template<class T>
inline void Add(const Array2<T>& A, const Array2<T>& B, const Array2<T>& C)
{
	for(int i=0; i < B.Nx(); i++) {
		Array1(T) Ai=A[i], Bi=B[i], Ci=C[i];
		for(int j=0; j < B.Ny(); j++) {
			Ai[j]=Bi[j]+Ci[j];
		}
	}
}

template<class T>
inline Array2<T> operator + (const Array2<T>& B, const Array2<T>& C)
{
	if(B.Nx() != C.Nx() || B.Ny() != C.Ny())
		msg(ERROR, "Incompatible matrices");
	
	Array2<T> A(B.Nx(),B.Ny());
	A.Freeze();
	Add(A,B,C);
	return A;
}

template<class T>
inline void Unary(const Array2<T>& A, const Array2<T>& B)
{
	int size=A.Size(); 
	Array1(double) a=A, b=B;
	for(int i=0; i < size; i++) a[i]=-b[i];
}

template<class T>
inline Array2<T> operator - (const Array2<T>& B)
{
	Array2<T> A(B.Nx(),B.Ny());
	A.Freeze();
	Unary(A,B);
	return A;
}

template<class T>
inline void Sub(Array2<T>& A, const Array2<T>& B, T C)
{
	A=B;
	A -= C;
}

template<class T>
inline Array2<T> operator - (const Array2<T>& B, T C)
{
	Array2<T> A(B.Nx(),B.Ny());
	A.Freeze();
	Sub(A,B,C);
	return A;
}

template<class T>
inline void Sub(const Array2<T>& A, T B, const Array2<T>& C)
{
	Unary(A,C);
	A += B;
}

template<class T>
inline Array2<T> operator - (T B, const Array2<T>& C)
{
	Array2<T> A(C.Nx(),C.Ny());
	A.Freeze();
	Sub(A,B,C);
	return A;
}

template<class T>
inline void Sub(const Array2<T>& A, const Array2<T>& B, const Array2<T>& C)
{
	for(int i=0; i < B.Nx(); i++) {
		Array1(T) Ai=A[i], Bi=B[i], Ci=C[i];
		for(int j=0; j < B.Ny(); j++) {
			Ai[j]=Bi[j]-Ci[j];
		}
	}
}

template<class T>
inline Array2<T> operator - (const Array2<T>& B, const Array2<T>& C)
{
	if(B.Nx() != C.Nx() || B.Ny() != C.Ny())
		msg(ERROR, "Incompatible matrices");
	
	Array2<T> A(B.Nx(),B.Ny());
	A.Freeze();
	Sub(A,B,C);
	return A;
}

#endif
