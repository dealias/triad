#ifndef __Array_h__
#define __Array_h__ 1

#include "iostream.h"

#ifndef __utils_h__
const char newl='\n';

template<class T>
inline T min(T a, T b)
{
	return (a < b) ? a : b;
}

template<class T>
inline T max(T a, T b)
{
	return (a > b) ? a : b;
}
#endif

template<class T>
class Array1 {
protected:
	int nx;
	T *v;
public:
	virtual int Size() const {return nx;}
	void Allocate(int nx0) {Dimension(nx0); v=new T[Size()];}
	void DeAllocate() {delete [] v;}
	void Dimension(int nx0) {nx=nx0;}
	void Dimension(int nx0, T *v0) {Dimension(nx0); v=v0;}
	
	Array1() {}
	virtual ~Array1() {}
	Array1(int nx0) {Allocate(nx0);}
	Array1(int nx0, T *v0) {Dimension(nx0,v0);}
	
	int Nx() const {return nx;}
	T *Base() const {return v;}
	T& operator [] (int ix) {return v[ix];}
	T& operator () (int ix) {return v[ix];}
	
	void Load(const T a) {int size=Size(); for(int i=0; i < size; i++) v[i]=a;}
	void Load(T *a) {set(v,a,Size());}
	void Store(T *a) {set(a,v,Size());}
	void Set(T *a) {v=a;}
	
	Array1<T>& operator = (T a) {Load(a); return *this;}
	Array1<T>& operator = (T *a) {Load(a); return *this;}
	Array1<T>& operator = (const Array1<T>& A) {Load(A.Base()); return *this;}
	
	Array1<T>& operator += (Array1<T>& A) {
		int size=Size(); for(int i=0; i < size; i++) v[i] += A(i);
		return *this;
	}
	Array1<T>& operator -= (Array1<T>& A) {
		int size=Size(); for(int i=0; i < size; i++) v[i] -= A(i);
		return *this;
	}
	
	Array1<T>& operator += (T a) {
		int size=Size(); for(int i=0; i < size; i++) v[i] += a;
		return *this;
	}
	Array1<T>& operator -= (T a) {
		int size=Size(); for(int i=0; i < size; i++) v[i] -= a;
		return *this;
	}
	Array1<T>& operator *= (T a) {
		int size=Size(); for(int i=0; i < size; i++) v[i] *= a;
		return *this;
	}
	Array1<T>& operator /= (T a) {
		int size=Size(); for(int i=0; i < size; i++) v[i] /= a;
		return *this;
	}
	
};

template<class T>
ostream& operator << (ostream& s, Array1<T>& A)
{
	T *p=A.Base();
	for(int i=0; i < A.Nx(); i++) {
		s << *(p++) << " ";
	}
	return s;
}

template<class T>
istream& operator >> (istream& s, Array1<T>& A)
{
	T *p=A.Base();
	for(int i=0; i < A.Nx(); i++) {
		s >> *(p++);
	}
	return s;
}

template<class T>
class Array2 : public Array1<T> {
protected:
	int ny;
public:
	void Dimension(int nx0, int ny0) {nx=nx0; ny=ny0;}
	void Dimension(int nx0, int ny0, T *v0) {Dimension(nx0,ny0); v=v0;}
	int Size() const {return nx*ny;}
	void Allocate(int nx0, int ny0) {Dimension(nx0,ny0); v=new T[Size()];}
	
	Array2() {}
	Array2(int ny0, T *v0) {ny=ny0; v=v0;}
	Array2(int nx0, int ny0) {Allocate(nx0,ny0);}
	Array2(int nx0, int ny0, T *v0) {Dimension(nx0,ny0,v0);}
	
	int Ny() const {return ny;}
	T* operator [] (int ix) const {return v+ix*ny;}
	T& operator () (int i) {return v[i];}
	T& operator () (int ix, int iy) const {return v[ix*ny+iy];}
	
	Array2<T>& operator = (T a) {Load(a); return *this;}
	Array2<T>& operator = (T *a) {Load(a); return *this;}
	Array2<T>& operator = (const Array2<T>& A) {Load(A.Base()); return *this;}
	
	Array2<T>& operator += (T a) {
		int ny1=ny+1, size=Size(); for(int i=0; i < size; i += ny1) v[i] += a;
		return *this;
	}
	Array2<T>& operator -= (T a) {
		int ny1=ny+1, size=Size(); for(int i=0; i < size; i += ny1) v[i] -= a;
		return *this;
	}
};

template<class T>
ostream& operator << (ostream& s, Array2<T>& A)
{
	T *p=A.Base();
	for(int i=0; i < A.Nx(); i++) {
		for(int j=0; j < A.Ny(); j++) {
			s << *(p++) << " ";
		}
		s << newl;
	}
	s << flush;
	return s;
}

template<class T>
istream& operator >> (istream& s, Array2<T>& A)
{
	T *p=A.Base();
	for(int i=0; i < A.Nx(); i++) {
		for(int j=0; j < A.Ny(); j++) {
			s >> *(p++);
		}
	}
	return s;
}

template<class T>
class Array3 : public Array2<T> {
protected:
	int nz;
	int nyz;
public:
	int Size() const {return nx*nyz;}
	void Allocate(int nx0, int ny0, int nz0) {
		Dimension(nx0,ny0,nz0);
		v=new T[Size()];
	}
	void Dimension(int nx0, int ny0, int nz0) {
		nx=nx0; ny=ny0; nz=nz0; nyz=ny*nz;
	}
	void Dimension(int nx0, int ny0, int nz0, T *v0) {
		Dimension(nx0,ny0,nz0); v=v0;
	}
	Array3() {}
	Array3(int nx0, int ny0, int nz0) {Allocate(nx0,ny0,nz0);}
	Array3(int nx0, int ny0, int nz0, T *v0) {Dimension(nx0,ny0,nz0,v0);}
	int Nz() const {return nz;}
	Array2<T> operator [] (int ix) const {return Array2<T>(nz,v+ix*nyz);}
	T& operator () (int i) {return v[i];}
#if 0	
	T& operator () (int ix, int iy, int iz) const {return v[(ix*ny+iy)*nz+iz];}
#else	
	T& operator () (int ix, int iy, int iz) const {return v[ix*nyz+iy*nz+iz];}
#endif	
	Array3<T>& operator = (T a) {Load(a); return *this;}
	Array3<T>& operator = (T *a) {Load(a); return *this;}
	Array3<T>& operator = (const Array3<T>& A) {Load(A.Base()); return *this;}
};

template<class T>
ostream& operator << (ostream& s, Array3<T>& A)
{
	T *p=A.Base();
	for(int i=0; i < A.Nx(); i++) {
		for(int j=0; j < A.Ny(); j++) {
			for(int k=0; k < A.Nz(); k++) {
				s << *(p++) << " ";
			}
			s << newl;
		}
		s << newl;
	}
	s << flush;
	return s;
}

template<class T>
istream& operator >> (istream& s, Array3<T>& A)
{
	T *p=A.Base();
	for(int i=0; i < A.Nx(); i++) {
		for(int j=0; j < A.Ny(); j++) {
			for(int k=0; k < A.Nz(); k++) {
				s >> *(p++);
			}
		}
	}
	return s;
}

#endif
