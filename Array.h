#ifndef __Array_h__
#define __Array_h__ 1

#include "iostream.h"

#ifndef __utils_h__
const char newl='\n';
template<class T> 
inline void set(T *to, const T * from, int n)
{
	memcpy(to,from,sizeof(*from)*n);
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
//	virtual ~Array1() {}
	Array1(T *v0) : v(v0) {}
	Array1(int nx0) {Allocate(nx0);}
	Array1(int nx0, T *v0) {Dimension(nx0,v0);}
	
	int Nx() const {return nx;}
	int N1() const {return nx;}
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
	Array2(int nz, T *v0) : ny(nz) {v=v0;}
	Array2(int nx0, int ny0) {Allocate(nx0,ny0);}
	Array2(int nx0, int ny0, T *v0) {Dimension(nx0,ny0,v0);}
	
	int Ny() const {return ny;}
	int N2() const {return ny;}
	Array1<T> operator [] (int ix) const {return Array1<T>(v+ix*ny);}
	T* operator [] (int ix) {return v+ix*ny;}
	T& operator () (int ix, int iy) const {return v[ix*ny+iy];}
	T& operator () (int i) {return v[i];}
	
	Array2<T>& operator = (T a) {Load(a); return *this;}
	Array2<T>& operator = (T *a) {Load(a); return *this;}
	Array2<T>& operator = (const Array2<T>& A) {Load(A.Base()); return *this;}
	
	Array2<T>& operator += (Array2<T>& A) {
		int size=Size(); for(int i=0; i < size; i++) v[i] += A(i);
		return *this;
	}
	Array2<T>& operator -= (Array2<T>& A) {
		int size=Size(); for(int i=0; i < size; i++) v[i] -= A(i);
		return *this;
	}
	
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
	int nyz;
	int nz;
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
	Array3(int nza, int na, T *v0) : nyz(nza), nz(na) {v=v0;}
	Array3(int nx0, int ny0, int nz0) {Allocate(nx0,ny0,nz0);}
	Array3(int nx0, int ny0, int nz0, T *v0) {Dimension(nx0,ny0,nz0,v0);}
	int Nz() const {return nz;}
	Array2<T> operator [] (int ix) const {return Array2<T>(nz,v+ix*nyz);}
	T& operator () (int ix, int iy, int iz) const {return v[ix*nyz+iy*nz+iz];}
	T& operator () (int i) {return v[i];}
	Array3<T>& operator = (T a) {Load(a); return *this;}
	Array3<T>& operator = (T *a) {Load(a); return *this;}
	Array3<T>& operator = (const Array3<T>& A) {Load(A.Base()); return *this;}
	
	Array3<T>& operator += (Array3<T>& A) {
		int size=Size(); for(int i=0; i < size; i++) v[i] += A(i);
		return *this;
	}
	Array3<T>& operator -= (Array3<T>& A) {
		int size=Size(); for(int i=0; i < size; i++) v[i] -= A(i);
		return *this;
	}
	
	Array3<T>& operator += (T a) {
		int ny1=ny+1, size=Size(); for(int i=0; i < size; i += ny1) v[i] += a;
		return *this;
	}
	Array3<T>& operator -= (T a) {
		int ny1=ny+1, size=Size(); for(int i=0; i < size; i += ny1) v[i] -= a;
		return *this;
	}
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

template<class T>
class Array4 : public Array3<T> {
protected:
	int nyza;
	int nza;
	int na;
public:
	int Size() const {return nx*nyza;}
	void Allocate(int nx0, int ny0, int nz0, int na0) {
		Dimension(nx0,ny0,nz0,na0);
		v=new T[Size()];
	}
	void Dimension(int nx0, int ny0, int nz0, int na0) {
		nx=nx0; ny=ny0; nz=nz0; na=na0; nza=nz*na; nyza=ny*nza;
	}
	void Dimension(int nx0, int ny0, int nz0, int a0, T *v0) {
		Dimension(nx0,ny0,nz0,a0); v=v0;
	}
	Array4() {}
	Array4(int nzab, int nab, int nb, T *v0) : nyza(nzab), nza(nab), na(nb) {
		v=v0;
	}
	Array4(int nx0, int ny0, int nz0, int na0) {Allocate(nx0,ny0,nz0,na0);}
	Array4(int nx0, int ny0, int nz0, int na0, T *v0) {
		Dimension(nx0,ny0,nz0,na0,v0);
	}
	int N4() const {return na;}
	int Na() const {return na;}
	Array3<T> operator [] (int ix) const {return Array3<T>(nza,na,v+ix*nyza);}
	T& operator () (int ix, int iy, int iz, int ia) const {
		return v[ix*nyza+iy*nza+iz*na+ia];
	}
	T& operator () (int i) {return v[i];}
	Array4<T>& operator = (T a) {Load(a); return *this;}
	Array4<T>& operator = (T *a) {Load(a); return *this;}
	Array4<T>& operator = (const Array4<T>& A) {Load(A.Base()); return *this;}
	
	Array4<T>& operator += (Array4<T>& A) {
		int size=Size(); for(int i=0; i < size; i++) v[i] += A(i);
		return *this;
	}
	Array4<T>& operator -= (Array4<T>& A) {
		int size=Size(); for(int i=0; i < size; i++) v[i] -= A(i);
		return *this;
	}
	
	Array4<T>& operator += (T a) {
		int ny1=ny+1, size=Size(); for(int i=0; i < size; i += ny1) v[i] += a;
		return *this;
	}
	Array4<T>& operator -= (T a) {
		int ny1=ny+1, size=Size(); for(int i=0; i < size; i += ny1) v[i] -= a;
		return *this;
	}
};

template<class T>
ostream& operator << (ostream& s, Array4<T>& A)
{
	T *p=A.Base();
	for(int i=0; i < A.Nx(); i++) {
		for(int j=0; j < A.Ny(); j++) {
			for(int k=0; k < A.Nz(); k++) {
				for(int l=0; l < A.Na(); l++) {
					s << *(p++) << " ";
				}
				s << newl;
			}
			s << newl;
		}
		s << newl;
	}	
	s << flush;
	return s;
}

template<class T>
istream& operator >> (istream& s, Array4<T>& A)
{
	T *p=A.Base();
	for(int i=0; i < A.Nx(); i++) {
		for(int j=0; j < A.Ny(); j++) {
			for(int k=0; k < A.Nz(); k++) {
				for(int l=0; l < A.Na(); l++) {
					s >> *(p++);
				}
			}
		}
	}
	return s;
}

#endif

