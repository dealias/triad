/* Array.h:  A high-performance multi-dimensional C++ array class
Copyright (C) 1998 John C. Bowman (bowman@math.ualberta.ca)

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

#ifndef __Array_h__
#define __Array_h__ 1

#define __ARRAY_H_VERSION__ 1.01

// Setting ARRAY_CHECK to 1 enables optional argument checking.

#ifndef ARRAY_CHECK
#define ARRAY_CHECK 0
#endif

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
	T *v;
	int nx;
	mutable int allocate;
	mutable int temporary;
public:
	virtual int Size() const {return nx;}
	int Size0() {
#if ARRAY_CHECK
		if(!allocate && Size() == 0)
			cout << "WARNING: Operation attempted on unallocated array." 
				 << endl;
#endif		
		return Size();
	}
	
	void Allocate(int nx0) {Dimension(nx0); v=new T[Size()]; allocate=1;}
	void Deallocate() const {delete [] v; allocate=0;}
	void Dimension(int nx0) {nx=nx0;}
	void Dimension(int nx0, T *v0) {Dimension(nx0); v=v0; allocate=0;}
	
	Array1() : nx(0), allocate(0), temporary(0) {}
	Array1(int nx0) : temporary(0) {Allocate(nx0);}
	Array1(int nx0, T *v0) : temporary(0) {Dimension(nx0,v0);}
	Array1(const Array1<T>& A) : 
		v(A.v), nx(A.nx), allocate(0), temporary(A.temporary){}
	virtual ~Array1() {if(allocate) Deallocate();}
	
	void Freeze() {allocate=0;}
	void Hold() {if(allocate) {temporary=1; allocate=0;}}
	void Purge() const {if(temporary) {Deallocate(); temporary=0;}}
	
	void check(int i, int n, int dim, int m=0) const {
#if ARRAY_CHECK		
		if(i < 0 || i >= n) {
			cout << "ERROR: Array" << dim << " index ";
			if(m) cout << m << " ";
			cout << "is out of bounds (" << i;
			if(i < 0) cout << " < " << 0;
			else cout << " > " << n-1;
			cout << ")." << endl;
			exit(1);
		}
#endif		
	}
	
	int Nx() const {return nx;}
	int N1() const {return nx;}
	Array1<T> operator () () const {return Array1<T>(Size(),v);}
	T& operator [] (int ix) const {check(ix,nx,1,1); return v[ix];}
	T& operator () (int ix) const {check(ix,nx,1,1); return v[ix];}
	operator T* () const {return v;}
	
	void Load(T a) {
		int size=Size0();
		for(int i=0; i < size; i++) v[i]=a;
	}
	void Load(T *a) {set(v,a,Size0());}
	void Store(T *a) {set(a,v,Size0());}
	void Set(T *a) {v=a; allocate=0;}
	
	Array1<T>& operator = (T a) {Load(a); return *this;}
	Array1<T>& operator = (T *a) {Load(a); return *this;}
	Array1<T>& operator = (const Array1<T>& A) {Load(A()); return *this;}
	
	Array1<T>& operator += (const Array1<T>& A) {
		int size=Size0(); for(int i=0; i < size; i++) v[i] += A()[i];
		return *this;
	}
	Array1<T>& operator -= (const Array1<T>& A) {
		int size=Size0(); for(int i=0; i < size; i++) v[i] -= A()[i];
		return *this;
	}
	
	Array1<T>& operator += (T a) {
		int size=Size0(); for(int i=0; i < size; i++) v[i] += a;
		return *this;
	}
	Array1<T>& operator -= (T a) {
		int size=Size0(); for(int i=0; i < size; i++) v[i] -= a;
		return *this;
	}
	Array1<T>& operator *= (T a) {
		int size=Size0(); for(int i=0; i < size; i++) v[i] *= a;
		return *this;
	}
	Array1<T>& operator /= (T a) {
		int size=Size0(); for(int i=0; i < size; i++) v[i] /= a;
		return *this;
	}
	
};

template<class T>
ostream& operator << (ostream& s, const Array1<T>& A)
{
	T *p=A();
	for(int i=0; i < A.Nx(); i++) {
		s << *(p++) << " ";
	}
	return s;
}

template<class T>
istream& operator >> (istream& s, const Array1<T>& A)
{
	T *p=A();
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
	int Size() const {return nx*ny;}
	void Dimension(int nx0, int ny0) {nx=nx0; ny=ny0;}
	void Dimension(int nx0, int ny0, T *v0) {
		Dimension(nx0,ny0);
		v=v0;
		allocate=0;
	}
	void Allocate(int nx0, int ny0) {
		Dimension(nx0,ny0);
		v=new T[Size()];
		allocate=1;
	}
	
	Array2() {}
	Array2(int nx0, int ny0) {Allocate(nx0,ny0);}
	Array2(int nx0, int ny0, T *v0) {Dimension(nx0,ny0,v0);}
	
	int Ny() const {return ny;}
	int N2() const {return ny;}
	Array1<T> operator () () const {return Array1<T>(Size(),v);}
	Array1<T> operator [] (int ix) const {
		check(ix,nx,2,1);
		return Array1<T>(ny,v+ix*ny);
	}
	T& operator () (int ix, int iy) const {
		check(ix,nx,2,1);
		check(iy,ny,2,2);
		return v[ix*ny+iy];
	}
	T& operator () (int i) const {
		check(i,Size(),2);
		return v[i];
	}
	
	Array2<T>& operator = (T a) {Load(a); return *this;}
	Array2<T>& operator = (T *a) {Load(a); return *this;}
	Array2<T>& operator = (const Array2<T>& A) {
		Load(A);
		A.Purge();
		return *this;
	}
	
	Array2<T>& operator += (Array2<T>& A) {
		int size=Size0(); for(int i=0; i < size; i++) v[i] += A()[i];
		return *this;
	}
	Array2<T>& operator -= (Array2<T>& A) {
		int size=Size0(); for(int i=0; i < size; i++) v[i] -= A()[i];
		return *this;
	}
	
	Array2<T>& operator += (T a) {
		int inc=ny+1, size=Size0();
		for(int i=0; i < size; i += inc) v[i] += a;
		return *this;
	}
	Array2<T>& operator -= (T a) {
		int inc=ny+1, size=Size0();
		for(int i=0; i < size; i += inc) v[i] -= a;
		return *this;
	}
};

template<class T>
ostream& operator << (ostream& s, const Array2<T>& A)
{
	T *p=A();
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
istream& operator >> (istream& s, const Array2<T>& A)
{
	T *p=A();
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
		allocate=1;
	}
	void Dimension(int nx0, int ny0, int nz0) {
		nx=nx0; ny=ny0; nz=nz0; nyz=ny*nz;
	}
	void Dimension(int nx0, int ny0, int nz0, T *v0) {
		Dimension(nx0,ny0,nz0);
		v=v0;
		allocate=0;
	}
	
	Array3() {}
	Array3(int nx0, int ny0, int nz0) {Allocate(nx0,ny0,nz0);}
	Array3(int nx0, int ny0, int nz0, T *v0) {Dimension(nx0,ny0,nz0,v0);}
	
	int Nz() const {return nz;}
	int N3() const {return nz;}
	Array1<T> operator () () const {return Array1<T>(Size(),v);}
	Array2<T> operator [] (int ix) const {
		check(ix,nx,3,1);
		return Array2<T>(ny,nz,v+ix*nyz);
	}
	T& operator () (int ix, int iy, int iz) const {
		check(ix,nx,3,1);
		check(iy,ny,3,2);
		check(iz,nz,3,3);
		return v[ix*nyz+iy*nz+iz];
	}
	T& operator () (int i) const {
		check(i,Size(),3);
		return v[i];
	}
	
	Array3<T>& operator = (T a) {Load(a); return *this;}
	Array3<T>& operator = (T *a) {Load(a); return *this;}
	Array3<T>& operator = (const Array3<T>& A) {Load(A()); return *this;}
	
	Array3<T>& operator += (Array3<T>& A) {
		int size=Size0(); for(int i=0; i < size; i++) v[i] += A()[i];
		return *this;
	}
	Array3<T>& operator -= (Array3<T>& A) {
		int size=Size0(); for(int i=0; i < size; i++) v[i] -= A()[i];
		return *this;
	}
	
	Array3<T>& operator += (T a) {
		int inc=nyz+nz+1, size=Size0();
		for(int i=0; i < size; i += inc) v[i] += a;
		return *this;
	}
	Array3<T>& operator -= (T a) {
		int inc=nyz+nz+1, size=Size0();
		for(int i=0; i < size; i += inc) v[i] -= a;
		return *this;
	}
};

template<class T>
ostream& operator << (ostream& s, const Array3<T>& A)
{
	T *p=A();
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
istream& operator >> (istream& s, const Array3<T>& A)
{
	T *p=A();
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
	int nyzw;
	int nzw;
	int nw;
public:
	int Size() const {return nx*nyzw;}
	void Allocate(int nx0, int ny0, int nz0, int nw0) {
		Dimension(nx0,ny0,nz0,nw0);
		v=new T[Size()];
		allocate=1;
	}
	void Dimension(int nx0, int ny0, int nz0, int nw0) {
		nx=nx0; ny=ny0; nz=nz0; nw=nw0; nzw=nz*nw; nyzw=ny*nzw;
	}
	void Dimension(int nx0, int ny0, int nz0, int w0, T *v0) {
		Dimension(nx0,ny0,nz0,w0);
		v=v0;
		allocate=0;
	}
	
	Array4() {}
	Array4(int nx0, int ny0, int nz0, int nw0) {Allocate(nx0,ny0,nz0,nw0);}
	Array4(int nx0, int ny0, int nz0, int nw0, T *v0) {
		Dimension(nx0,ny0,nz0,nw0,v0);
	}

	int Nw() const {return nw;}
	int N4() const {return nw;}
	Array1<T> operator () () const {return Array1<T>(Size(),v);}
	Array3<T> operator [] (int ix) const {
		check(ix,nx,3,1);
		return Array3<T>(ny,nz,nw,v+ix*nyzw);
	}
	T& operator () (int ix, int iy, int iz, int iw) const {
		check(ix,nx,4,1);
		check(iy,ny,4,2);
		check(iz,nz,4,3);
		check(iw,nw,4,4);
		return v[ix*nyzw+iy*nzw+iz*nw+iw];
	}
	T& operator () (int i) const {
		check(i,Size(),4);
		return v[i];
	}
	
	Array4<T>& operator = (T a) {Load(a); return *this;}
	Array4<T>& operator = (T *a) {Load(a); return *this;}
	Array4<T>& operator = (const Array4<T>& A) {Load(A()); return *this;}
	
	Array4<T>& operator += (Array4<T>& A) {
		int size=Size0(); for(int i=0; i < size; i++) v[i] += A()[i];
		return *this;
	}
	Array4<T>& operator -= (Array4<T>& A) {
		int size=Size0(); for(int i=0; i < size; i++) v[i] -= A()[i];
		return *this;
	}
	
	Array4<T>& operator += (T a) {
		int inc=nyzw+nzw+nw+1, size=Size0();
		for(int i=0; i < size; i += inc) v[i] += a;
		return *this;
	}
	Array4<T>& operator -= (T a) {
		int inc=nyzw+nzw+nw+1, size=Size0();
		for(int i=0; i < size; i += inc) v[i] -= w;
		return *this;
	}
};

template<class T>
ostream& operator << (ostream& s, const Array4<T>& A)
{
	T *p=A();
	for(int i=0; i < A.Nx(); i++) {
		for(int j=0; j < A.Ny(); j++) {
			for(int k=0; k < A.Nz(); k++) {
				for(int l=0; l < A.Nw(); l++) {
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
istream& operator >> (istream& s, const Array4<T>& A)
{
	T *a=A();
	for(int i=0; i < A.Size(); i++) s >> a[i];
	return s;
}

#if ARRAY_CHECK
#define Array1(T) Array1<T>
#define Array1a(T) Array1<T>
#define Array1b(T) Array1<T>
#define Array1c(T) Array1<T>
#define Array1d(T) Array1<T>
#define Array1e(T) Array1<T>
#else
#define Array1(T) typedef T* Tstar; Tstar
#define Array1a(T) typedef T* Tstara; Tstara
#define Array1b(T) typedef T* Tstarb; Tstarb
#define Array1c(T) typedef T* Tstarc; Tstarc
#define Array1d(T) typedef T* Tstard; Tstard
#define Array1e(T) typedef T* Tstare; Tstare
#endif

#endif

