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

#define __ARRAY_H_VERSION__ 1.03J

// Defining NDEBUG improves optimization but disables argument checking.

#ifdef NDEBUG
#define check(i,n,dim,m)
#else
#define check(i,n,dim,m) Check(i,n,dim,m)
#endif

#include "iostream.h"

#if __AIX
#define CONST const
#else
#define mutable
#define CONST
#endif

inline ostream& _newl(ostream& s) {s << '\n'; return s;}

template<class T>
class Array1 {
protected:
	T *v;
	int nx;
	mutable int state;
public:
    enum alloc_state {unallocated=0, allocated=1, temporary=2};
	virtual int Size() const {return nx;}
	int Size0() const {
		if(!test(allocated) && Size() == 0)
			cout << "WARNING: Operation attempted on unallocated array." 
				 << endl;
		return Size();
	}
	
	int test(int flag) const {return state & flag;}
	void clear(int flag) CONST {state &= ~flag;}
	void set(int flag) CONST {state |= flag;}
	void Allocate(int nx0) {Dimension(nx0); v=new T[Size()]; set(allocated);}
	void Deallocate() CONST {delete [] v; clear(allocated);}
	void Dimension(int nx0) {nx=nx0;}
	void Dimension(int nx0, T *v0) {Dimension(nx0); v=v0; clear(allocated);}

	Array1() : nx(0), state(unallocated) {}
	Array1(int nx0) : state(unallocated) {Allocate(nx0);}
	Array1(int nx0, T *v0) : state(unallocated) {Dimension(nx0,v0);}
	Array1(const Array1<T>& A) : v(A.v), nx(A.nx), state(A.test(temporary)) {}
	virtual ~Array1() {if(test(allocated)) Deallocate();}
	
	void Freeze() {state=unallocated;}
	void Hold() {if(test(allocated)) {state=temporary;}}
	void Purge() CONST {if(test(temporary)) {Deallocate(); state=unallocated;}}
#ifdef mutable
	void Purge() const {((Array1<T> *) this)->Purge();}
#endif
	
	void Check(int i, int n, int dim, int m) const {
		if(i < 0 || i >= n) {
			cout << newl << "ERROR: Array" << dim << " index ";
			if(m) cout << m << " ";
			cout << "is out of bounds (" << i;
			if(i < 0) cout << " < " << 0;
			else cout << " > " << n-1;
			cout << ")." << endl;
			exit(1);
		}
	}
	
	int Nx() const {return nx;}
	int N1() const {return nx;}
	T& operator [] (int ix) const {check(ix,nx,1,1); return v[ix];}
	T& operator () (int ix) const {check(ix,nx,1,1); return v[ix];}
	T* operator () () const {return v;}
	operator T* () const {return v;}
	
	Array1<T> operator + (int i) const {return Array1<T>(nx,v+i);}
	
	void Load(T a) const {
		int size=Size0();
		for(int i=0; i < size; i++) v[i]=a;
	}
	void Load(T *a) const {memcpy(v,a,sizeof(T)*Size0());}
	void Store(T *a) const {memcpy(a,v,sizeof(T)*Size0());}
	void Set(T *a) {v=a; clear(allocated);}
	istream& Input (istream &s) const {
		int size=Size0();
		for(int i=0; i < size; i++) s >> v[i];
		return s;
	}
	
	Array1<T>& operator = (T a) {Load(a); return *this;}
	Array1<T>& operator = (T *a) {Load(a); return *this;}
	Array1<T>& operator = (const Array1<T>& A) {Load(A()); return *this;}
	
	Array1<T>& operator += (const Array1<T>& A) {
		int size=Size0(); for(int i=0; i < size; i++) v[i] += A(i);
		return *this;
	}
	Array1<T>& operator -= (const Array1<T>& A) {
		int size=Size0(); for(int i=0; i < size; i++) v[i] -= A(i);
		return *this;
	}
	Array1<T>& operator *= (const Array1<T>& A) {
		int size=Size0(); for(int i=0; i < size; i++) v[i] *= A(i);
		return *this;
	}
	Array1<T>& operator /= (const Array1<T>& A) {
		int size=Size0(); for(int i=0; i < size; i++) v[i] /= A(i);
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
	return A.Input(s);
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
		clear(allocated);
	}
	void Allocate(int nx0, int ny0) {
		Dimension(nx0,ny0);
		v=new T[Size()];
		set(allocated);
	}
	
	Array2() {}
	Array2(int nx0, int ny0) {Allocate(nx0,ny0);}
	Array2(int nx0, int ny0, T *v0) {Dimension(nx0,ny0,v0);}
	
	int Ny() const {return ny;}
	int N2() const {return ny;}
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
		check(i,Size(),2,0);
		return v[i];
	}
	T* operator () () const {return v;}
	
	Array2<T>& operator = (T a) {Load(a); return *this;}
	Array2<T>& operator = (T *a) {Load(a); return *this;}
	Array2<T>& operator = (const Array2<T>& A) {
		Load(A());
		A.Purge();
		return *this;
	}
	
	Array2<T>& operator += (const Array2<T>& A) {
		int size=Size0(); for(int i=0; i < size; i++) v[i] += A(i);
		return *this;
	}
	Array2<T>& operator -= (const Array2<T>& A) {
		int size=Size0(); for(int i=0; i < size; i++) v[i] -= A(i);
		return *this;
	}
	Array2<T>& operator *= (const Array2<T>& A);
	
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
	Array2<T>& operator *= (T A);
};

template<class T>
ostream& operator << (ostream& s, const Array2<T>& A)
{
	T *p=A();
	for(int i=0; i < A.Nx(); i++) {
		for(int j=0; j < A.Ny(); j++) {
			s << *(p++) << " ";
		}
		s << _newl;
	}
	s << flush;
	return s;
}

template<class T>
istream& operator >> (istream& s, const Array2<T>& A)
{
	return A.Input(s);
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
		set(allocated);
	}
	void Dimension(int nx0, int ny0, int nz0) {
		nx=nx0; ny=ny0; nz=nz0; nyz=ny*nz;
	}
	void Dimension(int nx0, int ny0, int nz0, T *v0) {
		Dimension(nx0,ny0,nz0);
		v=v0;
		clear(allocated);
	}
	
	Array3() {}
	Array3(int nx0, int ny0, int nz0) {Allocate(nx0,ny0,nz0);}
	Array3(int nx0, int ny0, int nz0, T *v0) {Dimension(nx0,ny0,nz0,v0);}
	
	int Nz() const {return nz;}
	int N3() const {return nz;}
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
		check(i,Size(),3,0);
		return v[i];
	}
	T* operator () () const {return v;}
	
	Array3<T>& operator = (T a) {Load(a); return *this;}
	Array3<T>& operator = (T *a) {Load(a); return *this;}
	Array3<T>& operator = (const Array3<T>& A) {Load(A()); return *this;}
	
	Array3<T>& operator += (Array3<T>& A) {
		int size=Size0(); for(int i=0; i < size; i++) v[i] += A(i);
		return *this;
	}
	Array3<T>& operator -= (Array3<T>& A) {
		int size=Size0(); for(int i=0; i < size; i++) v[i] -= A(i);
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
			s << _newl;
		}
		s << _newl;
	}
	s << flush;
	return s;
}

template<class T>
istream& operator >> (istream& s, const Array3<T>& A)
{
	return A.Input(s);
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
		set(allocated);
	}
	void Dimension(int nx0, int ny0, int nz0, int nw0) {
		nx=nx0; ny=ny0; nz=nz0; nw=nw0; nzw=nz*nw; nyzw=ny*nzw;
	}
	void Dimension(int nx0, int ny0, int nz0, int w0, T *v0) {
		Dimension(nx0,ny0,nz0,w0);
		v=v0;
		clear(allocated);
	}
	
	Array4() {}
	Array4(int nx0, int ny0, int nz0, int nw0) {Allocate(nx0,ny0,nz0,nw0);}
	Array4(int nx0, int ny0, int nz0, int nw0, T *v0) {
		Dimension(nx0,ny0,nz0,nw0,v0);
	}

	int Nw() const {return nw;}
	int N4() const {return nw;}
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
		check(i,Size(),4,0);
		return v[i];
	}
	T* operator () () const {return v;}
	
	Array4<T>& operator = (T a) {Load(a); return *this;}
	Array4<T>& operator = (T *a) {Load(a); return *this;}
	Array4<T>& operator = (const Array4<T>& A) {Load(A()); return *this;}
	
	Array4<T>& operator += (Array4<T>& A) {
		int size=Size0(); for(int i=0; i < size; i++) v[i] += A(i);
		return *this;
	}
	Array4<T>& operator -= (Array4<T>& A) {
		int size=Size0(); for(int i=0; i < size; i++) v[i] -= A(i);
		return *this;
	}
	
	Array4<T>& operator += (T a) {
		int inc=nyzw+nzw+nw+1, size=Size0();
		for(int i=0; i < size; i += inc) v[i] += a;
		return *this;
	}
	Array4<T>& operator -= (T a) {
		int inc=nyzw+nzw+nw+1, size=Size0();
		for(int i=0; i < size; i += inc) v[i] -= a;
		return *this;
	}
};

template<class T>
ostream& operator << (ostream& s, const Array4<T>& A)
{
	T *p=A;
	for(int i=0; i < A.Nx(); i++) {
		for(int j=0; j < A.Ny(); j++) {
			for(int k=0; k < A.Nz(); k++) {
				for(int l=0; l < A.Nw(); l++) {
					s << *(p++) << " ";
				}
				s << _newl;
			}
			s << _newl;
		}
		s << _newl;
	}	
	s << flush;
	return s;
}

template<class T>
istream& operator >> (istream& s, const Array4<T>& A)
{
	return A.Input(s);
}

#ifdef NDEBUG
#define Array1(T) T*
#else
#define Array1(T) Array1<T>
#endif

#undef check
#undef CONST

#endif

