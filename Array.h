/* Array.h:  A high-performance multi-dimensional C++ array class
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

#ifndef __Array_h__
#define __Array_h__ 1

#define __ARRAY_H_VERSION__ 1.07J	

// Defining NDEBUG improves optimization but disables argument checking.

#ifdef NDEBUG
#define check(i,n,o,dim,m)
#else
#define check(i,n,o,dim,m) Check(i,n,o,dim,m)
#endif

#include <iostream.h>
#include <strstream.h>
#include <stdlib.h>

#if __AIX
#define CONST const
#else
#define mutable
#define CONST
#endif

inline ostream& _newl(ostream& s) {s << '\n'; return s;}

#ifndef __ARRAY_EXIT
#define __ARRAY_EXIT(x) cout << _newl << "ERROR: " << x << "." << endl; exit(1)
#endif

template<class T>
class Array1 {
protected:
	T *v; // Pointer to first element of contiguous memory block
	T *voff; // Offset pointer to memory block
	unsigned int nx;
	int ox;
	mutable int state;
public:
    enum alloc_state {unallocated=0, allocated=1, temporary=2};
	virtual unsigned int Size() const {return nx;}
	unsigned int Size0() const {
		if(!test(allocated) && Size() == 0)
			cout << "WARNING: Operation attempted on unallocated array." 
				 << endl;
		return Size();
	}
	
	int test(int flag) const {return state & flag;}
	void clear(int flag) CONST {state &= ~flag;}
	void set(int flag) CONST {state |= flag;}
	void Offsets() {
		voff=v-ox;
	}
	void Dimension(unsigned int nx0, int ox0=0) {
		nx=nx0;
		ox=ox0;
	}
	void Dimension(unsigned int nx0, T *v0, int ox0=0) {
		Dimension(nx0,ox0);
		v=v0;
		Offsets();
		clear(allocated);
	}
	void Allocate(unsigned int nx0, int ox0=0) {
		Dimension(nx0,ox0);
		v=new T[Size()];
		Offsets();
		set(allocated);
	}
	void Deallocate() CONST {delete [] v; clear(allocated);}

	Array1() : nx(0), ox(0), state(unallocated) {}
	Array1(unsigned int nx0, int ox0=0) : state(unallocated) {
		Allocate(nx0,ox0);
	}
	Array1(unsigned int nx0, T *v0, int ox0=0) : state(unallocated) {
		Dimension(nx0,v0,ox0);
	}
	Array1(const Array1<T>& A) : v(A.v), voff(A.voff), nx(A.nx), ox(A.ox),
		state(A.test(temporary)) {}
	virtual ~Array1() {if(test(allocated)) Deallocate();}
	
	void Freeze() {state=unallocated;}
	void Hold() {if(test(allocated)) {state=temporary;}}
	void Purge() CONST {if(test(temporary)) {Deallocate(); state=unallocated;}}
#ifdef mutable
	void Purge() const {((Array1<T> *) this)->Purge();}
#endif
	
	virtual void Check(int& i, int n, int o, 
					   unsigned int dim, unsigned int m) const {
		if(i < o || i >= n+o) {
			strstream buf;
			buf << "Array" << dim << " index ";
			if(m) buf << m << " ";
			buf << "is out of bounds (" << i;
			if(i < o) buf << " < " << o;
			else buf << " > " << n+o-1;
			buf << ")" << ends;
			__ARRAY_EXIT(buf.str());
		}
	}
	
	unsigned int Nx() const {return nx;}
	unsigned int N1() const {return nx;}
	T& operator [] (int ix) const {check(ix,nx,ox,1,1); return voff[ix];}
	T& operator () (int i) const {check(i,nx,0,1,1); return v[i];}
	T* operator () () const {return voff;}
	operator T* () const {return voff;}
	
	Array1<T> operator + (int i) const {return Array1<T>(nx-i,v+i,ox+i);}
	
	void Load(T a) const {
		unsigned int size=Size0();
		for(unsigned int i=0; i < size; i++) v[i]=a;
	}
	void Load(const T *a) const {memcpy(v,a,sizeof(T)*Size0());}
	void Store(T *a) const {memcpy(a,v,sizeof(T)*Size0());}
	void Set(T *a) {v=a; Offsets(); clear(allocated);}
	istream& Input (istream &s) const {
		unsigned int size=Size0();
		for(unsigned int i=0; i < size; i++) s >> v[i];
		return s;
	}
	
	Array1<T>& operator = (T a) {Load(a); return *this;}
	Array1<T>& operator = (const T *a) {Load(a); return *this;}
	Array1<T>& operator = (const Array1<T>& A) {
		Load(A());
		A.Purge();
		return *this;
	}
	
	Array1<T>& operator += (const Array1<T>& A) {
		unsigned int size=Size0(); for(unsigned int i=0; i < size; i++) v[i] += A(i);
		return *this;
	}
	Array1<T>& operator -= (const Array1<T>& A) {
		unsigned int size=Size0(); for(unsigned int i=0; i < size; i++) v[i] -= A(i);
		return *this;
	}
	Array1<T>& operator *= (const Array1<T>& A) {
		unsigned int size=Size0(); for(unsigned int i=0; i < size; i++) v[i] *= A(i);
		return *this;
	}
	Array1<T>& operator /= (const Array1<T>& A) {
		unsigned int size=Size0(); for(unsigned int i=0; i < size; i++) v[i] /= A(i);
		return *this;
	}
	
	Array1<T>& operator += (T a) {
		unsigned int size=Size0(); for(unsigned int i=0; i < size; i++) v[i] += a;
		return *this;
	}
	Array1<T>& operator -= (T a) {
		unsigned int size=Size0(); for(unsigned int i=0; i < size; i++) v[i] -= a;
		return *this;
	}
	Array1<T>& operator *= (T a) {
		unsigned int size=Size0(); for(unsigned int i=0; i < size; i++) v[i] *= a;
		return *this;
	}
	Array1<T>& operator /= (T a) {
		unsigned int size=Size0(); for(unsigned int i=0; i < size; i++) v[i] /= a;
		return *this;
	}
	
};

template<class T>
ostream& operator << (ostream& s, const Array1<T>& A)
{
	T *p=A();
	for(unsigned int i=0; i < A.Nx(); i++) {
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
	unsigned int ny;
	int oy;
	T *vtemp;
public:
	unsigned int Size() const {return nx*ny;}
	void Offsets() {
		vtemp=v-ox*ny;
		voff=vtemp-oy;
	}
	void Dimension(unsigned int nx0, unsigned int ny0, int ox0=0, int oy0=0) {
		nx=nx0; ny=ny0;
		ox=ox0; oy=oy0;
	}
	void Dimension(unsigned int nx0, unsigned int ny0, T *v0, int ox0=0,
				   int oy0=0) {
		Dimension(nx0,ny0,ox0,oy0);
		v=v0;
		Offsets();
		clear(allocated);
	}
	void Allocate(unsigned int nx0, unsigned int ny0, int ox0=0, int oy0=0) {
		Dimension(nx0,ny0,ox0,oy0);
		v=new T[Size()];
		Offsets();
		set(allocated);
	}

	Array2() : ny(0), oy(0) {}
	Array2(unsigned int nx0, unsigned int ny0, int ox0=0, int oy0=0) {
		Allocate(nx0,ny0,ox0,oy0);
	}
	Array2(unsigned int nx0, unsigned int ny0, T *v0, int ox0=0, int oy0=0) {
		Dimension(nx0,ny0,v0,ox0,oy0);
	}
	
	unsigned int Ny() const {return ny;}
	unsigned int N2() const {return ny;}
	Array1<T> operator [] (int ix) const {
		check(ix,nx,ox,2,1);
		return Array1<T>(ny,vtemp+ix*ny,oy);
	}
	T& operator () (int ix, int iy) const {
		check(ix,nx,ox,2,1);
		check(iy,ny,oy,2,2);
		return voff[ix*ny+iy];
	}
	T& operator () (int i) const {
		check(i,Size(),0,2,0);
		return v[i];
	}
	T* operator () () const {return voff;}
	
	Array2<T>& operator = (T a) {Load(a); return *this;}
	Array2<T>& operator = (T *a) {Load(a); return *this;}
	Array2<T>& operator = (const Array2<T>& A) {
		Load(A());
		A.Purge();
		return *this;
	}
	void Set(T *a) {v=a; Offsets(); clear(allocated);}
	
	Array2<T>& operator += (const Array2<T>& A) {
		unsigned int size=Size0();
		for(unsigned int i=0; i < size; i++) v[i] += A(i);
		return *this;
	}
	Array2<T>& operator -= (const Array2<T>& A) {
		unsigned int size=Size0();
		for(unsigned int i=0; i < size; i++) v[i] -= A(i);
		return *this;
	}
	Array2<T>& operator *= (const Array2<T>& A);
	
	Array2<T>& operator += (T a) {
		unsigned int inc=ny+1, size=Size0();
		for(unsigned int i=0; i < size; i += inc) v[i] += a;
		return *this;
	}
	Array2<T>& operator -= (T a) {
		unsigned int inc=ny+1, size=Size0();
		for(unsigned int i=0; i < size; i += inc) v[i] -= a;
		return *this;
	}
	Array2<T>& operator *= (T a) {
		unsigned int size=Size0();
		for(unsigned int i=0; i < size; i++) v[i] *= a;
		return *this;
	}
};

template<class T>
ostream& operator << (ostream& s, const Array2<T>& A)
{
	T *p=A();
	for(unsigned int i=0; i < A.Nx(); i++) {
		for(unsigned int j=0; j < A.Ny(); j++) {
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
	unsigned int nyz;
	unsigned int nz;
	int oz;
public:
	unsigned int Size() const {return nx*nyz;}
	void Offsets() {
		vtemp=v-ox*nyz;
		voff=vtemp-oy*nz-oz;
	}
	void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
				   int ox0=0, int oy0=0, int oz0=0) {
		nx=nx0; ny=ny0; nz=nz0;
		ox=ox0; oy=oy0; oz=oz0;
		nyz=ny*nz;
	}
	void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
				   T *v0, int ox0=0, int oy0=0, int oz0=0) {
		Dimension(nx0,ny0,nz0,ox0,oy0,oz0);
		v=v0;
		Offsets();
		clear(allocated);
	}
	void Allocate(unsigned int nx0, unsigned int ny0, unsigned int nz0,
				  int ox0=0, int oy0=0, int oz0=0) {
		Dimension(nx0,ny0,nz0,ox0,oy0,oz0);
		v=new T[Size()];
		Offsets();
		set(allocated);
	}
	
	Array3() : nz(0), oz(0) {}
	Array3(unsigned int nx0, unsigned int ny0, unsigned int nz0,
		   int ox0=0, int oy0=0, int oz0=0) {
		Allocate(nx0,ny0,nz0,ox0,oy0,oz0);
	}
	Array3(unsigned int nx0, unsigned int ny0, unsigned int nz0, T *v0,
		   int ox0=0, int oy0=0, int oz0=0) {
		Dimension(nx0,ny0,nz0,v0,ox0,oy0,oz0);
	}
	
	unsigned int Nz() const {return nz;}
	unsigned int N3() const {return nz;}
	Array2<T> operator [] (int ix) const {
		check(ix,nx,ox,3,1);
		return Array2<T>(ny,nz,vtemp+ix*nyz,oy,oz);
	}
	T& operator () (int ix, int iy, int iz) const {
		check(ix,nx,ox,3,1);
		check(iy,ny,oy,3,2);
		check(iz,nz,oz,3,3);
		return voff[ix*nyz+iy*nz+iz];
	}
	T& operator () (int i) const {
		check(i,Size(),0,3,0);
		return v[i];
	}
	T* operator () () const {return voff;}
	
	Array3<T>& operator = (T a) {Load(a); return *this;}
	Array3<T>& operator = (T *a) {Load(a); return *this;}
	Array3<T>& operator = (const Array3<T>& A) {
		Load(A());
		A.Purge(); 
		return *this;
	}
	void Set(T *a) {v=a; Offsets(); clear(allocated);}
	
	Array3<T>& operator += (Array3<T>& A) {
		unsigned int size=Size0();
		for(unsigned int i=0; i < size; i++) v[i] += A(i);
		return *this;
	}
	Array3<T>& operator -= (Array3<T>& A) {
		unsigned int size=Size0();
		for(unsigned int i=0; i < size; i++) v[i] -= A(i);
		return *this;
	}
	
	Array3<T>& operator += (T a) {
		unsigned int inc=nyz+nz+1, size=Size0();
		for(unsigned int i=0; i < size; i += inc) v[i] += a;
		return *this;
	}
	Array3<T>& operator -= (T a) {
		unsigned int inc=nyz+nz+1, size=Size0();
		for(unsigned int i=0; i < size; i += inc) v[i] -= a;
		return *this;
	}
};

template<class T>
ostream& operator << (ostream& s, const Array3<T>& A)
{
	T *p=A();
	for(unsigned int i=0; i < A.Nx(); i++) {
		for(unsigned int j=0; j < A.Ny(); j++) {
			for(unsigned int k=0; k < A.Nz(); k++) {
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
	unsigned int nyzw;
	unsigned int nzw;
	unsigned int nw;
	int ow;
public:
	unsigned int Size() const {return nx*nyzw;}
	void Offsets() {
		vtemp=v-ox*nyzw;
		voff=vtemp-oy*nzw-oz*nw-ow;
	}
	void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
				   unsigned int nw0, 
				   int ox0=0, int oy0=0, int oz0=0, int ow0=0) {
		nx=nx0; ny=ny0; nz=nz0; nw=nw0;
		ox=ox0; oy=oy0; oz=oz0; ow=ow0;
		nzw=nz*nw; nyzw=ny*nzw;
	}
	void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
				   unsigned int w0, T *v0,
				   int ox0=0, int oy0=0, int oz0=0, int ow0=0) {
		Dimension(nx0,ny0,nz0,w0,ox0,oy0,oz0,ow0);
		v=v0;
		Offsets();
		clear(allocated);
	}
	void Allocate(unsigned int nx0, unsigned int ny0, unsigned int nz0,
				  unsigned int nw0,
				  int ox0=0, int oy0=0, int oz0=0, int ow0=0) {
		Dimension(nx0,ny0,nz0,nw0,ox0,oy0,oz0,ow0);
		v=new T[Size()];
		Offsets();
		set(allocated);
	}
	
	Array4() : nw(0), ow(0) {}
	Array4(unsigned int nx0, unsigned int ny0, unsigned int nz0,
		   unsigned int nw0,
		   int ox0=0, int oy0=0, int oz0=0, int ow0=0) {
		Allocate(nx0,ny0,nz0,nw0,ox0,oy0,oz0,ow0);
	}
	Array4(unsigned int nx0, unsigned int ny0, unsigned int nz0,
		   unsigned int nw0, T *v0,
		   int ox0=0, int oy0=0, int oz0=0, int ow0=0) {
		Dimension(nx0,ny0,nz0,nw0,v0,ox0,oy0,oz0,ow0);
	}

	unsigned int Nw() const {return nw;}
	unsigned int N4() const {return nw;}
	Array3<T> operator [] (int ix) const {
		check(ix,nx,ox,3,1);
		return Array3<T>(ny,nz,nw,vtemp+ix*nyzw,oy,oz,ow);
	}
	T& operator () (int ix, int iy, int iz, int iw) const {
		check(ix,nx,ox,4,1);
		check(iy,ny,oy,4,2);
		check(iz,nz,oz,4,3);
		check(iw,nw,ow,4,4);
		return voff[ix*nyzw+iy*nzw+iz*nw+iw];
	}
	T& operator () (int i) const {
		check(i,Size(),0,4,0);
		return v[i];
	}
	T* operator () () const {return voff;}
	
	Array4<T>& operator = (T a) {Load(a); return *this;}
	Array4<T>& operator = (T *a) {Load(a); return *this;}
	Array4<T>& operator = (const Array4<T>& A) {
		Load(A());
		A.Purge();
		return *this;
	}
	
	void Set(T *a) {v=a; Offsets(); clear(allocated);}
	
	Array4<T>& operator += (Array4<T>& A) {
		unsigned int size=Size0();
		for(unsigned int i=0; i < size; i++) v[i] += A(i);
		return *this;
	}
	Array4<T>& operator -= (Array4<T>& A) {
		unsigned int size=Size0();
		for(unsigned int i=0; i < size; i++) v[i] -= A(i);
		return *this;
	}
	
	Array4<T>& operator += (T a) {
		unsigned int inc=nyzw+nzw+nw+1, size=Size0();
		for(unsigned int i=0; i < size; i += inc) v[i] += a;
		return *this;
	}
	Array4<T>& operator -= (T a) {
		unsigned int inc=nyzw+nzw+nw+1, size=Size0();
		for(unsigned int i=0; i < size; i += inc) v[i] -= a;
		return *this;
	}
};

template<class T>
ostream& operator << (ostream& s, const Array4<T>& A)
{
	T *p=A;
	for(unsigned int i=0; i < A.Nx(); i++) {
		for(unsigned int j=0; j < A.Ny(); j++) {
			for(unsigned int k=0; k < A.Nz(); k++) {
				for(unsigned int l=0; l < A.Nw(); l++) {
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

