/* Array.h:  A high-performance multi-dimensional C++ array class
Copyright (C) 2001-2004 John C. Bowman (bowman@math.ualberta.ca)

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

#define __ARRAY_H_VERSION__ 1.29

// Defining NDEBUG improves optimization but disables argument checking.
// Defining __NOARRAY2OPT inhibits special optimization of Array2[].

#include <iostream>
#include <sstream>
#include <unistd.h>
#include <climits>
#include <cstdlib>

#ifdef NDEBUG
#define __check(i,n,dim,m)
#define __checkSize()
#define __checkActivate(i) CheckActivate(i)
#define __NULLARRAY NULL;
#else
#define __check(i,n,dim,m) Check(i,n,dim,m)
#define __checkSize() CheckSize()
#define __checkActivate(i) Activate()
#define __NULLARRAY (void *) 0;
#ifndef __NOARRAY2OPT
#define __NOARRAY2OPT
#endif
#endif

using std::istream;
using std::ostream;
using std::ostringstream;
using std::cin;
using std::cout;
using std::endl;
using std::flush;

namespace Array {
  
inline ostream& _newl(ostream& s) {s << '\n'; return s;}

inline void ArrayExit(const char *x);
  
#ifndef __ExternalArrayExit
inline void ArrayExit(const char *x)
{
  cout << _newl << "ERROR: " << x << "." << endl;
  exit(1);
} 
#endif

template<class T>
class array1 {
 protected:
  T *v;
  unsigned int size;
  mutable int state;
 public:
  enum alloc_state {unallocated=0, allocated=1, temporary=2};
  virtual unsigned int Size() const {return size;}
  unsigned int CheckSize() const {
    if(!test(allocated) && size == 0)
      ArrayExit("Operation attempted on unallocated array"); 
    return size;
  }
	
  int test(int flag) const {return state & flag;}
  void clear(int flag) const {state &= ~flag;}
  void set(int flag) const {state |= flag;}
  void Activate() {
    v=new T[size];
    set(allocated);
  }
  void CheckActivate(int dim) {
    if (test(allocated)) {
      ostringstream buf;
      buf << "Reallocation of Array" << dim
	  << " attempted (must Deallocate first)";
      ArrayExit(buf.str().c_str());
    }
    Activate();
  }
  void Deallocate() const {
    if(test(allocated)) {delete [] v; clear(allocated);}
  }
  void Dimension(unsigned int nx0) {size=nx0;}
  void Dimension(unsigned int nx0, T *v0) {
    Dimension(nx0); v=v0; clear(allocated);
  }
  void Dimension(const array1<T>& A) {
    Dimension(A.size,A.v); state=A.test(temporary);
  }

  void Allocate(unsigned int nx0) {
    Dimension(nx0);
    __checkActivate(1);
  }
  
  array1() : size(0), state(unallocated) {}
  array1(const void *) : size(0), state(unallocated) {}
  array1(unsigned int nx0) : state(unallocated) {Allocate(nx0);}
  array1(unsigned int nx0, T *v0) : state(unallocated) {Dimension(nx0,v0);}
  array1(T *v0) : state(unallocated) {Dimension(INT_MAX,v0);}
  array1(const array1<T>& A) : v(A.v), size(A.size),
    state(A.test(temporary)) {}

  virtual ~array1() {Deallocate();}
	
  void Freeze() {state=unallocated;}
  void Hold() {if(test(allocated)) {state=temporary;}}
  void Purge() const {if(test(temporary)) {Deallocate(); state=unallocated;}}
	
  virtual void Check(int i, int n, unsigned int dim, unsigned int m,
		     int o=0) const {
    if(i < 0 || i >= n) {
      ostringstream buf;
      buf << "Array" << dim << " index ";
      if(m) buf << m << " ";
      buf << "is out of bounds (" << i+o;
      if(i < 0) buf << " < " << o;
      else buf << " > " << n+o-1;
      buf << ")";
      ArrayExit(buf.str().c_str());
    }
  }
	
  unsigned int Nx() const {return size;}
  
#ifdef NDEBUG
  typedef T *opt;
#else
  typedef array1<T> opt;
#endif
  
  T& operator [] (int ix) const {__check(ix,size,1,1); return v[ix];}
  T& operator () (int ix) const {__check(ix,size,1,1); return v[ix];}
  T* operator () () const {return v;}
  operator T* () const {return v;}
	
  array1<T> operator + (int i) const {return array1<T>(size-i,v+i);}
	
  void Load(T a) const {
    __checkSize();
    for(unsigned int i=0; i < size; i++) v[i]=a;
  }
  void Load(const T *a) const {memcpy(v,a,sizeof(T)*size);}
  void Store(T *a) const {memcpy(a,v,sizeof(T)*size);}
  void Set(T *a) {v=a; clear(allocated);}

  istream& Input (istream &s) const {
    __checkSize();
    for(unsigned int i=0; i < size; i++) s >> v[i];
    return s;
  }
	
  array1<T>& operator = (T a) {Load(a); return *this;}
  array1<T>& operator = (const T *a) {Load(a); return *this;}
  array1<T>& operator = (const array1<T>& A) {
    Load(A());
    A.Purge();
    return *this;
  }
	
  array1<T>& operator += (const array1<T>& A) {
    __checkSize();
    for(unsigned int i=0; i < size; i++) v[i] += A(i);
    return *this;
  }
  array1<T>& operator -= (const array1<T>& A) {
    __checkSize();
    for(unsigned int i=0; i < size; i++) v[i] -= A(i);
    return *this;
  }
  array1<T>& operator *= (const array1<T>& A) {
    __checkSize();
    for(unsigned int i=0; i < size; i++) v[i] *= A(i);
    return *this;
  }
  array1<T>& operator /= (const array1<T>& A) {
    __checkSize();
    for(unsigned int i=0; i < size; i++) v[i] /= A(i);
    return *this;
  }
	
  array1<T>& operator += (T a) {
    __checkSize();
    for(unsigned int i=0; i < size; i++) v[i] += a;
    return *this;
  }
  array1<T>& operator -= (T a) {
    __checkSize();
    for(unsigned int i=0; i < size; i++) v[i] -= a;
    return *this;
  }
  array1<T>& operator *= (T a) {
    __checkSize();
    for(unsigned int i=0; i < size; i++) v[i] *= a;
    return *this;
  }
  array1<T>& operator /= (T a) {
    __checkSize();
    T ainv=1.0/a;
    for(unsigned int i=0; i < size; i++) v[i] *= ainv;
    return *this;
  }
	
  double L1() const {
    __checkSize();
    double norm=0.0;
    for(unsigned int i=0; i < size; i++) norm += abs(v[i]);
    return norm/size;
  }
#ifdef __ArrayExtensions
  double Abs2() const {
    __checkSize();
    double norm=0.0;
    for(unsigned int i=0; i < size; i++) norm += abs2(v[i]);
    return norm;
  }
  double L2() const {
    return sqrt(Abs2()/size);
  }
  double LInfinity() const {
    __checkSize();
    double norm=0.0;
    for(unsigned int i=0; i < size; i++) norm=max(norm,abs(v[i]));
    return norm;
  }
  double LMinusInfinity() const {
    __checkSize();
    double norm=DBL_MAX;
    for(unsigned int i=0; i < size; i++) norm=min(norm,abs(v[i]));
    return norm;
  }
#endif	
};

template<class T>
void swaparray(T& A, T& B) {
  T C;
  C.Dimension(A);
  A.Dimension(B);
  B.Dimension(C);
}
  
template<class T>
void leftshiftarray(T& A, T& B, T& C) {
  T D;
  D.Dimension(A);
  A.Dimension(B);
  B.Dimension(C);
  C.Dimension(D);
}
  
template<class T>
ostream& operator << (ostream& s, const array1<T>& A)
{
  T *p=A();
  for(unsigned int i=0; i < A.Nx(); i++) {
    s << *(p++) << " ";
  }
  return s;
}

template<class T>
istream& operator >> (istream& s, const array1<T>& A)
{
  return A.Input(s);
}

template<class T>
class array2 : public array1<T> {
 protected:
  unsigned int nx;
  unsigned int ny;
 public:
  void Dimension(unsigned int nx0, unsigned int ny0) {
    nx=nx0; ny=ny0;
    size=nx*ny;
  }
  void Dimension(unsigned int nx0, unsigned int ny0, T *v0) {
    Dimension(nx0,ny0);
    v=v0;
    clear(allocated);
  }
  
  void Dimension(const array1<T> &A) {ArrayExit("Operation not implemented");} 
  
  void Allocate(unsigned int nx0, unsigned int ny0) {
    Dimension(nx0,ny0);
    __checkActivate(2);
  }
	
  array2() : nx(0), ny(0) {}
  array2(unsigned int nx0, unsigned int ny0) {Allocate(nx0,ny0);}
  array2(unsigned int nx0, unsigned int ny0, T *v0) {Dimension(nx0,ny0,v0);}
	
  unsigned int Nx() const {return nx;}
  unsigned int Ny() const {return ny;}

#ifndef __NOARRAY2OPT
  T *operator [] (int ix) const {
    return v+ix*ny;
  }
#else
  array1<T> operator [] (int ix) const {
    __check(ix,nx,2,1);
    return array1<T>(ny,v+ix*ny);
  }
#endif
  T& operator () (int ix, int iy) const {
    __check(ix,nx,2,1);
    __check(iy,ny,2,2);
    return v[ix*ny+iy];
  }
  T& operator () (int i) const {
    __check(i,size,2,0);
    return v[i];
  }
  T* operator () () const {return v;}
	
  array2<T>& operator = (T a) {Load(a); return *this;}
  array2<T>& operator = (T *a) {Load(a); return *this;}
  array2<T>& operator = (const array2<T>& A) {
    Load(A());
    A.Purge();
    return *this;
  }
	
  array2<T>& operator += (const array2<T>& A) {
    __checkSize();
    for(unsigned int i=0; i < size; i++) v[i] += A(i);
    return *this;
  }
  array2<T>& operator -= (const array2<T>& A) {
    __checkSize();
    for(unsigned int i=0; i < size; i++) v[i] -= A(i);
    return *this;
  }
  array2<T>& operator *= (const array2<T>& A);
	
  array2<T>& operator += (T a) {
    __checkSize();
    unsigned int inc=ny+1;
    for(unsigned int i=0; i < size; i += inc) v[i] += a;
    return *this;
  }
  array2<T>& operator -= (T a) {
    __checkSize();
    unsigned int inc=ny+1;
    for(unsigned int i=0; i < size; i += inc) v[i] -= a;
    return *this;
  }
  array2<T>& operator *= (T a) {
    __checkSize();
    for(unsigned int i=0; i < size; i++) v[i] *= a;
    return *this;
  }
};

template<class T>
ostream& operator << (ostream& s, const array2<T>& A)
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
istream& operator >> (istream& s, const array2<T>& A)
{
  return A.Input(s);
}

template<class T>
class array3 : public array1<T> {
 protected:
  unsigned int nx;
  unsigned int ny;
  unsigned int nz;
  unsigned int nyz;
 public:
  void Allocate(unsigned int nx0, unsigned int ny0, unsigned int nz0) {
    Dimension(nx0,ny0,nz0);
    __checkActivate(3);
  }
  void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0) {
    nx=nx0; ny=ny0; nz=nz0; nyz=ny*nz;
    size=nx*nyz;
  }
  void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0, T *v0) {
    Dimension(nx0,ny0,nz0);
    v=v0;
    clear(allocated);
  }
	
  array3() : nx(0), ny(0), nz(0), nyz(0) {}
  array3(unsigned int nx0, unsigned int ny0, unsigned int nz0) {
    Allocate(nx0,ny0,nz0);
  }
  array3(unsigned int nx0, unsigned int ny0, unsigned int nz0, T *v0) {
    Dimension(nx0,ny0,nz0,v0);
  }
	
  unsigned int Nx() const {return nx;}
  unsigned int Ny() const {return ny;}
  unsigned int Nz() const {return nz;}

  array2<T> operator [] (int ix) const {
    __check(ix,nx,3,1);
    return array2<T>(ny,nz,v+ix*nyz);
  }
  T& operator () (int ix, int iy, int iz) const {
    __check(ix,nx,3,1);
    __check(iy,ny,3,2);
    __check(iz,nz,3,3);
    return v[ix*nyz+iy*nz+iz];
  }
  T& operator () (int i) const {
    __check(i,size,3,0);
    return v[i];
  }
  T* operator () () const {return v;}
	
  array3<T>& operator = (T a) {Load(a); return *this;}
  array3<T>& operator = (T *a) {Load(a); return *this;}
  array3<T>& operator = (const array3<T>& A) {
    Load(A());
    A.Purge(); 
    return *this;
  }
	
  array3<T>& operator += (array3<T>& A) {
    __checkSize();
    for(unsigned int i=0; i < size; i++) v[i] += A(i);
    return *this;
  }
  array3<T>& operator -= (array3<T>& A) {
    __checkSize();
    for(unsigned int i=0; i < size; i++) v[i] -= A(i);
    return *this;
  }
	
  array3<T>& operator += (T a) {
    __checkSize();
    unsigned int inc=nyz+nz+1;
    for(unsigned int i=0; i < size; i += inc) v[i] += a;
    return *this;
  }
  array3<T>& operator -= (T a) {
    __checkSize();
    unsigned int inc=nyz+nz+1;
    for(unsigned int i=0; i < size; i += inc) v[i] -= a;
    return *this;
  }
};

template<class T>
ostream& operator << (ostream& s, const array3<T>& A)
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
istream& operator >> (istream& s, const array3<T>& A)
{
  return A.Input(s);
}

template<class T>
class array4 : public array1<T> {
 protected:
  unsigned int nx;
  unsigned int ny;
  unsigned int nz;
  unsigned int nw;
  unsigned int nyz;
  unsigned int nzw;
  unsigned int nyzw;
 public:
  void Allocate(unsigned int nx0, unsigned int ny0, unsigned int nz0,
		unsigned int nw0) {
    Dimension(nx0,ny0,nz0,nw0);
    __checkActivate(4);
  }
  void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
		 unsigned int nw0) {
    nx=nx0; ny=ny0; nz=nz0; nw=nw0; nzw=nz*nw; nyzw=ny*nzw;
    size=nx*nyzw;
  }
  void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
		 unsigned int nw0, T *v0) {
    Dimension(nx0,ny0,nz0,nw0);
    v=v0;
    clear(allocated);
  }
	
  array4() : nx(0), ny(0), nz(0), nw(0), nyz(0), nzw(0), nyzw(0) {}
  array4(unsigned int nx0, unsigned int ny0, unsigned int nz0,
	 unsigned int nw0) {Allocate(nx0,ny0,nz0,nw0);}
  array4(unsigned int nx0, unsigned int ny0, unsigned int nz0,
	 unsigned int nw0, T *v0) {
    Dimension(nx0,ny0,nz0,nw0,v0);
  }

  unsigned int Nx() const {return nx;}
  unsigned int Ny() const {return ny;}
  unsigned int Nz() const {return ny;}
  unsigned int N4() const {return nw;}

  array3<T> operator [] (int ix) const {
    __check(ix,nx,3,1);
    return array3<T>(ny,nz,nw,v+ix*nyzw);
  }
  T& operator () (int ix, int iy, int iz, int iw) const {
    __check(ix,nx,4,1);
    __check(iy,ny,4,2);
    __check(iz,nz,4,3);
    __check(iw,nw,4,4);
    return v[ix*nyzw+iy*nzw+iz*nw+iw];
  }
  T& operator () (int i) const {
    __check(i,size,4,0);
    return v[i];
  }
  T* operator () () const {return v;}
	
  array4<T>& operator = (T a) {Load(a); return *this;}
  array4<T>& operator = (T *a) {Load(a); return *this;}
  array4<T>& operator = (const array4<T>& A) {
    Load(A());
    A.Purge();
    return *this;
  }
	
  array4<T>& operator += (array4<T>& A) {
    __checkSize();
    for(unsigned int i=0; i < size; i++) v[i] += A(i);
    return *this;
  }
  array4<T>& operator -= (array4<T>& A) {
    __checkSize();
    for(unsigned int i=0; i < size; i++) v[i] -= A(i);
    return *this;
  }
	
  array4<T>& operator += (T a) {
    __checkSize();
    unsigned int inc=nyzw+nzw+nw+1;
    for(unsigned int i=0; i < size; i += inc) v[i] += a;
    return *this;
  }
  array4<T>& operator -= (T a) {
    __checkSize();
    unsigned int inc=nyzw+nzw+nw+1;
    for(unsigned int i=0; i < size; i += inc) v[i] -= a;
    return *this;
  }
};

template<class T>
ostream& operator << (ostream& s, const array4<T>& A)
{
  T *p=A;
  for(unsigned int i=0; i < A.Nx(); i++) {
    for(unsigned int j=0; j < A.Ny(); j++) {
      for(unsigned int k=0; k < A.Nz(); k++) {
	for(unsigned int l=0; l < A.N4(); l++) {
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
istream& operator >> (istream& s, const array4<T>& A)
{
  return A.Input(s);
}

template<class T>
class array5 : public array1<T> {
 protected:
  unsigned int nx;
  unsigned int ny;
  unsigned int nz;
  unsigned int nw;
  unsigned int nv;
  unsigned int nwv;
  unsigned int nzwv;
  unsigned int nyzwv;
 public:
  void Allocate(unsigned int nx0, unsigned int ny0, unsigned int nz0,
		unsigned int nw0, unsigned int nv0) {
    Dimension(nx0,ny0,nz0,nw0,nv0);
    __checkActivate(5);
  }
  void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
		 unsigned int nw0, unsigned int nv0) {
    nx=nx0; ny=ny0; nz=nz0; nw=nw0; nv=nv0; nwv=nw*nv; nzwv=nz*nwv;
    nyzwv=ny*nzwv;
    size=nx*nyzwv;
  }
  void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
		 unsigned int nw0, unsigned int nv0, T *v0) {
    Dimension(nx0,ny0,nz0,nw0,nv0);
    v=v0;
    clear(allocated);
  }
	
  array5() : nx(0), ny(0), nz(0), nw(0), nv(0), nwv(0), nzwv(0), nyzwv(0) {}
  array5(unsigned int nx0, unsigned int ny0, unsigned int nz0,
	 unsigned int nw0, unsigned int nv0) {Allocate(nx0,ny0,nz0,nw0,nv0);}
  array5(unsigned int nx0, unsigned int ny0, unsigned int nz0,
	 unsigned int nw0, unsigned int nv0, T *v0) {
    Dimension(nx0,ny0,nz0,nw0,nv0,nv0);
  }

  unsigned int Nx() const {return nx;}
  unsigned int Ny() const {return ny;}
  unsigned int Nz() const {return ny;}
  unsigned int N4() const {return nw;}
  unsigned int N5() const {return nv;}

  array4<T> operator [] (int ix) const {
    __check(ix,nx,4,1);
    return array4<T>(ny,nz,nw,nv,v+ix*nyzwv);
  }
  T& operator () (int ix, int iy, int iz, int iw, int iv) const {
    __check(ix,nx,5,1);
    __check(iy,ny,5,2);
    __check(iz,nz,5,3);
    __check(iw,nw,5,4);
    __check(iv,nv,5,5);
    return v[ix*nyzwv+iy*nzwv+iz*nwv+iw*nv+iv];
  }
  T& operator () (int i) const {
    __check(i,size,5,0);
    return v[i];
  }
  T* operator () () const {return v;}
	
  array5<T>& operator = (T a) {Load(a); return *this;}
  array5<T>& operator = (T *a) {Load(a); return *this;}
  array5<T>& operator = (const array5<T>& A) {
    Load(A());
    A.Purge();
    return *this;
  }
	
  array5<T>& operator += (array5<T>& A) {
    __checkSize();
    for(unsigned int i=0; i < size; i++) v[i] += A(i);
    return *this;
  }
  array5<T>& operator -= (array5<T>& A) {
    __checkSize();
    for(unsigned int i=0; i < size; i++) v[i] -= A(i);
    return *this;
  }
	
  array5<T>& operator += (T a) {
    __checkSize();
    unsigned int inc=nyzwv+nzwv+nwv+nv+1;
    for(unsigned int i=0; i < size; i += inc) v[i] += a;
    return *this;
  }
  array5<T>& operator -= (T a) {
    __checkSize();
    unsigned int inc=nyzwv+nzwv+nwv+nv+1;
    for(unsigned int i=0; i < size; i += inc) v[i] -= a;
    return *this;
  }
};

template<class T>
ostream& operator << (ostream& s, const array5<T>& A)
{
  T *p=A;
  for(unsigned int i=0; i < A.Nx(); i++) {
    for(unsigned int j=0; j < A.Ny(); j++) {
      for(unsigned int k=0; k < A.Nz(); k++) {
	for(unsigned int l=0; l < A.N4(); l++) {
	  for(unsigned int l=0; l < A.N5(); l++) {
	    s << *(p++) << " ";
	  }
	  s << _newl;
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
istream& operator >> (istream& s, const array5<T>& A)
{
  return A.Input(s);
}

#undef __check

#ifdef NDEBUG
#define __check(i,n,o,dim,m)
#else
#define __check(i,n,o,dim,m) Check(i-o,n,dim,m,o)
#endif

template<class T>
class Array1 : public array1<T> {
 protected:
  T *voff; // Offset pointer to memory block
  int ox;
 public:
  void Offsets() {
    voff=v-ox;
  }
  void Dimension(unsigned int nx0, int ox0=0) {
    size=nx0;
    ox=ox0;
  }
  void Dimension(unsigned int nx0, T *v0, int ox0=0) {
    Dimension(nx0,ox0);
    v=v0;
    Offsets();
    clear(allocated);
  }
  void Dimension(const Array1<T>& A) {
    Dimension(A.size,A.v,A.ox); state=A.test(temporary);
  }
  void Allocate(unsigned int nx0, int ox0=0) {
    Dimension(nx0,ox0);
    __checkActivate(1);
    Offsets();
  }
	
  Array1() : ox(0) {}
  Array1(unsigned int nx0, int ox0=0) {
    Allocate(nx0,ox0);
  }
  Array1(unsigned int nx0, T *v0, int ox0=0) {
    Dimension(nx0,v0,ox0);
  }
  Array1(T *v0) {
    Dimension(INT_MAX,v0);
  }

#ifdef NDEBUG
  typedef T *opt;
#else
  typedef Array1<T> opt;
#endif  
  
  T& operator [] (int ix) const {__check(ix,size,ox,1,1); return voff[ix];}
  T& operator () (int i) const {__check(i,size,0,1,1); return v[i];}
  T* operator () () const {return voff;}
  operator T* () const {return voff;}
	
  Array1<T> operator + (int i) const {return Array1<T>(size-i,v+i,ox);}
  void Set(T *a) {v=a; Offsets(); clear(allocated);}
	
  Array1<T>& operator = (T a) {Load(a); return *this;}
  Array1<T>& operator = (const T *a) {Load(a); return *this;}
  Array1<T>& operator = (const Array1<T>& A) {
    Load(A());
    A.Purge();
    return *this;
  }
  Array1<T>& operator = (const array1<T>& A) {
    Load(A());
    A.Purge();
    return *this;
  }
  
  int Ox() const {return ox;}
};

template<class T>
class Array2 : public array2<T> {
 protected:
  T *voff,*vtemp;
  int ox,oy;
 public:
  void Offsets() {
    vtemp=v-ox*(int) ny;
    voff=vtemp-oy;
  }
  void Dimension(unsigned int nx0, unsigned int ny0, int ox0=0, int oy0=0) {
    nx=nx0; ny=ny0;
    size=nx*ny;
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
    __checkActivate(2);
    Offsets();
  }

  Array2() : ox(0), oy(0) {}
  Array2(unsigned int nx0, unsigned int ny0, int ox0=0, int oy0=0) {
    Allocate(nx0,ny0,ox0,oy0);
  }
  Array2(unsigned int nx0, unsigned int ny0, T *v0, int ox0=0, int oy0=0) {
    Dimension(nx0,ny0,v0,ox0,oy0);
  }

#ifndef __NOARRAY2OPT
  T *operator [] (int ix) const {
    return voff+ix*(int) ny;
  }
#else
  Array1<T> operator [] (int ix) const {
    __check(ix,nx,ox,2,1);
    return Array1<T>(ny,vtemp+ix*(int) ny,oy);
  }
#endif
  
  T& operator () (int ix, int iy) const {
    __check(ix,nx,ox,2,1);
    __check(iy,ny,oy,2,2);
    return voff[ix*(int) ny+iy];
  }
  T& operator () (int i) const {
    __check(i,size,0,2,0);
    return v[i];
  }
  T* operator () () const {return voff;}
  void Set(T *a) {v=a; Offsets(); clear(allocated);}
	
  Array2<T>& operator = (T a) {Load(a); return *this;}
  Array2<T>& operator = (T *a) {Load(a); return *this;}
  Array2<T>& operator = (const Array2<T>& A) {
    Load(A());
    A.Purge();
    return *this;
  }
  Array2<T>& operator = (const array2<T>& A) {
    Load(A());
    A.Purge();
    return *this;
  }
	
  Array2<T>& operator *= (const Array2<T>& A);
	
  Array2<T>& operator *= (T a) {
    __checkSize();
    for(unsigned int i=0; i < size; i++) v[i] *= a;
    return *this;
  }
	
  int Ox() const {return ox;}
  int Oy() const {return oy;}
  
  void Identity() {
    Load((T) 0.0);
    __checkSize();
    unsigned int inc=ny+1;
    for(unsigned int i=0; i < size; i += inc) v[i]=1.0;
  }
};

template<class T>
class Array3 : public array3<T> {
 protected:
  T *voff,*vtemp;
  int ox,oy,oz;
 public:
  void Offsets() {
    vtemp=v-ox*(int) nyz;
    voff=vtemp-oy*(int) nz-oz;
  }
  void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
		 int ox0=0, int oy0=0, int oz0=0) {
    nx=nx0; ny=ny0; nz=nz0; nyz=ny*nz;
    size=nx*nyz;
    ox=ox0; oy=oy0; oz=oz0;
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
    __checkActivate(3);
    Offsets();
  }
	
  Array3() : ox(0), oy(0), oz(0) {}
  Array3(unsigned int nx0, unsigned int ny0, unsigned int nz0,
	 int ox0=0, int oy0=0, int oz0=0) {
    Allocate(nx0,ny0,nz0,ox0,oy0,oz0);
  }
  Array3(unsigned int nx0, unsigned int ny0, unsigned int nz0, T *v0,
	 int ox0=0, int oy0=0, int oz0=0) {
    Dimension(nx0,ny0,nz0,v0,ox0,oy0,oz0);
  }
	
  Array2<T> operator [] (int ix) const {
    __check(ix,nx,ox,3,1);
    return Array2<T>(ny,nz,vtemp+ix*(int) nyz,oy,oz);
  }
  T& operator () (int ix, int iy, int iz) const {
    __check(ix,nx,ox,3,1);
    __check(iy,ny,oy,3,2);
    __check(iz,nz,oz,3,3);
    return voff[ix*(int) nyz+iy*(int) nz+iz];
  }
  T& operator () (int i) const {
    __check(i,size,0,3,0);
    return v[i];
  }
  T* operator () () const {return voff;}
  void Set(T *a) {v=a; Offsets(); clear(allocated);}
	
  Array3<T>& operator = (T a) {Load(a); return *this;}
  Array3<T>& operator = (T *a) {Load(a); return *this;}
  Array3<T>& operator = (const Array3<T>& A) {
    Load(A());
    A.Purge(); 
    return *this;
  }
  Array3<T>& operator = (const array3<T>& A) {
    Load(A());
    A.Purge(); 
    return *this;
  }
  
  int Ox() const {return ox;}
  int Oy() const {return oy;}
  int Oz() const {return oz;}

};

template<class T>
class Array4 : public array4<T> {
 protected:
  T *voff,*vtemp;
  int ox,oy,oz,ow;
 public:
  void Offsets() {
    vtemp=v-ox*(int) nyzw;
    voff=vtemp-oy*(int) nzw-oz*(int) nw-ow;
  }
  void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
		 unsigned int nw0, 
		 int ox0=0, int oy0=0, int oz0=0, int ow0=0) {
    nx=nx0; ny=ny0; nz=nz0; nw=nw0; nzw=nz*nw; nyzw=ny*nzw;
    size=nx*nyzw;
    ox=ox0; oy=oy0; oz=oz0; ow=ow0;
  }
  void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
		 unsigned int nw0, T *v0,
		 int ox0=0, int oy0=0, int oz0=0, int ow0=0) {
    Dimension(nx0,ny0,nz0,nw0,ox0,oy0,oz0,ow0);
    v=v0;
    Offsets();
    clear(allocated);
  }
  void Allocate(unsigned int nx0, unsigned int ny0, unsigned int nz0,
		unsigned int nw0,
		int ox0=0, int oy0=0, int oz0=0, int ow0=0) {
    Dimension(nx0,ny0,nz0,nw0,ox0,oy0,oz0,ow0);
    __checkActivate(4); 
    Offsets();
  }
	
  Array4() : ox(0), oy(0), oz(0), ow(0) {}
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

  Array3<T> operator [] (int ix) const {
    __check(ix,nx,ox,3,1);
    return Array3<T>(ny,nz,nw,vtemp+ix*(int) nyzw,oy,oz,ow);
  }
  T& operator () (int ix, int iy, int iz, int iw) const {
    __check(ix,nx,ox,4,1);
    __check(iy,ny,oy,4,2);
    __check(iz,nz,oz,4,3);
    __check(iw,nw,ow,4,4);
    return voff[ix*(int) nyzw+iy*(int) nzw+iz*(int) nw+iw];
  }
  T& operator () (int i) const {
    __check(i,size,0,4,0);
    return v[i];
  }
  T* operator () () const {return voff;}
  void Set(T *a) {v=a; Offsets(); clear(allocated);}
	
  Array4<T>& operator = (T a) {Load(a); return *this;}
  Array4<T>& operator = (T *a) {Load(a); return *this;}
  Array4<T>& operator = (const Array4<T>& A) {
    Load(A());
    A.Purge();
    return *this;
  }
  Array4<T>& operator = (const array4<T>& A) {
    Load(A());
    A.Purge();
    return *this;
  }
  
  int Ox() const {return ox;}
  int Oy() const {return oy;}
  int Oz() const {return oz;}
  int O4() const {return ow;}
};

template<class T>
class Array5 : public array5<T> {
 protected:
  T *voff,*vtemp;
  int ox,oy,oz,ow,ov;
 public:
  void Offsets() {
    vtemp=v-ox*(int) nyzwv;
    voff=vtemp-oy*(int) nzwv-oz*(int) nwv-ow*(int) nv-ov;
  }
  void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
		 unsigned int nw0,  unsigned int nv0,
		 int ox0=0, int oy0=0, int oz0=0, int ow0=0, int ov0=0) {
    nx=nx0; ny=ny0; nz=nz0; nw=nw0; nv=nv0; nwv=nw*nv; nzwv=nz*nwv;
    nyzwv=ny*nzwv;
    size=nx*nyzwv;
    ox=ox0; oy=oy0; oz=oz0; ow=ow0; ov=ov0;
  }
  void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
		 unsigned int nw0, unsigned int nv0, T *v0,
		 int ox0=0, int oy0=0, int oz0=0, int ow0=0, int ov0=0) {
    Dimension(nx0,ny0,nz0,nw0,nv0,ox0,oy0,oz0,ow0,ov0);
    v=v0;
    Offsets();
    clear(allocated);
  }
  void Allocate(unsigned int nx0, unsigned int ny0, unsigned int nz0,
		unsigned int nw0, unsigned int nv0,
		int ox0=0, int oy0=0, int oz0=0, int ow0=0, int ov0=0) {
    Dimension(nx0,ny0,nz0,nw0,nv0,ox0,oy0,oz0,ow0,ov0);
    __checkActivate(5); 
    Offsets();
  }
	
  Array5() : ox(0), oy(0), oz(0), ow(0), ov(0) {}
  Array5(unsigned int nx0, unsigned int ny0, unsigned int nz0,
	 unsigned int nw0, unsigned int nv0,
	 int ox0=0, int oy0=0, int oz0=0, int ow0=0, int ov0=0) {
    Allocate(nx0,ny0,nz0,nw0,nv0,ox0,oy0,oz0,ow0,ov0);
  }
  Array5(unsigned int nx0, unsigned int ny0, unsigned int nz0,
	 unsigned int nw0, unsigned int nv0, T *v0,
	 int ox0=0, int oy0=0, int oz0=0, int ow0=0, int ov0=0) {
    Dimension(nx0,ny0,nz0,nw0,nv0,v0,ox0,oy0,oz0,ow0,ov0);
  }

  Array4<T> operator [] (int ix) const {
    __check(ix,nx,ox,4,1);
    return Array4<T>(ny,nz,nw,nv,vtemp+ix*(int) nyzwv,oy,oz,ow,ov);
  }
  T& operator () (int ix, int iy, int iz, int iw, int iv) const {
    __check(ix,nx,ox,5,1);
    __check(iy,ny,oy,5,2);
    __check(iz,nz,oz,5,3);
    __check(iw,nw,ow,5,4);
    __check(iv,nv,ov,5,5);
    return voff[ix*(int) nyzwv+iy*(int) nzwv+iz*(int) nwv+iw*(int) nv+iv];
  }
  T& operator () (int i) const {
    __check(i,size,0,5,0);
    return v[i];
  }
  T* operator () () const {return voff;}
  void Set(T *a) {v=a; Offsets(); clear(allocated);}
	
  Array5<T>& operator = (T a) {Load(a); return *this;}
  Array5<T>& operator = (T *a) {Load(a); return *this;}
  Array5<T>& operator = (const Array5<T>& A) {
    Load(A());
    A.Purge();
    return *this;
  }
  Array5<T>& operator = (const array5<T>& A) {
    Load(A());
    A.Purge();
    return *this;
  }
  int Ox() const {return ox;}
  int Oy() const {return oy;}
  int Oz() const {return oz;}
  int O4() const {return ow;}
  int O5() const {return ov;}
};

template<class T>
inline bool Active(array1<T>& A)
{
  return A.Size();
}

template<class T>
inline bool Active(T *A)
{
  return A;
}

template<class T>
inline void Set1(T *&A, T *v)
{
  A=v;
}

template<class T>
inline void Set1(array1<T>& A, T *v)
{
  A.Set(v);
}

template<class T>
inline void Set1(array1<T>& A, const array1<T>& B)
{
  A.Set(B());
}

template<class T>
inline void Dimension1(T *&A, unsigned int, T *v)
{
  A=v;
}

template<class T>
inline void Dimension1(array1<T>& A, unsigned int n, T *v)
{
  A.Dimension(n,v);
}

template<class T>
inline void Dimension1(T *&A, T *v)
{
  A=v;
}

template<class T>
inline void Dimension1(array1<T>& A, const array1<T>& B)
{
  A.Dimension(B);
}

template<class T>
inline void Dimension1(Array1<T>& A, unsigned int n, const array1<T>& B, int o)
{
  A.Dimension(n,B,o);
}

template<class T>
inline void Dimension1(Array1<T>& A, unsigned int n, T *v, int o)
{
  A.Dimension(n,v,o);
}

template<class T>
inline void Dimension1(T *&A, unsigned int, T *v, int o)
{
  A=v-o;
}

template<class T>
inline void Allocate1(T *&A, unsigned int n)
{
  A=new T[n];
}

template<class T>
inline void Allocate1(T *&A, unsigned int n, int o)
{
  A=new T[n]-o;
}

template<class T>
inline void Allocate1(array1<T>& A, unsigned int n)
{  
  A.Allocate(n);
}

template<class T>
inline void Allocate1(Array1<T>& A, unsigned int n, int o)
{  
  A.Allocate(n,o);
}

template<class T>
inline void Deallocate1(T *A)
{
  delete A;
}

template<class T>
inline void Deallocate1(T *A, int o)
{
  delete (A+o);
}

template<class T>
inline void Deallocate1(array1<T>& A)
{  
  A.Deallocate();
}

}

// Abbreviated form for optimized 1D arrays:
#define array1(T) array1<T>::opt
#define Array1(T) Array1<T>::opt

#undef __check
#undef __checkSize
#undef __checkActivate

#endif
