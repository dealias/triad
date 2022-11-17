/* Array2h.h:  A high-performance multi-dimensional Hermitian C++ array class
   Copyright (C) 2022 John C. Bowman

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

#ifndef __Array2h_h__
#define __Array2h_h__ 1

#define __ARRAY2h_H_VERSION__ 1.00

// Defining NDEBUG improves optimization.

#include "Array.h"

namespace Array {

#undef __check

#ifdef NDEBUG
#define __check(i,n,o,dim,m)
#define __checkActivate(i,align) this->Activate(align)
#else
#define __check(i,n,o,dim,m) this->Check(i-(o),n,dim,m,o)
#define __checkActivate(i,align) this->CheckActivate(i,align)
#endif

// 2*(nx+1)*ny Array corresponding to compact 2D Hermitian "half
// plane", excluding origin.
template<class T>
class array2h : public array1<T> {
protected:
  unsigned int nx;
  unsigned int ny;
public:
  using array1<T>::Dimension;

  void Dimension(unsigned int nx0, unsigned int ny0) {
    nx=nx0; ny=ny0;
    this->size=2*nx*ny-nx-ny;
  }
  void Dimension(unsigned int nx0, unsigned int ny0, T *v0) {
    Dimension(nx0,ny0);
    this->v=v0;
    this->clear(this->allocated);
  }
  void Dimension(const array1<T> &A) {ArrayExit("Operation not implemented");}

  void Allocate(unsigned int nx0, unsigned int ny0, size_t align=0) {
    Dimension(nx0,ny0);
    __checkActivate(2,align);
  }

  array2h() : nx(0), ny(0) {}
  array2h(unsigned int nx0, unsigned int ny0, size_t align=0) {
    Allocate(nx0,ny0,align);
  }
  array2h(unsigned int nx0, unsigned int ny0, T *v0) {Dimension(nx0,ny0,v0);}

  unsigned int Nx() const {return 2*nx-1;}
  unsigned int Ny() const {return ny;}

#ifndef __NOARRAY2OPT
  T *operator [] (int ix) const {
    return ix <= 0 ? this->v+(nx+ix)*(ny-1)-ny : this->v+(nx+ix)*ny-nx-ny;
  }
#else
  array1<T> operator [] (int ix) const {
    __check(ix,Nx(),-nx+1,2,1);
    return array1<T>(ny,ix <= 0 ?
                     this->v+(nx+ix)*(ny-1)-ny :
                     this->v+(nx+ix)*ny-nx-ny);
  }
#endif
  T& operator () (int ix, int iy) const {
    __check(ix,Nx(),-nx+1,2,1);
    __check(iy,ny,ix <= 0,2,2);
    return this->v[ix <= 0 ? (nx+ix)*(ny-1)-ny+iy : (nx+ix)*ny-nx-ny+iy];
  }
  T& operator () (int i) const {
    __check(i,this->size,0,2,0);
    return this->v[i];
  }
  T* operator () () const {return this->v;}

  array2h<T>& operator = (T a) {this->Load(a); return *this;}
  array2h<T>& operator = (T *a) {this->Load(a); return *this;}
  array2h<T>& operator = (const array2h<T>& A) {
    __checkEqual(Nx(),A.Nx(),2,1);
    __checkEqual(ny,A.Ny(),2,2);
    this->Load(A());
    A.Purge();
    return *this;
  }
};

#undef __check
#undef __checkActivate

}

#endif
