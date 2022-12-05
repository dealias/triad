/* ArrayH.h:  A high-performance multi-dimensional Hermitian C++ array class
   Copyright (C) 2011 John C. Bowman and Malcolm Roberts

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

#ifndef __ArrayH_h__
#define __ArrayH_h__ 1

#define __ARRAYH_H_VERSION__ 1

// Defining NDEBUG improves optimization.

#include "Array.h"

namespace Array {

template<class T>
class Array1H : public array1<T> {
public:
  Array1H() {}
  Array1H(size_t mx) {this->Allocate(mx);}
  Array1H(size_t mx, T *v0) {this->Dimension(mx,v0);}

  void set(int ix, T a) const {
    if(ix >= 0)
      array1<T>::operator()(ix)=a;
    else
      array1<T>::operator()(-ix)=conj(a);
  }

  T operator () (int ix) const {
    if(ix >= 0)
      return array1<T>::operator()(ix);
    else
      return conj(array1<T>::operator()(-ix));
  }
  T operator [] (int ix) const {
    return operator()(ix);
  }
  Array1H<T>& operator = (T a) {Load(a); return *this;}
};

#undef __check

#ifdef NDEBUG
#define __check(i,n,o,dim,m)
#else
#define __check(i,n,o,dim,m) this->Check(i-o,n,dim,m,o)
#endif

template<class T>
class Array1HH {
  typename array1<T>::opt pos,neg;
public:
  Array1HH(size_t mx, T *vpos, T *vneg) {
    Dimension(pos,mx,vpos);
    Dimension(neg,mx,vneg);
  }

  void set(int ix, T a) const {
    if(ix >= 0)
      pos(ix)=a;
    if(ix <= 0)
      neg(-ix)=conj(a);
  }

  T operator () (int ix) const {
    if(ix >= 0)
      return pos(ix);
    else
      return conj(neg(-ix));
  }
};

template<class T>
class Array2H : public Array2<T> {
public:
  Array2H() {}
  Array2H(size_t mx, size_t my) {
    this->Allocate(2*mx-1,my,1-mx,0);
  }
  Array2H(size_t mx, size_t my, T *v0) {
    this->Allocate(2*mx-1,my,1-mx,0,v0);
  }

  void set(int ix, int iy, T a) const {
    if(iy >= 0)
      Array2<T>::operator()(ix,iy)=a;
    if(iy <= 0)
      Array2<T>::operator()(-ix,-iy)=conj(a);
  }

  T operator () (int ix, int iy) const {
    if(iy >= 0)
      return Array2<T>::operator()(ix,iy);
    else
      return conj(Array2<T>::operator()(-ix,-iy));
  }

  Array1HH<T> operator [] (int ix) const {
    __check(ix,this->nx,this->ox,2,1);
    return Array1HH<T>(this->ny,
                       this->vtemp+ix*(int) this->ny,
                       this->vtemp-ix*(int) this->ny);
  }

  Array2H<T>& operator = (T a) {Load(a); return *this;}
};

}

#endif
