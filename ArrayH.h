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

#ifdef NDEBUG
#define __checkH(ix,iy)
#else
#define __checkH(ix,iy) this->CheckH(ix,iy)
#endif

#include "Array.h"

#ifdef NDEBUG
#define __check(i,n,dim,m)
#else
#define __check(i,n,dim,m) Check(i,n,dim,m)
#endif

namespace Array {
  
template<class T>
class Array1H : public array1<T> {
public:
  Array1H() {}
  Array1H(unsigned int mx) {this->Allocate(mx);}
  Array1H(unsigned int mx, T *v0) {this->Dimension(mx,v0);}
  
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
#define __check(i,n,o,dim,m) Check(i-o,n,dim,m,o)
#endif

template<class T>
class Array1HH {
  typename array1<T>::opt pos,neg;
public:
  Array1HH(unsigned int mx, T *vpos, T *vneg) {
    Dimension(pos,mx,vpos);
    Dimension(neg,mx,vneg);
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
  
  virtual void CheckH(int ix, int iy) const {
    //Check(ix,this->nx,2,1); // FIXME: segfaults
    if(abs(iy) >= this->ny) {
      std::ostringstream buf;
      buf << "Array2H index ";
      if(this->ny) buf << this->ny << " ";
      buf << "is out of bounds (" << iy;
      if(this->ny == 0) 
        buf << " index given to empty array";
      else {
        if(iy > 0) buf << ">" << this->ny;
        else buf << "<" << -(int) (this->ny-1);
      }
      buf << ")";
      ArrayExit(buf.str().c_str());
    }
  }

  Array2H() {}
  Array2H(unsigned int mx, unsigned int my) {
    this->Allocate(2*mx-1,my,1-mx,0);
  }
  Array2H(unsigned int mx, unsigned int my, T *v0) {
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
