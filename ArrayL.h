/* ArrayL.h:  A high-performance periodic multi-dimensional C++ array class
Copyright (C) 2001 John C. Bowman (bowman@math.ualberta.ca)

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

#ifndef __ArrayL_h__
#define __ArrayL_h__ 1

#define __ARRAYL_H_VERSION__ 1.22

// Lower triangular Array class

// Defining NDEBUG improves optimization.

#include "Array.h"

#ifdef NDEBUG
#define __check(i,n,dim,m)
#else
#define __check(i,n,dim,m) this->Check(i,n,dim,m)
#endif

namespace Array {
  
template<class T>
class array2L : public array1<T> {
  size_t n;
 public:
  void Dimension(size_t n0) {n=n0; this->size=n*(n+1)/2;}

  size_t index(int i, int j) const {return i*(i+1)/2+j;}
  size_t N() const {return n;}
  
  array2L() {}
  array2L(size_t n0, size_t align=0) {
    Dimension(n0);
    this->CheckActivate(align);
  }
  array2L(size_t n0, T *v0) {Dimension(n0); Dimension(n0,v0);}
  array2L(const array2L<T>& A) : n(A.n) {
    this->v=A.v; this->size=A.size; this->state=A.test(this->temporary);
  }
	
  array1<T> operator [] (size_t i) const {
    __check(i,n,1,1);
    return array1<T>(i+1,this->v+index(i,0));
  }
	
  array1<T> operator [] (int i) const {
    __check(i,n,1,1);
    return array1<T>(i+1,this->v+index(i,0));
  }
	
  T& operator () (size_t i, size_t j) const {
    __check(i,n,2,1);
    __check(j,i+1,2,1);
    return this->v[index(i,j)];
  }
  T* operator () () const {return this->v;}
	
  array2L<T>& operator = (T a) {Load(a); return *this;}
};
  
template<class T>
  std::ostream& operator << (std::ostream& s, const array2L<T>& A)
{
  T *p=A();
  for(size_t i=0; i < A.N(); i++) {
    for(size_t j=0; j <= i; j++) {
      s << *(p++) << " ";
    }
    s << _newl;
  }
  s << std::flush;
  return s;
}

}

#undef __check

#endif
