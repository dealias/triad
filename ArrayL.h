/* Array2n.h:  A high-performance periodic multi-dimensional C++ array class
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

#ifndef __ArrayL_h__
#define __ArrayL_h__ 1

#define __ARRAYL_H_VERSION__ 1.20

// Lower triangular Array class 

// Defining NDEBUG improves optimization.

#include <iostream.h>
#include <strstream.h>
#include "Array.h"

#ifdef NDEBUG
#define __check(i,n,dim,m)
#else
#define __check(i,n,dim,m) Check(i,n,dim,m)
#endif

typedef unsigned int uint;

namespace Array {
  
template<class T>
class array2L : public array1<T> {
  uint n;
 public:
  void Dimension(uint n0) {n=n0; size=n*(n+1)/2;}

  uint index(int i, int j) const {return i*(i+1)/2+j;}
  int N() const {return n;}
  
  array2L() {}
  array2L(uint n) {Allocate(n);}
  array2L(uint n, T *v0) {Dimension(n,v0);}
  array2L(const array2L<T>& A) {v=A.v; size=A.size; state=A.test(temporary);}
	
  array1<T> operator [] (uint i) const {
    __check(i,n,1,1);
    return array1<T>(i+1,v+index(i,0));
  }
	
  T& operator () (uint i, uint j) const {
    __check(i,n,2,1);
    __check(j,i+1,2,1);
    return v[index(i,j)];
  }
  T* operator () () const {return v;}
	
  array2L<T>& operator = (T a) {Load(a); return *this;}
};
  
template<class T>
ostream& operator << (ostream& s, const array2L<T>& A)
{
  T *p=A();
  for(unsigned int i=0; i < A.N(); i++) {
    for(unsigned int j=0; j <= i; j++) {
      s << *(p++) << " ";
    }
    s << _newl;
  }
  s << flush;
  return s;
}

}

#undef __check

#endif
