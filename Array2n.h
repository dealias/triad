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

#ifndef __Array2n_h__
#define __Array2n_h__ 1

#define __ARRAY2N_H_VERSION__ 1.11

// Defining NDEBUG improves optimization.

#include <iostream.h>
#include <strstream.h>
#include "Array.h"
#include "Arrayp.h"

#ifdef NDEBUG
#define __check(i,n,dim,m)
#else
#define __check(i,n,dim,m) Check(i,n,dim,m)
#endif

typedef unsigned int uint;
const uint one=1;

template<class T>
class Array2n : public array1<T>, public ArrayMod {
 public:
  void Dimension(uint n) {size=(one << n)-1;}

  Array2n() {}
  Array2n(uint n) {Allocate(n);}
  Array2n(uint n, T *v0) {Dimension(n,v0);}
  Array2n(const Array2n<T>& A) {v=A.v; size=A.size; state=A.test(temporary);}
	
  Array1<T> operator [] (uint j) const {
    uint m=one << j;
    __check(2*(m-1),size,1,1);
    return Array1<T>(m,v+m-1);
  }
	
  T& operator () (uint j, int k) const {
    uint m=one << j;
    __check(2*(m-1),size,1,1);
//    __check(k,m,1,1);
    Mod(k,m);
    return v[m+k-1];
  }
	
  Array2n<T>& operator = (T a) {Load(a); return *this;}
};

#undef __check

#endif
