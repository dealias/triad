/* Arrayp.h:  A high-performance periodic multi-dimensional C++ array class
Copyright (C) 1999 John C. Bowman

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

#ifndef __Arrayp_h__
#define __Arrayp_h__ 1

#define __ARRAYP_H_VERSION__ 1.21

// Defining NDEBUG improves optimization.

#include "Array.h"

namespace Array {
  
inline void Mod(int& i, size_t n) {
  i %= (int) n;
  if(i < 0) i += (int) n;
}

template<class T>
class array1p : public array1<T> {
 public:
  array1p() {}
  array1p(size_t nx0) {this->Allocate(nx0);}
  array1p(size_t nx0, T *v0) {this->Dimension(nx0,v0);}
  T& operator [] (int ix) const {Mod(ix,this->size); return this->v[ix];}
  T& operator () (int ix) const {Mod(ix,this->size); return this->v[ix];}
  array1p<T>& operator = (T a) {this->Load(a); return *this;}
};

template<class T>
class array2p : public array2<T> {
 public:
  array2p() {}
  array2p(size_t nx0, size_t ny0) {this->Allocate(nx0,ny0);}
  array2p(size_t nx0, size_t ny0, T *v0) {
    this->Dimension(nx0,ny0,v0);}
  array1p<T> operator [] (int ix) const {
    Mod(ix,this->nx);
    return array1p<T>(this->ny,this->v+ix*this->ny);
  }
  T& operator () (int ix, int iy) const {
    Mod(ix,this->nx);
    Mod(iy,this->ny);
    return this->v[ix*this->ny+iy];
  }
  T& operator () (int i) const {
    Mod(i,this->size);
    return this->v[i];
  }
  array2p<T>& operator = (T a) {this->Load(a); return *this;}
};

template<class T>
class array3p : public array3<T> {
 public:	
  array3p() {}
  array3p(size_t nx0, size_t ny0, size_t nz0) {
    this->Allocate(nx0,ny0,nz0);
  }
  array3p(size_t nx0, size_t ny0, size_t nz0, T *v0) {
    this->Dimension(nx0,ny0,nz0,v0);
  }
  array2p<T> operator [] (int ix) const {
    Mod(ix,this->nx);
    return array2p<T>(this->ny,this->nz,this->v+ix*this->nyz);
  }
  T& operator () (int ix, int iy, int iz) const {
    Mod(ix,this->nx);
    Mod(iy,this->ny);
    Mod(iz,this->nz);
    return this->v[ix*this->nyz+iy*this->nz+iz];
  }
  T& operator () (int i) const {
    Mod(i,this->size);
    return this->v[i];
  }
};

template<class T>
class array4p : public array4<T> {
 public:	
  array4p() {}
  array4p(size_t nx0, size_t ny0, size_t nz0,
	  size_t nw0) {
    this->Allocate(nx0,ny0,nz0,nw0);
  }
  array4p(size_t nx0, size_t ny0, size_t nz0,
	  size_t nw0, T *v0) {
    this->Dimension(nx0,ny0,nz0,nw0,v0);
  }
  array3p<T> operator [] (int ix) const {
    Mod(ix,this->nx);
    return array3p<T>(this->ny,this->nz,this->nw,this->v+ix*this->nyzw);
  }
  T& operator () (int ix, int iy, int iz, int iw) const {
    Mod(ix,this->nx);
    Mod(iy,this->ny);
    Mod(iz,this->nz);
    Mod(iw,this->nw);
    return this->v[ix*this->nyzw+iy*this->nzw+iz*this->nw+iw];
  }
  T& operator () (int i) const {
    Mod(i,this->size);
    return this->v[i];
  }
};

}

#endif
