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
class array1H : public array1<T> {
private:
  unsigned nx0m1;
public:
  array1H() {}
  array1H(unsigned int mx) {this->Allocate(mx);}
  array1H(unsigned int mx, T *v0) {this->Dimension(mx,v0);}
  
  void Check(int i, int n, unsigned int dim, unsigned int m) const {
    if(i <= -n || i >= n) {
      std::ostringstream buf;
      buf << "ArrayH" << dim << " index ";
      if(m) buf << m << " ";
      buf << "is out of bounds (" << i;
      if(n == 0) buf << " index given to empty array";
      else {
	if(i < 0) buf << " < " << 1-n;
	else buf << " > " << n-1;
      }
      buf << ")";
      ArrayExit(buf.str().c_str());
    }
  }
	
  T get(int ix) const {
    return ix >= 0 ? this->v[ix] : (ix < 0 ? conj(this->v[-ix]) : this->v[0].re);
  }
    
  T operator [] (int ix) const {
    __check(ix,this->size,1,1);
    return get(ix);
  }
  array1H<T>& operator = (T a) {Load(a); return *this;}
};
  
template<class T>
class array2H : public array2<T> {
private:
  unsigned xorigin;
    
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
  

public:
  array2H() {}
  array2H(unsigned int nx0, unsigned int ny0) {
    xorigin=(nx0-1)/2;
    this->Allocate(nx0,ny0);
  }
  array2H(unsigned int nx0, unsigned int ny0, T *v0) {
    xorigin=(nx0-1)/2;
    this->Dimension(nx0,ny0,v0);
  }
  array1H<T> operator [] (int ix) const {
    // FIXME: return cc for certain cases?
    return array1H<T>(this->ny,this->v+ix*this->ny);
  }
    
  T& operator () (int ix, int iy) const {
    // FIXME: return cc for certain cases?
    return this->v[ix*this->ny+iy];
  }
    
  T get(int ix, int iy) {
    __checkH(ix,iy);
    if(iy > 0)
      return this->v[ix*this->ny+iy];
    if(iy == 0) {
      if(ix < xorigin) return conj(this->v[(2*xorigin-ix)*this->ny]);
      if(ix > xorigin) return this->v[ix*this->ny];
      return this->v[(2*xorigin-ix)*this->ny].re;
    }
    return conj(this->v[(2*xorigin-ix)*this->ny-iy]);
  }

  void set(int ix, int iy, T a) {
    __checkH(ix,iy);
    if(iy > 0) {this->v[ix*this->ny+iy]=a; return;}
    if(iy == 0) {
      if(ix == xorigin) {
        this->v[(2*xorigin-ix)*this->ny] = a.re; // FIXME: should be .re
      } else {
        this->v[(2*xorigin-ix)*this->ny]=conj(a);
        this->v[ix*this->ny] = a;
      }
      return;
    }
    this->v[(2*xorigin-ix)*this->ny-iy]=conj(a);
  }

  array2H<T>& operator = (T a) {Load(a); return *this;}
};
  

}

#endif
