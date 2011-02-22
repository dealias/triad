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

namespace Array {
  
  bool noconj(int ix, int xorigin) {
    if(ix >= xorigin) return true;
    return false;
  }
  
  template<class T>
    class array1H : public array1<T> {
  private:
    unsigned xorigin;
   public:
    array1H() {}
    array1H(unsigned int nx0) {
      if(nx0 % 2 == 0) exit(1); // FIXME: improve error message
      xorigin=(nx0-1)/2;
      this->Allocate(nx0);
    }
    array1H(unsigned int nx0, T *v0) {
      xorigin=(nx0-1)/2;
      this->Dimension(nx0,v0);
    }
    T& operator [] (int ix) const {
      // FIXME: return conjugate if ix < xorigin?
      return this->v[ix];
    }
    
    //T& operator () (int ix) const {return this->v[ix];}
    
    T get(int ix) {
      if(ix == xorigin)
	return 0.5*(this->v[ix]+conj(this->v[ix]));
      if(noconj(ix,xorigin))
	return this->v[ix];
      return conj(this->v[2*xorigin-ix]);
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
