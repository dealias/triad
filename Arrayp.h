/* Arrayp.h:  A high-performance multi-dimensional C++ array class
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

#ifndef __Arrayp_h__
#define __Arrayp_h__ 1

#define __ARRAYP_H_VERSION__ 1.07

// Defining NDEBUG improves optimization.

#include <iostream.h>
#include <strstream.h>
#include "Array.h"

class ArrayMod {
public:	
	void Mod(int& i, int n) const {
		if(i >= n) i -= n;
		else if(i < 0) i += n;
	
		if(i >= n) i %= n;
		else if(i < 0) {i=n-(-i % n); if(i == n) i=0;}
	}
};

template<class T>
class Array1p : public Array1<T>, public ArrayMod {
public:
	Array1p() {}
	Array1p(unsigned int nx0) {Allocate(nx0);}
	Array1p(unsigned int nx0, T *v0) {Dimension(nx0,v0);}
	Array1p(const Array1p<T>& A) {v=A.v; nx=A.nx; state=A.test(temporary);}
	T& operator [] (int ix) const {Mod(ix,nx); return v[ix];}
	T& operator () (int ix) const {Mod(ix,nx); return v[ix];}
};

template<class T>
class Array2p : public Array2<T>, public ArrayMod {
public:
	Array2p() {}
	Array2p(unsigned int nx0, unsigned int ny0) {Allocate(nx0,ny0);}
	Array2p(unsigned int nx0, unsigned int ny0, T *v0) {Dimension(nx0,ny0,v0);}
	Array1p<T> operator [] (int ix) const {
		Mod(ix,nx);
		return Array1p<T>(ny,v+ix*ny);
	}
	T& operator () (int ix, int iy) const {
		Mod(ix,nx);
		Mod(iy,ny);
		return v[ix*ny+iy];
	}
	T& operator () (int i) const {
		Mod(i,Size());
		return v[i];
	}
};

template<class T>
class Array3p : public Array3<T>, public ArrayMod {
public:	
	Array3p() {}
	Array3p(unsigned int nx0, unsigned int ny0, unsigned int nz0) {Allocate(nx0,ny0,nz0);}
	Array3p(unsigned int nx0, unsigned int ny0, unsigned int nz0, T *v0) {Dimension(nx0,ny0,nz0,v0);}
	Array2p<T> operator [] (int ix) const {
		Mod(ix,nx);
		return Array2p<T>(ny,nz,v+ix*nyz);
	}
	T& operator () (int ix, int iy, int iz) const {
		Mod(ix,nx);
		Mod(iy,ny);
		Mod(iz,nz);
		return v[ix*nyz+iy*nz+iz];
	}
	T& operator () (int i) const {
		Mod(i,Size());
		return v[i];
	}
};

template<class T>
class Array4p : public Array4<T>, public ArrayMod {
public:	
	Array4p() {}
	Array4p(unsigned int nx0, unsigned int ny0, unsigned int nz0, unsigned int nw0) {Allocate(nx0,ny0,nz0,nw0);}
	Array4p(unsigned int nx0, unsigned int ny0, unsigned int nz0, unsigned int nw0, T *v0) {
		Dimension(nx0,ny0,nz0,nw0,v0);
	}
	Array3p<T> operator [] (int ix) const {
		Mod(ix,nx);
		return Array3p<T>(ny,nz,nw,v+ix*nyzw);
	}
	T& operator () (int ix, int iy, int iz, int iw) const {
		Mod(ix,nx);
		Mod(iy,ny);
		Mod(iz,nz);
		Mod(iw,nw);
		return v[ix*nyzw+iy*nzw+iz*nw+iw];
	}
	T& operator () (int i) const {
		Mod(i,Size());
		return v[i];
	}
};

#ifdef NDEBUG
#define Array1p(T) T*
#else
#define Array1p(T) Array1p<T>
#endif

#endif

