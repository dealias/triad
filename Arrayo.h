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

#ifndef __Arrayo_h__
#define __Arrayo_h__ 1

#define __ARRAYO_H_VERSION__ 1.08

// Defining NDEBUG improves optimization.

#include <iostream.h>
#include <strstream.h>
#include "Array.h"

#ifdef NDEBUG
#define check(i,n,o,dim,m)
#else
#define check(i,n,o,dim,m) Check(i-o,n,dim,m,o)
#endif

template<class T>
class Array1o : public Array1<T> {
protected:
	T *voff; // Offset pointer to memory block
	int ox;
public:
	void Offsets() {
		voff=v-ox;
	}
	void Dimension(unsigned int nx0, int ox0=0) {
		nx=nx0;
		ox=ox0;
	}
	void Dimension(unsigned int nx0, T *v0, int ox0=0) {
		Dimension(nx0,ox0);
		v=v0;
		Offsets();
		clear(allocated);
	}
	void Allocate(unsigned int nx0, int ox0=0) {
		Dimension(nx0,ox0);
		v=new T[Size()];
		Offsets();
		set(allocated);
	}
	
	Array1o() : ox(0) {}
	Array1o(unsigned int nx0, int ox0=0) {
		Allocate(nx0,ox0);
	}
	Array1o(unsigned int nx0, T *v0, int ox0=0) {
		Dimension(nx0,v0,ox0);
	}

	T& operator [] (int ix) const {check(ix,nx,ox,1,1); return voff[ix];}
	T& operator () (int i) const {check(i,nx,0,1,1); return v[i];}
	T* operator () () const {return voff;}
	operator T* () const {return voff;}
	
	Array1o<T> operator + (int i) const {return Array1o<T>(nx-i,v+i,ox+i);}
	void Set(T *a) {v=a; Offsets(); clear(allocated);}
	
	Array1o<T>& operator = (T a) {Load(a); return *this;}
	Array1o<T>& operator = (const T *a) {Load(a); return *this;}
	Array1o<T>& operator = (const Array1<T>& A) {
		Load(A());
		A.Purge();
		return *this;
	}
};

template<class T>
class Array2o : public Array2<T> {
protected:
	T *voff,*vtemp;
	int ox,oy;
public:
	void Offsets() {
		vtemp=v-ox*ny;
		voff=vtemp-oy;
	}
	void Dimension(unsigned int nx0, unsigned int ny0, int ox0=0, int oy0=0) {
		nx=nx0; ny=ny0;
		ox=ox0; oy=oy0;
	}
	void Dimension(unsigned int nx0, unsigned int ny0, T *v0, int ox0=0,
				   int oy0=0) {
		Dimension(nx0,ny0,ox0,oy0);
		v=v0;
		Offsets();
		clear(allocated);
	}
	void Allocate(unsigned int nx0, unsigned int ny0, int ox0=0, int oy0=0) {
		Dimension(nx0,ny0,ox0,oy0);
		v=new T[Size()];
		Offsets();
		set(allocated);
	}

	Array2o() : oy(0) {}
	Array2o(unsigned int nx0, unsigned int ny0, int ox0=0, int oy0=0) {
		Allocate(nx0,ny0,ox0,oy0);
	}
	Array2o(unsigned int nx0, unsigned int ny0, T *v0, int ox0=0, int oy0=0) {
		Dimension(nx0,ny0,v0,ox0,oy0);
	}
	
	Array1o<T> operator [] (int ix) const {
		check(ix,nx,ox,2,1);
		return Array1o<T>(ny,vtemp+ix*ny,oy);
	}
	T& operator () (int ix, int iy) const {
		check(ix,nx,ox,2,1);
		check(iy,ny,oy,2,2);
		return voff[ix*ny+iy];
	}
	T& operator () (int i) const {
		check(i,Size(),0,2,0);
		return v[i];
	}
	T* operator () () const {return voff;}
	void Set(T *a) {v=a; Offsets(); clear(allocated);}
	
	Array2o<T>& operator = (T a) {Load(a); return *this;}
	Array2o<T>& operator = (T *a) {Load(a); return *this;}
	Array2o<T>& operator = (const Array2<T>& A) {
		Load(A());
		A.Purge();
		return *this;
	}
};

template<class T>
class Array3o : public Array3<T> {
protected:
	T *voff,*vtemp;
	int ox,oy,oz;
public:
	void Offsets() {
		vtemp=v-ox*nyz;
		voff=vtemp-oy*nz-oz;
	}
	void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
				   int ox0=0, int oy0=0, int oz0=0) {
		nx=nx0; ny=ny0; nz=nz0;
		ox=ox0; oy=oy0; oz=oz0;
		nyz=ny*nz;
	}
	void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
				   T *v0, int ox0=0, int oy0=0, int oz0=0) {
		Dimension(nx0,ny0,nz0,ox0,oy0,oz0);
		v=v0;
		Offsets();
		clear(allocated);
	}
	void Allocate(unsigned int nx0, unsigned int ny0, unsigned int nz0,
				  int ox0=0, int oy0=0, int oz0=0) {
		Dimension(nx0,ny0,nz0,ox0,oy0,oz0);
		v=new T[Size()];
		Offsets();
		set(allocated);
	}
	
	Array3o() : oz(0) {}
	Array3o(unsigned int nx0, unsigned int ny0, unsigned int nz0,
		   int ox0=0, int oy0=0, int oz0=0) {
		Allocate(nx0,ny0,nz0,ox0,oy0,oz0);
	}
	Array3o(unsigned int nx0, unsigned int ny0, unsigned int nz0, T *v0,
		   int ox0=0, int oy0=0, int oz0=0) {
		Dimension(nx0,ny0,nz0,v0,ox0,oy0,oz0);
	}
	
	Array2o<T> operator [] (int ix) const {
		check(ix,nx,ox,3,1);
		return Array2o<T>(ny,nz,vtemp+ix*nyz,oy,oz);
	}
	T& operator () (int ix, int iy, int iz) const {
		check(ix,nx,ox,3,1);
		check(iy,ny,oy,3,2);
		check(iz,nz,oz,3,3);
		return voff[ix*nyz+iy*nz+iz];
	}
	T& operator () (int i) const {
		check(i,Size(),0,3,0);
		return v[i];
	}
	T* operator () () const {return voff;}
	void Set(T *a) {v=a; Offsets(); clear(allocated);}
	
	Array3o<T>& operator = (T a) {Load(a); return *this;}
	Array3o<T>& operator = (T *a) {Load(a); return *this;}
	Array3o<T>& operator = (const Array3<T>& A) {
		Load(A());
		A.Purge(); 
		return *this;
	}
};

template<class T>
class Array4o : public Array4<T> {
protected:
	T *voff,*vtemp;
	int ox,oy,oz,ow;
public:
	void Offsets() {
		vtemp=v-ox*nyzw;
		voff=vtemp-oy*nzw-oz*nw-ow;
	}
	void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
				   unsigned int nw0, 
				   int ox0=0, int oy0=0, int oz0=0, int ow0=0) {
		nx=nx0; ny=ny0; nz=nz0; nw=nw0;
		ox=ox0; oy=oy0; oz=oz0; ow=ow0;
		nzw=nz*nw; nyzw=ny*nzw;
	}
	void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
				   unsigned int w0, T *v0,
				   int ox0=0, int oy0=0, int oz0=0, int ow0=0) {
		Dimension(nx0,ny0,nz0,w0,ox0,oy0,oz0,ow0);
		v=v0;
		Offsets();
		clear(allocated);
	}
	void Allocate(unsigned int nx0, unsigned int ny0, unsigned int nz0,
				  unsigned int nw0,
				  int ox0=0, int oy0=0, int oz0=0, int ow0=0) {
		Dimension(nx0,ny0,nz0,nw0,ox0,oy0,oz0,ow0);
		v=new T[Size()];
		Offsets();
		set(allocated);
	}
	
	Array4o() : ow(0) {}
	Array4o(unsigned int nx0, unsigned int ny0, unsigned int nz0,
		   unsigned int nw0,
		   int ox0=0, int oy0=0, int oz0=0, int ow0=0) {
		Allocate(nx0,ny0,nz0,nw0,ox0,oy0,oz0,ow0);
	}
	Array4o(unsigned int nx0, unsigned int ny0, unsigned int nz0,
		   unsigned int nw0, T *v0,
		   int ox0=0, int oy0=0, int oz0=0, int ow0=0) {
		Dimension(nx0,ny0,nz0,nw0,v0,ox0,oy0,oz0,ow0);
	}

	Array3o<T> operator [] (int ix) const {
		check(ix,nx,ox,3,1);
		return Array3o<T>(ny,nz,nw,vtemp+ix*nyzw,oy,oz,ow);
	}
	T& operator () (int ix, int iy, int iz, int iw) const {
		check(ix,nx,ox,4,1);
		check(iy,ny,oy,4,2);
		check(iz,nz,oz,4,3);
		check(iw,nw,ow,4,4);
		return voff[ix*nyzw+iy*nzw+iz*nw+iw];
	}
	T& operator () (int i) const {
		check(i,Size(),0,4,0);
		return v[i];
	}
	T* operator () () const {return voff;}
	void Set(T *a) {v=a; Offsets(); clear(allocated);}
	
	Array4o<T>& operator = (T a) {Load(a); return *this;}
	Array4o<T>& operator = (T *a) {Load(a); return *this;}
	Array4o<T>& operator = (const Array4<T>& A) {
		Load(A());
		A.Purge();
		return *this;
	}
};

#ifdef NDEBUG
#define Array1o(T) T*
#else
#define Array1o(T) Array1o<T>
#endif

#undef check

#endif

