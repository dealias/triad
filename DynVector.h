/* Dynvector.h:  A simple dynamic vector class
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

#ifndef __DynVector_h__
#define __DynVector_h__ 1

#define __DYNVECTOR_H_VERSION__ 1.00

#include <cstdlib>
#include <iostream>

using std::cout;
using std::endl;
using std::ostream;

#ifndef __ExternalDynVectorExit
inline void DynVectorExit(char *x)
{
  cout << endl << "ERROR: " << x << "." << endl;
  exit(1);
} 
#endif

template<class T>
class DynVector
{
protected:	
  T *v;
  mutable unsigned int size;
  mutable unsigned int alloc;
  mutable int state;
public:
  enum alloc_state {unallocated=0, allocated=1, temporary=2};
  void Allocate(unsigned int s) {
    v=new T[alloc=s]; size=0; set(allocated);
  }
  void Deallocate() const {
    if(alloc && test(allocated)) delete [] v;
    size=0; alloc=0; clear(allocated);
  }
	
  int test(int flag) const {return state & flag;}
  void clear(int flag) const {state &= ~flag;}
  void set(int flag) const {state |= flag;}
	
  DynVector() : v(NULL), size(0), alloc(0), state(unallocated) {}
  void SetDynVector(const DynVector<T>& A) {v=A.v; size=A.size;
                                     alloc=A.alloc; state(A.test(temporary));}
  DynVector(const DynVector<T>& A) {SetDynVector(A);}
  DynVector(unsigned int s) {Allocate(s);}
  ~DynVector() {Deallocate();}

  void Freeze() {state=unallocated;}
  void Hold() {if(test(allocated)) {state=temporary;}}
  void Purge() const {
    if(test(temporary)) {Deallocate(); state=unallocated;}
  }
	
  unsigned int Alloc() const {return alloc;}
  unsigned int Size() const {return size;}
  
  void Realloc(unsigned int i) {
    if (i == 0 && alloc && test(allocated)) delete [] v;
    else if(i > alloc) {
      T *v0=v;
      v=new T[i];
      if (size) {
	for(unsigned int j=0; j < size; j++) v[j]=v0[j];
	if(test(allocated)) delete [] v0;
      }
    }
    alloc=i;
    if(size > alloc) size=alloc;
  }
  
  void Size(unsigned int i) {
    if(i > alloc) Realloc(i);
    size=i;
  }
	
  int max(int a, int b)
  {
    return (a > b) ? a : b;
  }

  T& operator [] (unsigned int i) {
    if (i >= alloc) Realloc(max(i+1,2*alloc));
    if (i >= size) size=i+1;
    return v[i];
  }
	
  T& operator [] (int i) {
    return (*this)[(unsigned int) i];
  }
	
  T *operator + (unsigned int i) {return v+i;}
  T *operator + (int i) {return v+i;}
  T *operator () () const {return v;}
  operator T* () const {return v;}
	
  void Push(const T& value) {
    if (size == alloc) Realloc(alloc ? 2*alloc : 1);
    v[size++]=value;
  }
	
  int Pop(T& value) {
    if(size) {
      size--;
      value=v[size];
      return 0;
    } else return -1;
  }
  
  int Pop() {
    if(size) {
      size--;
      return 0;
    } else return -1;
  }
  
  // Pop v[i], close up, and return v[i].
  int Pop(T& value, unsigned int i) {
    if(size) {
      value=v[i];
     for (unsigned int j=i+1; j < size; j++) v[j-1]=v[j];
      size--;
      return 0;
    } else return -1;
  }
  
  void Expand(unsigned int i) {if (i > alloc) Realloc(i);}

  void Load(T a) const {for(unsigned int i=0; i < size; i++) v[i]=a;}
  void Load(const T *a) const {memcpy(v,a,size*sizeof(T));}
  void Store(T *a) const {memcpy(a,v,size*sizeof(T));}
  void Set(T *a, unsigned int n) {v=a; alloc=n; clear(allocated);}
	
  DynVector<T>& operator = (T a) {Load(a); return *this;}
  DynVector<T>& operator = (const T *a) {Load(a); return *this;}
  DynVector<T>& operator = (const DynVector<T>& A) {
    Load(A()); 
    A.Purge();
    return *this;
  }

  ostream& List(ostream& os, int nperline=1) {
    for(unsigned int i=0; i < size-1;) {
      os << v[i];
      if((++i % nperline) == 0) os << endl;
    }
    os << v[size-1] << endl;
    return os;
  }
};

template<class T>
class StackVector : public DynVector<T> {
public:
  StackVector() {}
  StackVector(const DynVector<T>& A) {DynVector<T>::SetDynVector(A);}
  StackVector(unsigned int s) {Allocate(s);}
  
  T& operator [] (unsigned int i) {
    if (i >= size) DynVectorExit("Attempt to access past end of StackVector");
    return v[i];
  }
	
  T& operator [] (int i) {
    return (*this)[(unsigned int) i];
  }
	
  StackVector<T>& operator = (T a) {Load(a); return *this;}
  StackVector<T>& operator = (const T *a) {Load(a); return *this;}
  StackVector<T>& operator = (const StackVector<T>& A) {
    Load(A()); 
    A.Purge();
    return *this;
  }
};

#endif
