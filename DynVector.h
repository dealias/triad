#ifndef __DynVector_h__
#define __DynVector_h__ 1

#include <stdlib.h>
#include <iostream.h>

template<class T>
class DynVector
{
protected:	
  T *v;
  mutable unsigned int size;
  mutable unsigned int alloc;
  mutable int state;
public:
  enum alloc_state {normal=0, temporary=1};
  void Allocate(unsigned int s) {v=new T[alloc=s]; size=0;}
  void Deallocate() const {
    if(alloc) delete [] v;
    size=0; alloc=0;
  }
	
  int test(int flag) const {return state & flag;}
  void clear(int flag) const {state &= ~flag;}
  void set(int flag) const {state |= flag;}
	
  DynVector() : v(NULL), size(0), alloc(0), state(normal) {}
  DynVector(const DynVector<T>& A) : v(A.v), size(A.size),
				     alloc(A.alloc), state(A.test(temporary)) {}
  DynVector(unsigned int s) {Allocate(s);}
  ~DynVector() {Deallocate();}

  void Freeze() {state=normal;}
  void Hold() {if(alloc) {state=temporary;}}
  void Purge() const {if(test(temporary)) {Deallocate(); state=normal;}}
	
  unsigned int Alloc() const {return alloc;}
  unsigned int Size() const {return size;}
	
  void Resize(unsigned int i) {
    if (i == 0 && alloc) Deallocate();
    else if(i > alloc) {
      T *v0=v;
      v=new T[i];
      if (size) {
	for(unsigned int j=0; j < size; j++) v[j]=v0[j];
	if(alloc) delete [] v0;
      }
    }
    alloc=i;
    if(size > alloc) size=alloc;
  }

  void SetTop(unsigned int i){
    if (alloc < i) Resize(i);
    size=i;
  }

  T& operator [] (unsigned int i) {
    if (i >= alloc) Resize(max(i+1,2*alloc));
    if (i >= size) size=i+1;
    return v[i];
  }
	
  T *operator + (unsigned int i) {return v+i;}
  T *operator + (int i) {return v+i;}
  T *operator () () const {return v;}
  operator T* () const {return v;}
	
  void Push(const T& value) {
    if (size == alloc) Resize(alloc ? 2*alloc : 1);
    v[size++]=value;
  }
	
  void Pop() {
    if(size) size--;
  }
  
  void Pop(unsigned int i) {
    if(size) {
      for (unsigned int j=i; j < size; j++) v[j-1]=v[j];
      size--;
    }
  }

  void Expand(unsigned int i) {if (i > alloc) Resize(i);}

  void Load(T a) const {for(unsigned int i=0; i < size; i++) v[i]=a;}
  void Load(const T *a) const {memcpy(v,a,size*sizeof(T));}
  void Store(T *a) const {memcpy(a,v,size*sizeof(T));}
  void Set(T *a) {v=a; alloc=0;}
	
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
    os << v[size-1] << newl;
    return os;
  }
};

#endif
