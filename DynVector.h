#ifndef __DynVector_h__
#define __DynVector_h__ 1

#include "utils.h"

template<class T>
class DynVector
{
	T *v;
	int sz;

public:
	void Allocate(int s) {v=new T[sz=s];}
	void Deallocate() {delete [] v; sz=0;}
	
	DynVector() : v(NULL), sz(0) {}
	DynVector(int s) {Allocate(s);}
	~DynVector() {Deallocate();}

	int Size() const {return sz;}
	T& operator [] (int i) {
		if (i >= sz) v=new(v,sz=max(i+1,2*sz)) (T);
		return v[i];
	}
	T *operator + (int i) {return v+i;}
	T *operator () () const {return v;}
	operator T* () const {return v;}
	
	void Resize(int i) {v=new(v,sz=i) (T);}
	
	void Expand(int i) {if (i > sz) Resize(i);}

	DynVector<T> operator = (const T *A) {
		memcpy(v,A,sz*sizeof(T));
		return *this;
	}
};

#endif
