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
	void DeAllocate() {if(v) {delete [] v; v=NULL; sz=0;}}
	
	DynVector() {v=NULL; sz=0;}
	DynVector(int s) {Allocate(s);}
	~DynVector() {DeAllocate();}

	int Size() const {return sz;}
	T *Base() const {return v;}
	T& operator [] (int i) {
		if (i >= sz) v=new(v,sz=max(i+1,2*sz)) (T);
		return v[i];
	}
	T *operator + (int i) {return v+i;}
	
	void Resize(int i) {v=new(v,sz=i) (T);}

	DynVector<T> operator = (const T *A) {
		memcpy(v,A,sz*sizeof(T));
		return *this;
	}
};

#endif
