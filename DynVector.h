#ifndef __DynVector_h__
#define __DynVector_h__ 1

#include "utils.h"

template<class T>
class DynVector
{
	T *v;
	int sz;

public:
	DynVector() {v=NULL; sz=0;}
	DynVector(int s) {v=new T[sz=s];}
	~DynVector() {delete [] v; sz=0;}

	int Size() const {return sz;}
	T *Base() const {return v;}
	T& operator [] (int i) {
		if (i >= sz) v=new(v,sz=max(i+1,2*sz)) (T);
		return v[i];
	}
	
	void Resize(int i) {v=new(v,sz=i) (T);}

	DynVector<T> operator = (const T *A) {
		memcpy(v,A,sz*sizeof(T));
		return *this;
	}
};

template<class T>
ostream& operator << (ostream& s, DynVector<T>& A)
{
	int i, sz=A.Size();
	for(i=0; i < sz; i++) s << A[i] << " ";
	return s;
}

#endif
