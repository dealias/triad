#ifndef __DynVector_h__
#define __DynVector_h__ 1

#include <stdlib.h>

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
	
	void Resize(int i) {
		if (i == 0) free(v);
		else if(i > sz) {
			T *v0=v;
			v=new T[i];
			for(int j=0; j < sz; j++) v[j]=v0[j];
			delete [] v0;
		}
		sz=i;
	}
	
	T& operator [] (int i) {
		if (i >= sz) Resize(max(i+1,2*sz));
		return v[i];
	}
	
	T *operator + (int i) {return v+i;}
	T *operator () () const {return v;}
	operator T* () const {return v;}
	
	void Expand(int i) {if (i > sz) Resize(i);}

	DynVector<T> operator = (const T *A) {
		memcpy(v,A,sz*sizeof(T));
		return *this;
	}
};

#endif
