#ifndef __DynVector_h__
#define __DynVector_h__ 1

#include <stdlib.h>

template<class T>
class DynVector
{
	T *v;
	int sz;
	int top;
public:
	void Allocate(int s) {v=new T[sz=s]; top=-1;}
	void Deallocate() {delete [] v; sz=0;}
	
	DynVector() : v(NULL), sz(0), top(-1) {}
	DynVector(int s) {Allocate(s);}
	~DynVector() {Deallocate();}

	int Capacity() const {return sz;}
	int Size() const {return top;}
	
	void Resize(int i) {
 		if (i == 0 && sz > 0) delete [] v;
		else if(i > sz) {
			T *v0=v;
			v=new T[i];
			if (sz) {
				for(int j=0; j < sz; j++) v[j]=v0[j];
				delete [] v0;
			}
		}
		sz=i;
		if(top >= sz) top=sz-1;
	}
	
	T& operator [] (int i) {
		if (i >= sz) Resize(max(i+1,2*sz));
		if (i > top) top=i;
		return v[i];
	}
	
	T *operator + (int i) {return v+i;}
	T *operator () () const {return v;}
	operator T* () const {return v;}
	
	void PushBack(T value) {
		top++;
		if (top >= sz) Resize(max(top+1,2*sz));
		v[top]=value;
	}
	
	void Expand(int i) {if (i > sz) Resize(i);}

	DynVector<T> operator = (const T *A) {
		memcpy(v,A,top*sizeof(T));
		return *this;
	}
};

#endif
