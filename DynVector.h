#ifndef __DynVector_h__
#define __DynVector_h__ 1

#include <stdlib.h>

template<class T>
class DynVector
{
	T *v;
	unsigned int alloc;
	unsigned int size;
public:
	void Allocate(unsigned int s) {v=new T[alloc=s]; size=0;}
	void Deallocate() {delete [] v; alloc=0;}
	
	DynVector() : v(NULL), alloc(0), size(0) {}
	DynVector(unsigned int s) {Allocate(s);}
	~DynVector() {Deallocate();}

	unsigned int Capacity() const {return alloc;}
	unsigned int Size() const {return size;}
	
	void Resize(unsigned int i) {
 		if (i == 0 && alloc > 0) delete [] v;
		else if(i > alloc) {
			T *v0=v;
			v=new T[i];
			if (alloc) {
				for(unsigned int j=0; j < alloc; j++) v[j]=v0[j];
				delete [] v0;
			}
		}
		alloc=i;
		if(size > alloc) size=alloc;
	}
	
	T& operator [] (unsigned int i) {
		if (i >= alloc) Resize(max(i+1,2*alloc));
		if (i >= size) size=i+1;
		return v[i];
	}
	
	T *operator + (unsigned int i) {return v+i;}
	T *operator () () const {return v;}
	operator T* () const {return v;}
	
	void PushBack(T value) {
		if (size == alloc) Resize(2*alloc+1);
		v[size++]=value;
	}
	
	void Expand(unsigned int i) {if (i > alloc) Resize(i);}

	DynVector<T> operator = (const T *A) {
		memcpy(v,A,size*sizeof(T));
		return *this;
	}
};

#endif
