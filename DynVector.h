#ifndef __DynVector_h__
#define __DynVector_h__ 1

/* Hello */

#include <stdlib.h>

template<class T>
class DynVector
{
	T *v;
	unsigned int size;
	unsigned int alloc;
	int allocated;
public:
	void Allocate(unsigned int s) {v=new T[alloc=s]; size=0; allocated=1;}
	void Deallocate() {delete [] v; size=0; alloc=0; allocated=0;}
	
	DynVector() : v(NULL), size(0), alloc(0), allocated(0) {}
	DynVector(const DynVector& A) : v(A.v), size(A.size),
		alloc(A.alloc), allocated(0) {}
	DynVector(unsigned int s) {Allocate(s);}
	~DynVector() {if(allocated) Deallocate();}

	unsigned int Alloc() const {return alloc;}
	unsigned int Size() const {return size;}
	
	void Resize(unsigned int i) {
 		if (i == 0 && alloc) delete [] v;
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
	
	void Push(const T& value) {
		if (size == alloc) Resize(alloc ? 2*alloc : 1);
		v[size++]=value;
	}
	
	void Pop() {
		if(size) size--;
	}

	void Expand(unsigned int i) {if (i > alloc) Resize(i);}

	DynVector<T> operator = (const T *A) {
		memcpy(v,A,size*sizeof(T));
		return *this;
	}
};

#endif
