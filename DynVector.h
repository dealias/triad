#ifndef __DynVector_h__
#define __DynVector_h__ 1

#include <stdlib.h>

#if __AIX
#define CONST const
#else
#define mutable
#define CONST
#endif

template<class T>
class DynVector
{
	T *v;
	unsigned int size;
	unsigned int alloc;
	mutable int state;
public:
    enum alloc_state {unallocated=0, allocated=1, temporary=2};
	void Allocate(unsigned int s) {v=new T[alloc=s]; size=0; set(allocated);}
	void Deallocate() {
		if(test(allocated)) delete [] v;
		size=0; alloc=0; clear(allocated);
	}
	
	int test(int flag) const {return state & flag;}
	void clear(int flag) CONST {state &= ~flag;}
	void set(int flag) CONST {state |= flag;}
	
	DynVector() : v(NULL), size(0), alloc(0), state(unallocated) {}
	DynVector(const DynVector& A) : v(A.v), size(A.size),
		alloc(A.alloc), state(A.test(temporary)) {}
	DynVector(unsigned int s) {Allocate(s);}
	~DynVector() {Deallocate();}

	void Freeze() {state=unallocated;}
	void Hold() {if(test(allocated)) {state=temporary;}}
	void Purge() CONST {if(test(temporary)) {Deallocate(); state=unallocated;}}
#ifdef mutable
	void Purge() const {((DynVector<T> *) this)->Purge();}
#endif
	
	unsigned int Alloc() const {return alloc;}
	unsigned int Size() const {return size;}
	
	void Resize(unsigned int i) {
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

	void Load(T a) const {for(unsigned int i=0; i < size; i++) v[i]=a;}
	void Load(T *a) const {memcpy(v,a,size*sizeof(T));}
	void Store(T *a) const {memcpy(a,v,size*sizeof(T));}
	void Set(T *a) {v=a; clear(allocated);}
	
	DynVector<T>& operator = (T a) {Load(a); return *this;}
	DynVector<T>& operator = (const T *a) {Load(a); return *this;}
	DynVector<T>& operator = (const DynVector<T>& A) {
		Load(A()); 
		A.Purge();
		return *this;
	}
	
};

#undef CONST
#endif
