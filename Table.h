#ifndef __Table_h__
#define __Table_h__ 1

#include <iostream.h>
#include <stdio.h>

#include "DynVector.h"

template<class B>
class EntryBase {
protected:	
	char *key;
public:
	char *Key() {return key;}
	virtual B* New()=0;
};

template<class B>
class Table {
	char *name;
	Compare_t *Compare;
	KeyCompare_t *KeyCompare;
	DynVector<EntryBase<B> *>list;
	unsigned int n;
public:
	static int DefaultCompare(const void *a, const void *b) {
		return strcasecmp((*(EntryBase<B> **)a)->Key(),
						  (*(EntryBase<B> **)b)->Key());
	}
	static int DefaultKeyCompare(const void *key, const void *p,
								 const size_t n) {
		return strcasecmpn((char *) key, (*(EntryBase<B> **)p)->Key(), n);
	}
	
	Table(char *name0, Compare_t Compare0, KeyCompare_t KeyCompare0) {
		name=name0; n=0; Compare=Compare0; KeyCompare=KeyCompare0;
	}
	Table(char *name0) {
		name=name0; n=0; Compare=DefaultCompare; KeyCompare=DefaultKeyCompare;
	}
	void Add(EntryBase<B> *ptr) {list[n++]=ptr;}
	unsigned int Size() {return n;}
	EntryBase<B> **Base() {return list;}
	EntryBase<B> *Entry(int i) {return list[i];}
	void List(ostream& os) {
		for(unsigned int i=0; i < Size(); i++) os << Entry(i)->Key() << newl;
		os << flush;
	}
	B *Locate (const char *& key) {
		int match_type;
		EntryBase<B> *e;
		B *p;
		qsort(Base(),Size(),sizeof(*Base()),Compare);
		e=*(EntryBase<B> **) bsearch2(key,Base(),Size(),sizeof(*Base()),
									  KeyCompare,&match_type);
		check_match(match_type,name,key);
		p=e->New();
		if(*p->Name()) {
			key=e->Key();
			cout << newl << upcase(name) << ": " << p->Name() << endl;
		}
		return p;
	}
	B *New() {return NULL;}
};

template<class T, class B>
class Entry : public EntryBase<B> {
public:
	Entry(char *key0,Table<B> *t) {key=key0; t->Add(this);}
	B *New() {return new T;}
};

#endif
