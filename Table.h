#ifndef __Table_h__
#define __Table_h__ 1

#include <iostream>
#include <cstdio>
#include <string>

using std::string;

#include "DynVector.h"

template<class B>
class EntryBase {
 protected:	
  const char *key;
 public:
  const char *Key() {return key;}
  virtual B* New()=0;
};

template<class B>
class Table {
  const char *name;
  Compare_t *Compare;
  KeyCompare_t *KeyCompare;
  unsigned int n;
  DynVector<EntryBase<B> *>list;
 public:
  static int DefaultCompare(const void *a, const void *b) {
    return strcasecmp((*(EntryBase<B> **)a)->Key(),
		      (*(EntryBase<B> **)b)->Key());
  }
  static int DefaultKeyCompare(const void *key, const void *p,
			       const size_t n) {
    return strcasecmpn((char *) key, (*(EntryBase<B> **)p)->Key(), n);
  }
	
  Table(const char *name, Compare_t Compare, KeyCompare_t KeyCompare) :
    name(name), Compare(Compare), KeyCompare(KeyCompare), n(0) {}
  Table(const char *name) : name(name), Compare(DefaultCompare),
			    KeyCompare(DefaultKeyCompare), n(0) {}
  void Add(EntryBase<B> *ptr) {list[n++]=ptr;}
  unsigned int Size() {return n;}
  EntryBase<B> **Base() {return list;}
  EntryBase<B> *Entry(unsigned int i) {return list[i];}
  void List(ostream& os) {
    for(unsigned int i=0; i < Size(); i++) os << Entry(i)->Key() << newl;
    os << flush;
  }
  B *Locate (const char *& key) {
    int match_type;
    EntryBase<B> *e;
    B *p;
    qsort(Base(),Size(),sizeof(*Base()),Compare);
    
    for(;;) {
	e=*(EntryBase<B> **) bsearch2(key,Base(),Size(),sizeof(*Base()),
				      KeyCompare,&match_type);
	if(check_match(match_type,name,key,2)) break;
	
	cout << "Recognized " << name << " values:" << endl << endl;
	List(cout);
	cout << newl << name << "=";
      
	string s;
	getline(cin,s);
	if(s.c_str()) key=strdup(s.c_str());
    }
    
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
  Entry(const char *key0, Table<B> *t) {this->key=key0; t->Add(this);}
  B *New() {return new T;}
};

template<class T, class B, class P>
class entry : public EntryBase<B> {
  P *p;
public:
  entry(const char *key0, Table<B> *t, P *p) : p(p) {
  this->key=key0; t->Add(this);
}
  B *New() {return new T(p);}
};

#endif
