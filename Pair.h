#ifndef __Pair_h__
#define __Pair_h__ 1

#include "types.h"
#include "DynVector.h"

class Pair {
public:
	Var *p, *q;
	Var psipq;
	Var *Store(Var *p0, Var *q0) {p=p0; q=q0; return Index();}
	Var *Index() {return &psipq;}
};

#if(0)
class Triad {
public:
	Var *pq;
	Mc Mkpq;
	void Store(Var *pqindex, Mc value) {pq=pqindex; Mkpq=value;}
};
#endif

extern int Npsi;
extern int Npair;
extern int Ntriad;

extern int reality;

extern Var *psibuffer;
extern Pair *pair;
extern DynVector<Mc> Mkpq;
extern DynVector<Var *> ptrpq;
extern Mc *MkpqBase,**MkpqStop;
extern Var **ptrpqBase;

#endif
