#ifndef __Pair_h__
#define __Pair_h__ 1

#include "types.h"
#include "DynVector.h"

class Pair {
public:
	Var *p, *q;
	void Store(Var *p0, Var *q0) {p=p0; q=q0;}
};

class Triad {
public:
	Var *pq;
	Mc Mkpq;
	void Store(Var *pqindex, Mc value) {pq=pqindex; Mkpq=value;}
};

extern int Npsi;
extern int Npair;
extern int Ntriad;

extern int reality;

extern Var *psibuffer,*pqbuffer;
extern DynVector<Pair> pair;
extern DynVector<Triad> triad;
extern Pair *pairBase;
extern Triad *triadBase;
extern Triad **triadStop;

#endif
