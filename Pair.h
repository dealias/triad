#ifndef __Pair_h__
#define __Pair_h__ 1

#include "types.h"
#include "DynVector.h"

class Triad {
public:
	Var *pq;
	Mc Mkpq;
	void Store(Var *pqindex, Mc value) {pq=pqindex; Mkpq=value;}
};

extern int Npsi;
extern int Ntriad;

extern int reality;

extern Var *psibuffer,*psibufferStop,*pqbuffer;
extern DynVector<Triad> triad;
extern Triad *triadBase;
extern Triad **triadStop;

#endif
