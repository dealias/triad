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

struct TriadLimits {
	Triad *start;
	Triad *stop;
};

extern int Npsi;   /* number of explictly evolved modes */
extern int NpsiR;  /* number of reflected modes */
extern int Ntotal; /* total number of (evolved+reflected) modes */
extern int Ntriad;

extern int reality;

extern Var *psibuffer,*psibufferR,*psibufferStop,*pqbuffer;
extern Var **pqIndex;
extern int *qStart;
extern DynVector<Triad> triad;
extern TriadLimits *triadLimits;

#endif
