#ifndef __Triad_h__
#define __Triad_h__ 1

#if(MCREAL)
typedef Real Mc;
typedef float McWeight;
#else
typedef Complex Mc;
typedef Complex McWeight;
#endif

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

extern int Ntriad;
extern DynVector<Triad> triad;
extern TriadLimits *triadLimits;

extern Var *pqbuffer;
extern Var **pqIndex;
extern int *qStart;

#endif
