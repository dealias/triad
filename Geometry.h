#ifndef __Geometry_h__
#define __Geometry_h__ 1

#include "Mode.h"

#if _CFRONT
#include "options.h" // Work around CRAY template instantiation problem
#endif
#include "kernel.h"

#if(MCREAL)
typedef Real Mc;
typedef float McWeight;
#else
typedef Complex Mc;
typedef Complex McWeight;
#endif

extern Var *psibuffer,*psibufferR,*psibufferStop;
extern int reality; // Reality condition flag 

class GeometryBase {
protected:	
	int Nmode; // number of unreflected (explicitly evolved) modes
	int n; // total number of modes, including reflected modes
	int nindependent; // total number of independent modes
public:	
	int TotalNumber() {return n;}
	int IndependentNumber() {return nindependent;}
	
	virtual char *Name()=0;
	virtual int Valid(char *)=0;
	virtual void MakeBins()=0;
	virtual void List(ostream &)=0;
	virtual void ListTriads(ostream &) {}
	virtual void Initialize()=0;
	
	virtual Nu Linear(int)=0;
	virtual Real Forcing(int)=0;
	
	virtual Mode& ModeOf(int)=0;
	virtual Real Area(int)=0;
	
	Real X(int i) {return ModeOf(i).X();}
	Real Y(int i) {return ModeOf(i).Y();}
	Real K2(int i) {return ModeOf(i).K2();}
	Real K(int i) {return ModeOf(i).K();}
	Real Th(int i) {return ModeOf(i).Th();}
	
	virtual Real Normalization(int)=0;
	
	int Create() {
		MakeBins();

		if(verbose > 2) {
			cout.precision(3);
			List(cout);
			cout.precision(REAL_DIG);
		}
	
		Initialize();
		return Nmode;
	}
};

extern GeometryBase *Geometry;
Compare_t GeometryCompare;
KeyCompare_t GeometryKeyCompare;


#endif
