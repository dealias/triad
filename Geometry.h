#ifndef __Geometry_h__
#define __Geometry_h__ 1

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
	
	virtual Nu Linearity(int)=0;
	virtual Real Forcing(int)=0;
	
	virtual Real Area(int)=0;
	virtual Real K(int)=0;
	virtual Real K2(int)=0;
	virtual Real Th(int)=0;
	virtual Real X(int)=0;
	virtual Real Y(int)=0;
	
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
