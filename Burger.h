#ifndef __Burger_h__
#define __Burger_h__ 1

#include "kernel.h"
#include "Geometry.h"
#include "Partition.h"
#include "Linearity.h"
#include "Basis.h"

extern int Nmode;   // number of explictly evolved modes
extern int NmodeR;  // number of reflected modes
extern int Ntotal; // total number of (evolved+reflected) modes

class BurgerVocabulary : public VocabularyBase {
public:
	char *Name();
	char *Abbrev();
	BurgerVocabulary();
	Table<GeometryBase> *GeometryTable;
	Table<LinearityBase> *LinearityTable;
	GeometryBase *NewGeometry(char *& key) {
		return GeometryTable->Locate(key);
	}
	LinearityBase *NewLinearity(char *& key) {
		return LinearityTable->Locate(key);
	}
};

class Burger : public ProblemBase {
public:
	Burger() {}
	void InitialConditions();
	void Initialize();
	void OpenOutput();
	void Output(int it);
	void FinalOutput();
	void LinearSrc(Var *, Var *, double);
};

class PS : public Burger {
public:
	PS() {
		if(!reality) msg(ERROR,"Pseudospectral approximation needs reality=1");
	}
	char *Name() {return "Pseudospectral";}
	void NonLinearSrc(Var *source, Var *psi, double);
};



#endif
