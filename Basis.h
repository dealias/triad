#ifndef __Basis_h__
#define __Basis_h__ 1

#include "Geometry.h"
#include "DynVector.h"
#include "Linearity.h"

#define BASIS(key) {new Entry<Basis<key>,GeometryBase> (#key,GeometryTable);}

extern int nfft;
extern int Nxb,Nxb1,Nyb,Nyp;

extern Var *psix,*psiy,*vort;
extern Real *knorm2,*kfactor;

template<class T>
class Basis : public GeometryBase {
	T *mode; // pointer to table of modes
	T low; // lower limits of grid
	T high; // upper limits of grid
public:
	char *Name();
	int Valid(char *s) {
		return (strcmp(s,"Convolution")==0 || strcmp(s,"PS")==0);
	}

	void MakeBins();
	INLINE void List(ostream &);
	void Initialize();
	
	int InGrid(T &);
	
	Mode& ModeOf(int k) {return mode[k];}
	Real Area(int) {return 1.0;}
	
// Factor which converts |y|^2 to energy:
	Real Normalization(int);
	
	INLINE Nu Linear(int);
	INLINE Real Forcing(int);
};

template<class T>
INLINE void Basis<T>::List(ostream &os)
{
	os << "         " << Name() << " Mode Geometry:" << endl;
	for(int i=0; i < n; i++) os << mode[i] << newl;
	os << flush;
}

template <class T>
INLINE Nu Basis<T>::Linear(int i)
{
	Nu nu;
	Linearity->Evaluate(Geometry->ModeOf(i),nu);
	return nu;
}

void ForcingAt(int i, Real &force);

template <class T>
INLINE Real Basis<T>::Forcing(int i)
{
	Real force;
	ForcingAt(i,force);
	return force;
}

#endif
