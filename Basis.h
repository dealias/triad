#ifndef __Basis_h__
#define __Basis_h__ 1

#include "Geometry.h"
#include "DynVector.h"

#define BASIS(key) {new Entry<Basis<key>,GeometryBase> (#key,GeometryTable);}

extern int nfft;
extern int Nxb,Nyb,Nyp;

extern Var *psix,*psiy,*vort;
extern Real *kfactor;

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
	INLINE void Initialize();
	
	int InGrid(T &);
	Real Area(int) {return 1.0;}
	
	Real K(int k) {return mode[k].K();}
	Real K2(int k) {return mode[k].K2();}
	Real Th(int k) {return mode[k].Th();}
	Real X(int k) {return mode[k].X();}
	Real Y(int k) {return mode[k].Y();}
	
// Factor which converts |y|^2 to energy:
	Real Normalization(int);
	
	inline Nu Linearity(int);
	inline Real Forcing(int);
};

template<class T>
INLINE void Basis<T>::List(ostream &os)
{
	os << "         " << Name() << " Mode Geometry:" << endl;
	for(int i=0; i < n; i++) os << mode[i] << newl;
	os << flush;
}

template<class T>
INLINE void Basis<T>::Initialize()
{
	kfactor=new Real[Nmode];
	
	if(strcmp(Problem->Abbrev(),"PS") == 0) {
		cout << endl << "ALLOCATING FFT BUFFERS (" << Nxb << "X" << Nyp << ")."
			 << endl;
		psix=new Var[nfft];
		psiy=new Var[nfft];
		vort=new Var[nfft];
		Real scale=Nxb*Nyb;
		for(int k=0; k < Nmode; k++) kfactor[k]=1.0/(scale*mode[k].K2());
	} else {
		psibuffer=new Var[n];
		psibufferR=(reality ? psibuffer+Nmode : psibuffer);
		for(int k=0; k < Nmode; k++) kfactor[k]=1.0/mode[k].K2();
	}
}

void LinearityAt(int i, Nu& nu);

template <class T>
INLINE Nu Basis<T>::Linearity(int i)
{
	Nu nu;
	LinearityAt(i,nu);
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
