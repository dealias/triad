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
	INLINE void Initialize();
	
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

template<class T>
INLINE void Basis<T>::Initialize()
{
	knorm2=new Real[Nmode];
	kfactor=new Real[Nmode];
	
	if(strcmp(Problem->Abbrev(),"PS") == 0) {
		cout << endl << "ALLOCATING FFT BUFFERS (" << Nxb << " x " << Nyp
			 << ")." << endl;
		psix=new Var[nfft];
		psiy=new Var[nfft];
		vort=new Var[nfft];
		Real scale=Nxb*Nyb;
		for(int k=0; k < Nmode; k++) {
			knorm2[k]=mode[k].K2();
			kfactor[k]=1.0/(scale*Normalization(k));
		}
	} else {
		psibuffer=new Var[n];
		psibufferR=(reality ? psibuffer+Nmode : psibuffer);
		for(int k=0; k < Nmode; k++) 
			kfactor[k]=1.0/Normalization(k);
	}
}

template <class T>
INLINE Nu Basis<T>::Linear(int i)
{
	Nu nu;
	Linearity->Evaluate(Polar(Geometry->K(i),Geometry->Th(i)),nu);
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
