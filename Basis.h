#ifndef __Basis_h__
#define __Basis_h__ 1

#include "Geometry.h"
#include "DynVector.h"

#define BASIS(key) {new Entry<Basis<key>,GeometryBase> (#key,GeometryTable);}

extern Real *kinv2;

template<class T>
class Basis : public GeometryBase {
	T *mode; // pointer to table of modes
	T low; // lower limits of grid
	T high; // upper limits of grid
public:
	char *Name();
	int ValidApproximation(char *s) {
		return (strcmp(s,"Convolution")==0 || strcmp(s,"PS")==0);
	}

	void MakeBins();
	void List(ostream &);
	void Initialize();
	
	int InGrid(T &);
	Real Area(int) {return 1.0;}
	
	Real K(int k) {return mode[k].K();}
	Real K2(int k) {return mode[k].K2();}
	Real Th(int k) {return mode[k].Th();}
	Real X(int k) {return mode[k].X();}
	Real Y(int k) {return mode[k].Y();}
	
// Factor which converts |y|^2 to energy:
	Real Normalization(int);
	
	Nu Linearity(int);
	Real Forcing(int);
};


template<class T>
void Basis<T>::List(ostream &os)
{
	os << "         " << Name() << " Mode Geometry:" << endl;
	for(int i=0; i < n; i++) os << mode[i] << endl;
	cout << endl;
}
template<class T>
void Basis<T>::Initialize()
{
	int n2=1;
	for(log2n=1; Nmode+1 > 4*n2/9; log2n++, n2 *= 2);
	cout << n2 << " FFT COMPONENTS ALLOCATED." << endl;
	convolution0=new Var[n2+1];
	convolution=convolution0+1;
	psibuffer0=new Var[n2+1];
	psibuffer=psibuffer0+1;
	psitemp=new Var[Nmode];
			
	kinv2=new Real[Nmode];
	for(int k=0; k < Nmode; k++) kinv2[k]=1.0/mode[k].K2();
}

void LinearityAt(int i, Nu& nu);

template <class T>
Nu Basis<T>::Linearity(int i)
{
	Nu nu;
	LinearityAt(i,nu);
	return nu;
}

void ForcingAt(int i, Real &force);

template <class T>
Real Basis<T>::Forcing(int i)
{
	Real force;
	ForcingAt(i,force);
	return force;
}

#endif
