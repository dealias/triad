#ifndef __Geometry_h__
#define __Geometry_h__ 1

#include "kernel.h"
#include "utils.h"

extern int reality; // Reality condition flag 
extern Var *psibuffer,*psibuffer0,*psibufferR,*psibufferStop;
extern Var *convolution,*convolution0;
extern int pseudospectral;
extern unsigned int log2n; // Number of FFT levels

class GeometryBase {
protected:	
	int Nmode; // number of unreflected (explicitly evolved) modes
	int n; // total number of modes, including reflected modes
	int nindependent; // total number of independent modes
public:	
	int TotalNumber() {return n;}
	int IndependentNumber() {return nindependent;}
	
	virtual char *Name()=0;
	virtual int ValidApproximation(char *)=0;
	virtual void MakeBins()=0;
	virtual void List(ostream &)=0;
	virtual void ListTriads()=0;
	virtual void ComputeTriads()=0;
	
	virtual Nu Linearity(int)=0;
	
	virtual Real Area(int)=0;
	virtual Real K(int)=0;
	virtual Real K2(int)=0;
	virtual Real Th(int)=0;
	virtual Real Kx(int)=0;
	virtual Real Ky(int)=0;
	
	virtual Real Normalization(int)=0;
	
	int Create() {
		MakeBins();

		if(pseudospectral) {
			for(log2n=0; Nmode+1 > (1 << log2n)/3; log2n++);
			int n2=1 << (log2n-1);
			cout << n2 << " FFT COMPONENTS ALLOCATED." << endl;
			convolution0=new Var[n2];
			convolution=convolution0+1;
			psibuffer0=new Var[n2];
			psibuffer=psibuffer0+1;
		} else {
			psibuffer=new Var[n];
			psibufferR=(reality ? psibuffer+Nmode : psibuffer);
			psibufferStop=psibuffer+n;
		}
		
		if(verbose > 2) {
			cout.precision(3);
			List(cout);
			cout.precision(REAL_DIG);
		}
	
		ComputeTriads();
		if(verbose > 2) ListTriads();
		return Nmode;
	}
};

extern GeometryBase *Geometry;
Compare_t GeometryCompare;
KeyCompare_t GeometryKeyCompare;


#endif
