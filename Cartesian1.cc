#include "options.h"
#include "kernel.h"
#include "Geometry.h"
#include "Basis.h"
#include "Cartesian1.h"

char *Basis<Cartesian1>::Name() {return "Cartesian1";}

// Cartesian1 vocabulary
int Nx=17; // Number of modes in x-direction
int *ModeBin;

int Nx0,NPad,xoffset;
int Nevolved;
unsigned int log2Nxb;
int nfft; // Total number of FFT elements;
int Nxb;
Cartesian1 *Cartesian1Mode;

Real krmin=1.0;
Real krmin2=0.0;

void Basis<Cartesian1>::MakeBins()
{
	int i;
	
	if(reality && ((Nx/2)*2 == Nx))
			msg(ERROR,"The reality condition requires that Nx be odd");
	
	int halfNx=Nx/2;
	low=Cartesian1(halfNx-Nx+1);
	high=Cartesian1(halfNx);
		
	n=Nx-1;
	mode=Cartesian1Mode=new Cartesian1[n];
	Nmode=(reality ? n/2 : n);
	nindependent=(reality || n % 2) ? Nmode : n/2;
	
	Nevolved=high.Column();
	
	Cartesian1 mode0=Cartesian1(0);
	for(i=0; i < n; i++) mode[i]=mode0;
		
	for(i=low.Column(); i <= high.Column(); i++) {
		Cartesian1 m(i);
		if(m != mode0) mode[m.ModeIndex()]=m;
	}
	
	for(i=0; i < n; i++) {
		if(mode[i] == mode0) msg(ERROR,"Zero mode (%d) encountered",i);
	}
	
	set_fft_parameters();

	return;
}

void set_fft_parameters()
{	
	int nminx=(3*Nx-1)/2;
	for(log2Nxb=0; nminx > (1 << log2Nxb); log2Nxb++);
	Nxb=1 << log2Nxb;
	xoffset=Nxb/2;
	nfft=Nxb;
	Nx0=(Nx-1)/2;
	NPad=Nxb-Nx;
}

