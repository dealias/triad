#include "options.h"
#include "kernel.h"
#include "Geometry.h"
#include "Basis.h"
#include "Cartesian.h"

char *Basis<Cartesian>::Name() {return "Cartesian";}

// Cartesian vocabulary
int Nx=17; // Number of modes in x-direction
int Ny=17; // Number of modes in y-direction
int *ModeBin;

int Nx0,NRows,NPad,NPadTop,xoffset;
int Nevolved;
unsigned int log2Nxb,log2Nyb;
int nfft; // Total number of FFT elements;
int Nxb,Nxb1,Nyb,Nyp;
Cartesian *CartesianMode;

Real krmin=1.0;
Real krmin2=0.0;

void Basis<Cartesian>::MakeBins()
{
	int i,j;
	
	if(reality && ((Nx/2)*2 == Nx || (Ny/2)*2 == Ny))
			msg(ERROR,"The reality condition requires that Nx and Ny be odd");
	
	int halfNx=Nx/2;
	int halfNy=Ny/2;
	low=Cartesian(halfNx-Nx+1,halfNy-Ny+1);
	high=Cartesian(halfNx,halfNy);
		
	n=Nx*Ny-1;
	mode=CartesianMode=new Cartesian[n];
	Nmode=(reality ? n/2 : n);
	nindependent=(reality || n % 2) ? Nmode : n/2;
	
	NRows=high.Row()+1;
	Nevolved=high.Row()*Nx+high.Column();
	
	Cartesian mode0=Cartesian(0,0);
	for(i=0; i < n; i++) mode[i]=mode0;
		
	for(j=low.Row(); j <= high.Row(); j++) {
		for(i=low.Column(); i <= high.Column(); i++) {
			Cartesian m(i,j);
			if(m != mode0) mode[m.ModeIndex()]=m;
		}
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
	int nminy=(3*Ny-1)/2;
	for(log2Nxb=0; nminx > (1 << log2Nxb); log2Nxb++);
	for(log2Nyb=0; nminy > (1 << log2Nyb); log2Nyb++);
	Nxb=1 << log2Nxb;
	Nyb=1 << log2Nyb;
	Nyp=Nyb/2+1;
	xoffset=Nxb/2;
	Nxb1=Nxb;
	
#if _CRAY // Avoid memory bank conflicts
	Nxb1 += 1;
#endif	
	nfft=Nxb1*Nyp;
	Nx0=(Nx-1)/2;
	NPad=Nxb1-Nx;
	NPadTop=(Nyp-(Ny+1)/2)*Nxb1+Nxb1-((Nx+1)/2+xoffset);
}

void DiscretePad(Var *to, Var *from, Real *norm)
{
	int k=0;
	Var *tostop=to+xoffset-Nx0;
#pragma ivdep		
	for(; to < tostop; to++) *to=0.0;
	tostop += Nx0;
	Var *s=from+Nx0;
#pragma ivdep		
	for(; to < tostop; to++) *to=conj(*(--s));
	*(to++)=0.0;
	tostop=to+Nx0;
#pragma ivdep		
	for(; to < tostop; to++) {
		int index=ModeBin[k++];
		if(index >= 0) *to=from[index]*norm[index];
		else *to=0.0;
	}
    for(int j=0; j < NRows-1; j++) {
        tostop += NPad;
#pragma ivdep		
        for(; to < tostop; to++) *to=0.0;
        tostop += Nx;
#pragma ivdep		
		for(; to < tostop; to++) {
			int index=ModeBin[k++];
			if(index >= 0) *to=from[index]*norm[index];
			else *to=0.0;
		}
    }
	tostop += NPadTop;
#pragma ivdep		
	for(; to < tostop; to++) *to=0.0;
}

#if _CRAYMVP
void CartesianPad(Var *to_, Var *from)
{
	Var *to=to_;
	Var *tostop=to+xoffset-Nx0;
#pragma ivdep		
	for(; to < tostop; to++) *to=0.0;
	tostop += Nx0;
	const Var *s=from+Nx0;
#pragma ivdep		
	for(; to < tostop; to++) *to=conj(*(--s));
	*(to++)=0.0;
	tostop=to+Nx0;
#pragma ivdep		
	for(; to < tostop; to++) *to=*(from++);
    for(int j=0; j < NRows-1; j++) {
        tostop += NPad;
#pragma ivdep		
        for(; to < tostop; to++) *to=0.0;
        tostop += Nx;
#pragma ivdep		
		for(; to < tostop; to++) *to=*(from++);
    }
	tostop += NPadTop;
#pragma ivdep		
	for(; to < tostop; to++) *to=0.0;
}

void CartesianUnPad(Var *to_, Var *from)
{
	from += xoffset+1;
	Var *to=to_;
	Var *tostop=to+Nx0;
#pragma ivdep		
	for(; to < tostop; to++) *to=*(from++);
    for(int j=1; j < NRows; j++) {
		from += NPad;
        tostop += Nx;
#pragma ivdep		
		for(; to < tostop; to++) *to=*(from++);
    }
}
#endif
