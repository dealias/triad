#include "options.h"
#include "kernel.h"
#include "Geometry.h"
#include "Basis.h"
#include "Cartesian.h"

char *Basis<Cartesian>::Name() {return "Cartesian";}

// Cartesian vocabulary
int Nx=17; // Number of modes in x-direction
int Ny=17; // Number of modes in y-direction

int Nx0,NRows,NPad,NPadTop,xoffset;
unsigned int log2Nxb,log2Nyb;
int nfft; // Total number of FFT elements;
int Nxb,Nyb,Nyp;
Cartesian *CartesianMode;

void Basis<Cartesian>::MakeBins()
{
	Cartesian *p;
	int i,j;
	
	if(reality && ((Nx/2)*2 == Nx || (Ny/2)*2 == Ny))
			msg(ERROR,"The reality condition requires that Nx and Ny be odd");
	
	int halfNx=Nx/2;
	int halfNy=Ny/2;
	low=Cartesian(halfNx-Nx+1,halfNy-Ny+1);
	high=Cartesian(halfNx,halfNy);
		
	n=Nx*Ny-1;
	p=mode=CartesianMode=new Cartesian[n];
	Nmode=(reality ? n/2 : n);
	nindependent=(reality || n % 2) ? Nmode : n/2;
	
	NRows=high.Row()+1;
	
	for(j=0; j <= high.Row(); j++) { // Evolved modes
		for(i=((j == 0) ? 1 : low.Column()); i <= high.Column(); i++)
			*(p++)=Cartesian(i,j);
		}
	
	int nminx=(3*Nx-1)/2;
	int nminy=(3*Ny-1)/2;
	for(log2Nxb=0; nminx > (1 << log2Nxb); log2Nxb++);
	for(log2Nyb=0; nminy > (1 << log2Nyb); log2Nyb++);
	Nxb=1 << log2Nxb;
	Nyb=1 << log2Nyb;
	Nyp=(Nyb/2+1);
	xoffset=Nxb/2;
	int Nxb1=Nxb;
	
#if _CRAY // Avoid memory bank conflicts
	Nxb1 += 1;
#endif	
	nfft=Nxb1*Nyp;
	Nx0=(Nx-1)/2;
	NPad=Nxb1-Nx;
	NPadTop=(Nyp-(Ny+1)/2)*Nxb1+Nxb1-((Nx+1)/2+xoffset);
	
	for(j=0; j >= low.Row(); j--) // Reflected modes
		for(i=((j == 0) ? -1 : high.Column()); i >= low.Column(); i--)
			*(p++)=Cartesian(i,j);
	
	if(p-mode != n) 
		msg(ERROR,"Calculated number and actual number of modes disagree."); 

	return;
}

#if _CRAY
void CartesianPad(Var *to_, Var * from)
{
	Var *to=to_+xoffset;
	*(to++)=0.0;
	Var *tostop=to+Nx0;
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
