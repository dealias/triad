#include "kernel.h"
#include "Geometry.h"
#include "Basis.h"
#include "Cartesian.h"

char *Basis<Cartesian>::Name() {return "Cartesian";}

// Cartesian vocabulary
int Nx=17; // Number of modes in x-direction
int Ny=17; // Number of modes in y-direction

int NRows,NPad;
int *RowBoundary;
Var *ZeroBuffer; 
unsigned int log2n; // Number of FFT levels
int Npsibuffer;

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
	RowBoundary=new int[NRows+1];
	
	for(j=0; j <= high.Row(); j++) { // Evolved modes
		RowBoundary[j]=p-mode;
		for(i=((j == 0) ? 1 : low.Column()); i <= high.Column(); i++)
			*(p++)=Cartesian(i,j);
		}
	
	RowBoundary[NRows]=p-mode;
	NPad=RowBoundary[1]-RowBoundary[0];
	ZeroBuffer=new Var[NPad];
	for(i=0; i < NPad; i++) ZeroBuffer[i]=0.0;
	Npsibuffer=Nmode+NRows*NPad;
	int ntotal=(3*Nx-1)*(3*Ny-1)/4;
	for(log2n=0; ntotal > (1 << log2n); log2n++);
	
	for(j=0; j >= low.Row(); j--) // Reflected modes
		for(i=((j == 0) ? -1 : high.Column()); i >= low.Column(); i--)
			*(p++)=Cartesian(i,j);
	
	if(p-mode != n) 
		msg(ERROR,"Calculated number and actual number of modes disagree."); 

	return;
}

#if _CRAY
void CartesianPad(Var * restrict to_, Var * from)
{
	Var *to=to_;
    for(int i=0; i < NRows; i++) {
        Var *tostop=to+RowBoundary[i+1]-RowBoundary[i];
#pragma ivdep		
        for(; to < tostop; to++) *to=*(from++);
        tostop += NPad;
        Var *zero=ZeroBuffer;
#pragma ivdep		
        for(; to < tostop; to++) *to=*(zero++);
    }
}

void CartesianUnPad(Var * restrict to_, Var * from)
{
	Var *to=to_;
    for(int i=0; i < NRows; i++) {
        int ncol=RowBoundary[i+1]-RowBoundary[i];
        Var *tostop=to+ncol;
#pragma ivdep		
        for(; to < tostop; to++) *to=*(from++);
        from += NPad;
    }
}
#endif
