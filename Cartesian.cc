#include "kernel.h"
#include "Geometry.h"
#include "Basis.h"
#include "Cartesian.h"

char *Basis<Cartesian>::Name() {return "Cartesian";}

// Cartesian vocabulary
int Nx=17; // Number of modes in x-direction
int Ny=17; // Number of modes in y-direction

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
	
	for(j=0; j <= high.Row(); j++) // Evolved modes
		for(i=((j == 0) ? 1 : low.Column()); i <= high.Column(); i++)
			*(p++)=Cartesian(i,j);
	
	for(j=0; j >= low.Row(); j--) // Reflected modes
		for(i=((j == 0) ? -1 : high.Column()); i >= low.Column(); i--)
			*(p++)=Cartesian(i,j);
	
	if(p-mode != n) 
		msg(ERROR,"Calculated number and actual number of modes disagree."); 

	return;
}
