#include "kernel.h"
#include "Geometry.h"
#include "Basis.h"
#include "Cartesian.h"

char *Basis<Cartesian>::Name() {return "Cartesian";}

// Cartesian vocabulary
extern int Nx;
extern int Ny;

void Basis<Cartesian>::MakeBins()
{
	Cartesian *p;
	int i,j;
	
	SetPeriod(Cartesian(Nx,Ny));
	
	n=(2*period.x+1)*(2*period.y+1)-1;
	p=mode=new Cartesian[n];
	Nmode=(reality ? n/2 : n);
	nindependent=n/2;
	
	for(j=period.y; j >= 0; j--) // Evolved modes
		for(i=((j == 0) ? 1 : -period.x); i <= period.x; i++)
			*(p++)=Cartesian(i,j);
	
	for(j=-period.y; j <= 0; j++) // Reflected modes
		for(i=((j == 0) ? -1 : period.x); i >= -period.x; i--)
			*(p++)=Cartesian(i,j);
	
	if(p-mode != n) 
		msg(ERROR,"Calculated number and actual number of modes disagree."); 

	return;
}
