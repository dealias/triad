#include <iostream>
#include "../xstream.h"

main()
{
	int Nx=255, Ny=255;

	oxstream fout("data");
	fout << Nx << Ny << 1;
	
	if(fout) 
	  for(unsigned j=0; j < Ny; j++) 
		for(unsigned i=0; i < Nx; i++) 
		    fout << (float) (i+Nx*j) << endl;
	fout.close();
}
