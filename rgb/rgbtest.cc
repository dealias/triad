#include <iostream>
#include "../xstream.h"

main()
{
	int Nx=255, Ny=255;

	oxstream fout("data");
	fout << Nx << Ny << 1;
	
	if(fout) 
	  for(unsigned j=0; j < N; j++) 
		for(unsigned i=0; i < N; i++) 
		    fout << (float) (i+j*N) << endl;
	fout.close();
}
