#include <iostream>
#include "../xstream.h"

main()
{
	int Nx=255, Ny=255;

	oxstream fout("data");
	fout << Nx << Ny << 1;
	
	if(fout) 
	  for(unsigned j=Ny-1; j >= 1; j--) 
	    for(unsigned i=0; i < Nx; i++) {
	      fout << (float) (i+Nx*j);
	      std::cout << (float) (i+Nx*j) << std::endl;
	    }
	fout.close();
}
