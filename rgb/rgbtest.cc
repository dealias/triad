#include <iostream>
#include "../xstream.h"

main()
{
	int Nx=256, Ny=256;

	oxstream fout("data");
	fout << Nx << Ny << 1;
	
	if(fout) 
	  for(unsigned j=Ny-1; j >= 0; j--) 
	    for(unsigned i=0; i < Nx; i++) {
	      fout << (float) (i+Nx*j);
	      std::cout << (float) (i+Nx*j) << std::endl;
	    }
	fout.close();
}
