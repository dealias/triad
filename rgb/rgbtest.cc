#include <iostream>
#include "../xstream.h"

main()
{
	int Nx=256, Ny=256;

	oxstream fout("data");
	fout << Nx << Ny << 1;
	
	if(fout) 
	  for(int j=Ny-1; j >= 0; j--) {
	    for(unsigned j=0; j < Ny; j++) {
	    for(unsigned i=0; i < Nx; i++) {
	      fout << (float) (i+Nx*j);
	  }
	fout.close();
}
