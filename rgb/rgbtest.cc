#include <iostream>
#include "../xstream.h"

main()
{
	int N=255;

	oxstream fout("data");
	fout << N << N << 1;
	
	if(fout) 
	  for(unsigned j=0; j < N; j++) 
		for(unsigned i=0; i < N; i++) 
		    fout << (float) (i+j*N) << endl;
	fout.close();
}
