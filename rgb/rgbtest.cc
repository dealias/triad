#include <iostream>
#include "xstream"

main()
{
	int N=255;

	oxstream fout("data");
	fout << N << N << 1;
	
	if(fout) 
	  for(j=0; j < N; j++) 
		for(i=0; i < N; i++) 
		    fout << (float) (i+j*N) << endl;
	fout.close();
}
