#include <iostream>
#include "xstream"

main()
{
	int N=255;
	
	fC << N << N << 1;

	oxstream fout("data");
	if(fout) 
	  for(j=0; j < N; j++) 
		for(i=0; i < N; i++) 
		    fout << (float) (i+j*N) << endl;
	fout.close();
}
