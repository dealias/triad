#include <iostream.h>
#if _CRAY
#include <sys/machd.h>
#endif
#include "utils.h"

extern "C" TREMAIN(const double& seconds);

extern double polltime;
extern int restart;

static double tlimit=0.0;
static double last_seconds=0.0;
static double cpu[ncputime];

int poll()
{
	double seconds=0.0;
	if(!tlimit) TREMAIN(tlimit);
	
	cputime(cpu);
	for(int i=0; i < ncputime; i++) seconds += cpu[i];

	if(2.0*seconds-last_seconds < (restart ? 0.90 : 0.95)*tlimit) {
		last_seconds=seconds;
		return 0;
	}
	cout << endl << "CPU time limit is approaching; exit forced." << endl;
	return 1;
}
