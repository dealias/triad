#include <iostream.h>
#if _CRAY
#include <sys/machd.h>
#endif

extern "C" TREMAIN(const double& seconds);

extern double polltime;

static double tlimit=0.0;
static double last_seconds=0.0;
static double cpu[ncputime];

int poll()
{
	double seconds=0.0;
	if(!tlimit) TREMAIN(tlimit);
	
	cputime(cpu);
	for(i=0; i < ncputime; i++) seconds += cpu[i];

	if(2*seconds-last_seconds < 0.95*tlimit) {
		last_seconds=seconds;
		return 0;
	}
	cout << endl << "CPU time limit is approaching; exit forced." << endl;
	return 1;
}
