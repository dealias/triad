#include <iostream.h>
#include <sys/machd.h>

extern "C" TREMAIN(const double& seconds);
extern "C" SECOND(const double& seconds);

extern double polltime;

static double tlimit=0.0;

int poll()
{
	double seconds;
	if(!tlimit) TREMAIN(tlimit);
	SECOND(seconds);
	if(seconds+NCPU*polltime < tlimit) return 0;
	cout << endl << "CPU time limit is approaching; exit forced." << endl;
	return 1;
}
