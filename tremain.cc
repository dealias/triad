#include <iostream.h>
#if _CRAY
#include <sys/machd.h>
#endif

extern "C" TREMAIN(const double& seconds);
extern "C" SECOND(const double& seconds);

extern double polltime;

static double tlimit=0.0;
static double last_seconds=0.0;

int poll()
{
	double seconds;
	if(!tlimit) TREMAIN(tlimit);
	SECOND(seconds);
	if(2*seconds-last_seconds < 0.95tlimit) {
		last_seconds=seconds;
		return 0;
	}
	cout << endl << "CPU time limit is approaching; exit forced." << endl;
	return 1;
}
