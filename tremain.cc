extern "C" TREMAIN(const double& seconds);

extern double polltime;

int poll()
{
	double seconds;
	TREMAIN(seconds);
	if(seconds > 10.0*polltime) return 0;
	return 1; // Process is running out of CPU time; force a graceful exit.
}
