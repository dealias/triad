extern "C" TREMAIN(const double& seconds);

extern double polltime;

int poll()
{
	double seconds;
	TREMAIN(seconds);
	if(seconds > 2.0*polltime) return 1;
	return 0; // Process is running out of CPU time; force a graceful exit.
}
