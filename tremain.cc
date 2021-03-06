#include <iostream.h>
#include "utils.h"

extern "C" TREMAIN(const double& seconds);

extern double polltime;

static double tlimit=0.0;
static double last_seconds=0.0;
static double cpu[ncputime];

void poll()
{
  double seconds=0.0;
  if(!tlimit) TREMAIN(tlimit);
	
  cputime(cpu);
  for(int i=0; i < ncputime; i++) seconds += cpu[i];

  if(2.0*seconds-last_seconds < 0.95*tlimit) {
    last_seconds=seconds;
    return 0;
  }
  msg(WARNING_GLOBAL, "CPU time limit is approaching; exit forced");
  exit(CONTINUE);
}
