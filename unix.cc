#ifdef _ANSI_C_SOURCE
#undef _ANSI_C_SOURCE
#endif

#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <pwd.h>
#include <sys/times.h>
#include <time.h>
#include <string.h>

#if _CRAY
#include <sys/mtimes.h>
static struct mtms buf;
static int parallel=-1;
#endif

extern char* run;
static const double init_time=time(NULL);

void cputime(double *cpu)
{
#if _CRAY
	if(parallel == -1) {
		parallel=strcmp(getenv("NCPUS"),"1");
		if(parallel) mtimes(&buf);
	}
	if(parallel) {
		cpu[0] = buf.mtms_mutime ? ((double) buf.mtms_mutime[0])/CLK_TCK : 0.0;
		cpu[1] = 0.0;
		cpu[2] = 0.0;
		return;
	}
#endif		
	struct tms buf;
	times(&buf);
	cpu[0] = ((double) buf.tms_utime)/CLK_TCK;
	cpu[1] = ((double) buf.tms_cutime)/CLK_TCK;
	cpu[2] = ((double) (buf.tms_stime+buf.tms_cstime))/CLK_TCK;
}

// Don't notify user about runs shorter than this many seconds.
static const double longrun=500.0;

void mailuser(char *text)
{
	if(time(NULL)-init_time < longrun) return;
	char *user=getpwuid(getuid())->pw_name;
	char *buf=new char[50+strlen(run),strlen(text)+strlen(user)];
	sprintf(buf,"mail -s 'Run %s %s.' %s < /dev/null > /dev/null",
			run,text,user);
	system(buf);
}
