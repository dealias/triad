#ifdef _ANSI_C_SOURCE
#undef _ANSI_C_SOURCE
#endif

#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <pwd.h>
#if(_CRAY)
#include <sys/mtimes.h>
#else
#include <sys/times.h>
#endif
#include <time.h>
#include <string.h>

extern char* run;
static const double init_time=time(NULL);

#if(_CRAY)
void cputime(double *cpu)
{
	struct mtms buf;
	mtimes(&buf);
	short ncpu=buf.mtms_conn;
	time_t child_mutime=0.0;
	
	for(int i=1; i < ncpu; i++) child_mutime += buf.mtms_mutime[i];
	
	cpu[0] = ((double) buf.mtms_mutime[0])/CLK_TCK;
	cpu[1] = ((double) child_mutime)/CLK_TCK;
	cpu[2] = 0.0;
}
#else
void cputime(double *cpu)
{
	struct tms buf;
	times(&buf);
	cpu[0] = ((double) buf.tms_utime)/CLK_TCK;
	cpu[1] = ((double) buf.tms_cutime)/CLK_TCK;
	cpu[2] = ((double) (buf.tms_stime+buf.tms_cstime))/CLK_TCK;
}
#endif

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
