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
#include <sys/types.h>
#include <sys/mtimes.h>
#include <sys/machd.h>
static struct mtms mbuf;
static int parallel=-1;
#endif

extern char* run;
static const double init_time=time(NULL);

void cputime(double *cpu)
{
#if _CRAY
	if(parallel == -1) parallel=strcmp(getenv("NCPUS"),"1");
	if(parallel) {
		time_t task=0,child=0;
		
		mtimes(&mbuf);
		double update=0;
		while(update != mbuf.mtms_update) {
			update=mbuf.mtms_update;
			for(int i=0; i < NCPU; i++) {
				task += mbuf.mtms_mutime[i];
				child += i*mbuf.mtms_mutime[i];
			}
		}
		
		cpu[0] = ((double) task)/CLK_TCK;
		cpu[1] = ((double) child)/CLK_TCK;
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

#if _CRAY
static char command[]="mailx";
#else
static char command[]="mail";
#endif

void mailuser(char *text)
{
	if(time(NULL)-init_time < longrun) return;
	char *user=getpwuid(getuid())->pw_name;
	char *buf=new char[50+strlen(run),strlen(text)+strlen(user)];
	sprintf(buf,"%s -s 'Run %s %s.' %s < /dev/null > /dev/null",command,
			run,text,user);
	system(buf);
}
