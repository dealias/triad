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
#include <sys/jtab.h>
#endif

extern char* run;
static const double init_time=time(NULL);

void cputime(double *cpu)
{
#if _CRAY
	struct jtab jbuf;
	getjtab(&jbuf);
	
	cpu[0] = ((double) jbuf.j_ucputime)/CLK_TCK;
	cpu[1] = 0.0;
	cpu[2] = ((double) jbuf.j_scputime)/CLK_TCK;
#else
	struct tms buf;
	times(&buf);
	cpu[0] = ((double) buf.tms_utime)/CLK_TCK;
	cpu[1] = ((double) buf.tms_cutime)/CLK_TCK;
	cpu[2] = ((double) (buf.tms_stime+buf.tms_cstime))/CLK_TCK;
#endif
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
