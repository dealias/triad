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
#include <strstream.h>
#include <malloc.h>
#include <sys/utsname.h>

#if _CRAYMVP
#include <sys/types.h>
#include <sys/jtab.h>
#endif

extern char* run;
static const double init_time=time(NULL);

void cputime(double *cpu)
{
#if _CRAYMVP
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
static char mail_cmd[]="mailx";
#else
static char mail_cmd[]="mail";
#endif

// Mail a message to the user.
void mailuser(char *text)
{
	if(time(NULL)-init_time < longrun) return;
	strstream buf;
	buf << mail_cmd << " -s 'Run " << run << " " << text << ".' "
		<< getpwuid(getuid())->pw_name << " < /dev/null > /dev/null" << ends;
	system(buf.str());
}

static char rmrf_cmd[]="rm -rf";

// Recursively delete a directory and all of its contents.
void remove_dir(char *text)
{
	strstream buf;
	buf << rmrf_cmd << " " << text << ends;
	system(buf.str());
}

extern "C" int getdomainname(char *name, size_t len);

struct utsname platform;

const int ndomain=65;
char domain[ndomain];

char *machine()
{
	uname(&platform);
	getdomainname(domain,ndomain);
	if(strcmp(domain,"(none)") == 0) *domain=0;
	strstream buf;
	buf << platform.nodename << ((*domain) ? "." : "") << domain << " ("
		<< platform.machine << ")" << ends;
	buf.rdbuf()->freeze();
	return buf.str();
}
