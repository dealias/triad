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

#if _CRAY && !_CRAYMPP
#include <sys/types.h>
#include <sys/jtab.h>
#endif

extern char* run;
static const double init_time=time(NULL);

void cputime(double *cpu)
{
#if _CRAY && !_CRAYMPP
	struct jtab jbuf;
	getjtab(&jbuf);
	
	cpu[0] = ((double) jbuf.j_ucputime)/CLOCKS_PER_SEC;
	cpu[1] = 0.0;
	cpu[2] = ((double) jbuf.j_scputime)/CLOCKS_PER_SEC;
#else
#if 0	
	struct tms buf;
	times(&buf);
	cpu[0] = ((double) buf.tms_utime)/CLOCKS_PER_SEC;
	cpu[1] = ((double) buf.tms_cutime)/CLOCKS_PER_SEC;
	cpu[2] = ((double) (buf.tms_stime+buf.tms_cstime))/CLOCKS_PER_SEC;
#else	
	static const double clockinv=1.0/CLOCKS_PER_SEC;
	cpu[0]=clock()*clockinv;
	cpu[1]=cpu[2]=0.0;
#endif
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
void mailuser(const char *text)
{
	if(time(NULL)-init_time < longrun) return;
	strstream buf;
	buf << mail_cmd << " -s 'Run " << run << " " << text << ".' "
		<< getpwuid(getuid())->pw_name << " < /dev/null > /dev/null" << ends;
	system(buf.str());
}

static char rmrf_cmd[]="rm -rf";

// Recursively delete a directory and all of its contents.
void remove_dir(const char *text)
{
	strstream buf;
	buf << rmrf_cmd << " " << text << ends;
	system(buf.str());
}

static char copy_cmd[]="cp -rf";

// Copy a file or directory to a new location
int copy(const char *oldname, const char *newname)
{
	strstream buf;
	buf << copy_cmd << " " << oldname << " " << newname << ends;
	return system(buf.str());
}

extern "C" int getdomainname(char *name, size_t len);

static struct utsname platform;
static const int ndomain=65;
static char domain[ndomain];

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

const int ndate=64;
static char time_date[ndate]="";

char *date()
{
	const time_t bintime=time(NULL);
	strftime(time_date,ndate,"%a %b %d %H:%M:%S %Z %Y",localtime(&bintime));
	return time_date;
}

char *tempdir()
{
	return getenv("TMP_DIR");
}
