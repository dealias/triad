#ifdef sun
typedef void stack_t;
#endif

#ifdef _ANSI_C_SOURCE
#undef _ANSI_C_SOURCE
#endif

#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include <pwd.h>
#include <sys/times.h>

#ifdef __STRICT_ANSI__
#undef __STRICT_ANSI__
#endif
#include <time.h>
#define __STRICT_ANSI__

#include <string>
#include <sstream>
#include <cstdlib>
#include <sys/utsname.h>

#if defined(_CRAY) && !defined(_CRAYMPP)
#include <sys/types.h>
#include <sys/jtab.h>
#endif

#include "utils.h"

using std::ostringstream;

extern char* run;
static const double init_time=time(NULL);
static const double ticktime=1.0/sysconf(_SC_CLK_TCK);

void cputime(double *cpu)
{
#if defined(_CRAY) && !defined(_CRAYMPP)
  struct jtab jbuf;
  getjtab(&jbuf);
	
  cpu[0] = ((double) jbuf.j_ucputime)*ticktime;
  cpu[1] = 0.0;
  cpu[2] = ((double) jbuf.j_scputime)*ticktime;
#else
  struct tms buf;
  times(&buf);
  cpu[0] = ((double) buf.tms_utime)*ticktime;
  cpu[1] = ((double) buf.tms_cutime)*ticktime;
  cpu[2] = ((double) (buf.tms_stime+buf.tms_cstime))*ticktime;
#endif
}

// Don't notify user about runs shorter than this many seconds.
static const double longrun=500.0;

#if defined(_CRAY)
static char mail_cmd[]="mailx";
#else
static char mail_cmd[]="mail";
#endif

// Mail a message to the user.
void mailuser(const char *text)
{
  if(time(NULL)-init_time < longrun) return;
  ostringstream buf;
  buf << mail_cmd << " -s 'Run " << run << " " << text << ".' "
      << getpwuid(getuid())->pw_name << " < /dev/null > /dev/null" << ends;
  system(buf.str().c_str());
}

static char rmrf_cmd[]="rm -rf";

// Recursively delete a directory and all of its contents.
void remove_dir(const char *text)
{
  ostringstream buf;
  buf << rmrf_cmd << " " << text << ends;
  system(buf.str().c_str());
}

static char copy_cmd[]="cp -rf";

// Copy a file or directory to a new location
int copy(const char *oldname, const char *newname)
{
  ostringstream buf;
  buf << copy_cmd << " " << oldname << " " << newname << ends;
  return system(buf.str().c_str());
}

static struct utsname platform;
static const int ndomain=65;
static char domain[ndomain];

const char *machine()
{
  uname(&platform);
  getdomainname(domain,ndomain);
  if(strcmp(domain,"(none)") == 0) *domain=0;
  ostringstream buf;
  buf << platform.nodename << ((*domain) ? "." : "") << domain << " ("
      << platform.machine << ")" << ends;
  return strdup(buf.str().c_str());
}

const int ndate=64;
static char time_date[ndate]="";

const char *date()
{
  const time_t bintime=time(NULL);
  strftime(time_date,ndate,"%a %b %d %H:%M:%S %Z %Y",localtime(&bintime));
  return time_date;
}

const char *tempdir()
{
  return getenv("TMP_DIR");
}
