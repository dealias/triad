#include "utils.h"

#include <ctype.h>
#include <iostream.h>
#include <errno.h>
#include <time.h>
#include <strstream.h>

#if __unix
#include <unistd.h>
#include <signal.h>
#endif

char beep='\a';

int abort_flag=1; 	// If nonzero, abort program after a fatal error.
int beep_enabled=1; // If nonzero, enable terminal beep for errors.
int msg_override=0;
void (*inform)(const char *)=NULL;

#if __unix
void wakeup(int)
{
  cout << "waking up!" << endl;
}
#endif

void msg(int severity, const char *file, int line, const char *format,...)
{
  int tty_override=0;
  char c;
  va_list vargs;
	
  cout << endl;
  if(beep_enabled) {
    if(severity == WARNING_) beep_enabled=0;
    cout << beep;
  }
	
  if(severity == OVERRIDE_ || severity == RETRY_) {
    if(msg_override) {severity=WARNING_; tty_override=0;}
    else {
#if __unix
      int errno_save=errno;
      if(isatty(STDIN_FILENO)) {
	if(severity == OVERRIDE_) tty_override=1;
	severity=WARNING_;
      } else {
	severity=ERROR_;
	errno=errno_save;
      }
#else      
      severity=ERROR_;
#endif	
    }
  }

  if(severity == ERROR_) cout << "ERROR: ";
  else cout << "WARNING: ";
  cout << flush;
	
  va_start(vargs,format);
  vform(format,vargs);
  va_end(vargs);
	
  if(*file) cout << " (\"" << file << "\":" << line << ")";
  cout << "." << endl;
	
  if(errno && errno < sys_nerr) cout << sys_errlist[errno] << "." << endl;
	
  if(tty_override) {
    cout << "Override (y/n)? ";
    cin >> c;
    if(c != 'Y' && c !='y') severity=ERROR_;
  }

  if((severity == ERROR_ && abort_flag) || severity == SLEEP_) {
    if(inform) {
      strstream buf,vbuf;
      va_start(vargs,format);
      vform(format,vargs,vbuf);
      va_end(vargs);
      vbuf << ends;
      buf << ((severity == SLEEP_ && __unix) ? "paused" : "terminated")
	  << " due to error";
      if(*file) buf << " from \"" << file << "\":" << line;
      buf << "." << newl << vbuf.str();
      if(errno && errno < sys_nerr)
	buf << " " << sys_errlist[errno] << ".";
      buf << ends;
      (*inform)(buf.str());
    }
    if(severity == SLEEP_ && __unix) {
#if __unix			
      cout << "Going to sleep..." << flush;
      signal(SIGCONT,wakeup);
      pause();
      signal(SIGCONT,SIG_DFL);
#endif			
    } else {
      exit(FATAL);
      cout << endl;
    }
  }
	
  cout << flush;
  errno=0;
}
