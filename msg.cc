#include "utils.h"

#include <ctype.h>
#include <iostream.h>
#include <errno.h>
#include <time.h>
#include <strstream.h>

#if __unix
#include <unistd.h>
#endif

char beep='\a';

extern int sys_nerr;
extern const char *const sys_errlist[];

int abort_flag=1; 	// If nonzero, abort program after a fatal error.
int beep_enabled=1; // If nonzero, enable terminal beep for errors.
int msg_override=0;
void (*inform)(char *)=NULL;

void msg(int severity, char *file, int line, char *format,...)
{
	int tty_override=0;
	char c;
	va_list vargs;
	
	cout << endl;
	if(beep_enabled) {
		if(severity == WARNING_) beep_enabled=0;
		cout << beep;
	}
	
	if(severity == OVERRIDE_) {
		if(msg_override) {severity=WARNING_; tty_override=0;}
		else severity=ERROR_;
#if __unix
		int errno_save=errno;
		if(isatty(STDIN_FILENO)) {severity=WARNING_; tty_override=1;}
		else errno=errno_save;
#endif	
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
			buf << ((severity == SLEEP_) ?	"paused" : "terminated")
				<< " due to error";
			if(*file) buf << " from \"" << file << "\":" << line;
			buf << "." << newl << vbuf.str() << ends;
			(*inform)(buf.str());
		}
		cout << endl;
		if(severity == SLEEP_) pause();
		else exit(FATAL);
	}
	
	cout << flush;
	errno=0;
}
