#include "utils.h"

#include <ctype.h>
#include <iostream.h>
#include <errno.h>
#include <time.h>

char beep='\a';

extern int sys_nerr;
#if  !__alpha__
extern char *sys_errlist[];
#endif

int abort_flag=1; 	// If nonzero, abort program after a fatal error.
int beep_enabled=1; // If nonzero, enable terminal beeping during errors.

#ifdef __GNUC__ 	
#define	vform(os,format,vargs) os.vform(format,vargs);
#else
#define vform(os,format,vargs) \
{os.flush(); vprintf(format,vargs); fflush(stdout);}
#endif			

#if __unix
#include <unistd.h>
#endif

#if __i386__
#include <fpu_control.h>
// Setup FPU for exceptions on overflow, zero divide and NaN.
fpu_control_t __fpu_control=
_FPU_EXTENDED | _FPU_RC_NEAREST | _FPU_MASK_DM | _FPU_MASK_UM | _FPU_MASK_PM;
#endif

void msg(int fatal, char *file, int line, char *format,...)
{
	int tty_override=0;
	char c;
	va_list vargs;

	cout << endl;
	if(beep_enabled) {beep_enabled=0; cout << beep;}
	
	if(fatal == -1) {
#if __unix
		int errno_save=errno;
		if(isatty(STDIN_FILENO)) {fatal=0; tty_override=1;}
		else errno=errno_save;
#endif	
	}

	if(fatal) cout << "ERROR: ";
	else cout << "WARNING: ";
	
	va_start(vargs,format);
	vform(cout,format,vargs);
	va_end(vargs);
	
	if(*file) cout << " (\"" << file << "\":" << line << ")" << ".";
	cout << endl;
	
	if(errno && errno < sys_nerr) cout << sys_errlist[errno] << "." << endl;
	
	if(tty_override) {
		cout << "Override (y/n)? ";
		cin >> c;
		if(c != 'Y' && c !='y') fatal=1;
	}

	if(fatal && abort_flag) {
		cout << endl;
		exit(FATAL);
	}
	cout << flush;
}
