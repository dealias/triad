#include <ctype.h>
#include <iostream.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>

#include "utils.h"

char* run=NULL;
Exit_code exit_signal=COMPLETE;

const double pi=PI;
const double twopi=2.0*pi;
const double twopi2=twopi*twopi;

char beep='\a';

extern int override;
extern int sys_nerr;
extern char *sys_errlist[];

char *upcase(const char *s, char *s2)
{
	char *p,*pstop;
	unsigned n;

	n = strlen(s);
	if(s2==NULL) s2=new char[n+1];
	pstop=s2+n;
	
	for(p=s2; p < pstop; p++) *p = toupper(*s++);
	*p='\0';

	return s2;
}

char *downcase(const char *s, char *s2)
{
	char *p,*pstop;
	unsigned n;

	n = strlen(s);
	if(s2==NULL) s2=new char[n+1];
	pstop=s2+n;
	
	for(p=s2; p < pstop; p++) *p = tolower(*s++);
	*p='\0';

	return s2;
}

char *dashify(const char *s, char *s2)
{
	return convert(s, '_', '-', s2);
}

char *undashify(const char *s, char *s2)
{
	return convert(s, '-', '_', s2);
}

extern "C" char *strdup(const char *s);

char *convert(const char *s, char from, char to, char *s2)
{
	char *p;
	if(s2==NULL) {s2=strdup(s);}
	p=s2;
	while((p=strchr(p,from))) {*p++=to;}
	return s2;
}

int RealCompare(const Real *a, const Real *b)
{
	return isgn(a-b);
}

// Generalized strncmp utility.
int strcmpn(const char *s1, const char *s2, size_t n)
{
	if(n==0) return strcmp(s1,s2);
	else return strncmp(s1,s2,n);
}

// Generalized strncasecmp utility.
int strcasecmpn(const char *s1, const char *s2, size_t n)
{
	if(n==0) return strcasecmp(s1,s2);
	else return strncasecmp(s1,s2,n);
}

// Like bsearch but also set flag match_type to:
//
//		-1 if no match is found,
//       0 if match is ambiguous,
//       1 if match is an unambiguous abbreviation,
//       2 if match is exact.
//
void *bsearch2(register const void *key,
			   register const void *base,
			   size_t nmemb,
			   register size_t size,
			   int (*compar)(const void *, const void *, const size_t),
			   int *match_type)
{
	register size_t l, u, idx, n;
	register const void *p, *p0;
	register int comparison;

	l = 0;
	u = nmemb;
	*match_type=-1;
	
	while (l < u) {
		idx = (l + u) / 2;
		p = (void *) (((const char *) base) + (idx * size));
		comparison = (*compar)((char *) key,p,0);
		if(comparison == 0) {
			*match_type=2;	// Exact match.
			return (void *) p;
		}
			
		n=strlen((char *) key);
			
		if((*compar)(key,p,n) == 0) {
			// Strings agree to n characters; an abbreviation has been found.
			*match_type=1;
				
			// Check for unambiguous match:
			// If the key also matches the previous entry to n characters,
			// it's ambiguous. 
			if(idx > 0) {
				p0 = (void *) (((const char *) base) + ((idx-1) * size));
				if((*compar)((char *) key,p0,n) == 0) *match_type=0; 
			}
				
			// If the key also matches the next entry to n characters,
			// it's ambiguous.
			if(idx < nmemb-1) {
				p0 = (void *) (((const char *) base) + ((idx+1) * size));
				if((*compar)((char *) key,p0,n) == 0) *match_type=0; 
			}
				
			if(*match_type) return (void *) p;
		}

		if (comparison < 0)
			u = idx;
		else if (comparison > 0)
			l = idx + 1;
	}

	// No match found.
	return (void *) base;
}

void check_match(int match_type, char *object, char *key)
{
	if(match_type == -1) msg(ABORT, "Unknown %s: %s",object,key);
	if(match_type == 0)  msg(ABORT, "Ambiguous %s: %s",object,key);
	return;
}

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
		if(override) {fatal=0; tty_override=0;}
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
		char *buf=new char[70+strlen(file)+strlen(format)];
		sprintf(buf,"terminated with error from \"%s\":%d; %s",
				file,line,format);
		mailuser(buf);
		exit(FATAL);
	}
	cout << flush;
}

Complex atoc(const char *s)
{
	double re,im;
	const char *ptr;
	const char sep=',';

	while(isspace(*s)) s++;
	if(*s != '(') return atof(s);
	if(!(ptr=strchr(s,sep))) return 0.0;
	re=atof(++s);
	if(!(strrchr(s,')'))) return 0.0;
	im=atof(++ptr);
	return Complex(re,im);
}

char *atos(const char *s)
{
	const int n=100;
	char *s2;
	
	s2=new char[n];
	strncpy(s2,s,n);
	return s2;
}

void crand_gauss(Real *w)
{
	double factor,r2,v1,v2;
	static int flag=0;
	static double save;
			  
	if (flag) {
		flag=0;
		*w=save;
	} else {
		flag=1;
		do {
			v1=2.0*drand()-1.0;
			v2=2.0*drand()-1.0;
			r2=v1*v1+v2*v2;
		} while (r2 >= 1.0 || r2 == 0.0);
		factor=sqrt(-2.0*log(r2)/r2);
		*w=v1*factor;
		save=v2*factor;
	}
}

void crand_gauss(Complex *w)
{
	double r2,v1,v2;
	do {
		v1=2.0*drand()-1.0;
		v2=2.0*drand()-1.0;
		r2=v1*v1+v2*v2;
	} while (r2 >= 1.0 || r2 == 0.0);
	*w=Complex(v1,v2)*sqrt(-log(r2)/r2);
}

char *output_filename(char *basename, char *suffix)
{
	char *filename;
	suffix=downcase(suffix);
	filename=new char[strlen(basename)+strlen(suffix)+1];
	strcpy(filename,basename);
	strcat(filename,suffix);
	return filename;
}

static Complex *base;
static Real out_re(int i) {return base[i].re;}
static Real out_im(int i) {return base[i].im;}

void out_real(ostream& os, Real *f, char *format, int n, int nperline)
{
	char *text=new char[strlen(format)+1];
	sprintf(text,format,"");
	out_curve(os,f,text,n,nperline);
}

void out_real(ostream& os, Complex *f, char *format, int n, int nperline)
{
	char *text=new char[strlen(format)+4];
	base=f;
	sprintf(text,format,".re");
	out_function(os,out_re,text,n,nperline);
	sprintf(text,format,".im");
	out_function(os,out_im,text,n,nperline);
}
