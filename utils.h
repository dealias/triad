#ifndef __utils_h__
#define __utils_h__ 1

#include "xstream.h"
#include <fstream.h>
#include <stddef.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>

#if _MAC
#include ":arch:mac.h"
#else // !_MAC

#include "arch/unix.h"

#endif // _MAC

#include <new.h>
void *operator new(size_t size, int);
void *operator new(size_t size, void *ptr, int new_len);
size_t memory();

#include "precision.h"
#include "Complex.h"
#include "pow.h"

#define __ARRAY_EXIT(x) msg(ERROR,x)
	
extern const double pi;
extern const double twopi;
extern const double twopi2;

extern char beep;

inline ostream& newl(ostream& s) {s << '\n'; return s;}
inline oxstream& newl(oxstream& s) {return s;}

inline Real sgn1(const Real x)
{
	return x >= 0.0 ? 1.0 : -1.0;
}

inline int isgn(const Real x)
{
	return x == 0.0 ? 0 : (x > 0.0 ? 1 : -1);
}

// x mod n in [0,n) where n > 0 (implementation-independent definition)
inline int mod(int x, int n)
{
	return x < 0 ? (n-(-x % n)) % n : x % n;
}

inline Real dmod(Real x, Real n) {
	return x < 0 ? x-floor(x/n)*n : fmod(x,n);
}

template<class T> 
inline void set(T *to, const T * from, unsigned int n)
{
	memcpy(to,from,sizeof(*from)*n);
}

#ifndef __SGI_STL_INTERNAL_ALGOBASE_H
template<class T>
inline void swap(T& p, T& q)
{
	T temp=p; p=q; q=temp;
}
#endif

template<class T>
inline void sort2(T& p, T& q, int& sign)
{
	if(p > q) {swap(p,q); sign *= -1;}
}

inline Real min(Real a, Real b)
{
	return (a < b) ? a : b;
}

inline Real max(Real a, Real b)
{
	return (a > b) ? a : b;
}

inline int min(int a, int b)
{
	return (a < b) ? a : b;
}

inline int max(int a, int b)
{
	return (a > b) ? a : b;
}

inline unsigned int min(unsigned int a, unsigned int b)
{
	return (a < b) ? a : b;
}

inline unsigned int max(unsigned int a, unsigned int b)
{
	return (a > b) ? a : b;
}

char *upcase(const char *s, char *s2=NULL);
char *downcase(const char *s, char *s2=NULL);
char *dashify(const char *s, char *s2=NULL);
char *undashify(const char *s, char *s2=NULL);
char *convert(const char *s, char from, char to, char *s2);

const char null[1]="";
extern "C" char *strdup(const char *);

extern "C" int strcasecmp (const char *s1, const char *s2);
extern "C" int strncasecmp (const char *s1, const char *s2, size_t n);

int strcmpn(const char *s1, const char *s2, size_t n);
int strcasecmpn(const char *s1, const char *s2, size_t n);

typedef int Compare_t(const void *, const void *);
typedef int KeyCompare_t(const void *, const void *, const size_t);

void *bsearch2(register const void *key,
			   register const void *base,
			   size_t nmemb,
			   register size_t size,
			   int (*compar)(const void *, const void *, const size_t),
			   int *match_type);

int check_match(int match_type, const char *object, const char *s, int warn=1);
			   
#ifdef __GNUC__ 	
inline void vform(const char *format, va_list& vargs, ostream& os=cout)
{
	os.vform(format,vargs);
	os.flush();
}
#else
// Lacking vsnprintf, formatting of arguments is available only for stdout.
inline void vform(const char *format, va_list& vargs, ostream& os)
{
	os << format; 
	os.flush();
}
inline void vform(const char *format, va_list& vargs)
{
	vprintf(format,vargs);
	fflush(stdout);
}
#endif

void msg(int severity, const char *file, int line, const char *format,...);

extern int abort_flag;
extern int beep_enabled;
extern int msg_override;
extern void (*inform)(const char *);

enum ErrorCode {WARNING_,OVERRIDE_,SLEEP_,ERROR_};

#define WARNING WARNING_,__FILE__,__LINE__
#define OVERRIDE OVERRIDE_,__FILE__,__LINE__
#define SLEEP SLEEP_,__FILE__,__LINE__
#define ERROR ERROR_,__FILE__,__LINE__

#define WARNING_GLOBAL WARNING_,"",0
#define OVERRIDE_GLOBAL OVERRIDE_,"",0
#define SLEEP_GLOBAL SLEEP_,"",0
#define ERROR_GLOBAL ERROR_,"",0

enum ExitCode {FATAL=-1,CONTINUE,COMPLETE};
extern ExitCode exit_signal;

void mailuser(const char *text);
void remove_dir(const char *text);
int copy(const char *oldname, const char *newname);

char *machine(), *date(), *tempdir();

const int ncputime=3;
void cputime(double *cpu);

char *output_filename(char *basename, char *suffix);

inline Real divide0(Real x, Real y)
{
	return (y ? x/y : 0.0);
}

inline Real max(Real x)
{
	return x;
	
}

inline Real min(Real x)
{
	return x;
}

#if !_AIX || __GNUC__
inline Real abs(Real x)
{
	return fabs(x);
}
#endif

inline Real abs2(Real x)
{
	return x*x;
}

inline Real abs2(const Complex& x)
{
	return x.re*x.re+x.im*x.im;
}

inline Real norm2(Real x)
{
	return abs2(x);
}

inline Real norm2(const Complex& x)
{
	return abs2(x);
}

inline Real norm(Real x)
{
	return abs(x);
}

inline Real product(Real x, Real y)
{
	return x*y;
}

inline Complex product(const Complex& x, const Complex& y)
{
	return Complex(x.re*y.re,x.im*y.im);
}

inline Real realproduct(Real x, Real y)
{
	return x*y;
}

inline Real realproduct(const Complex& x, const Complex& y)
{
	return x.re*y.re+x.im*y.im;
}

inline void conjugate(Real& x, Real y)
{
	x=y;
}

inline void conjugate(Complex& x, const Complex& y)
{
	x.re=y.re;
	x.im=-y.im;
}

inline Real conj(Real x)
{
	return x;
}

inline Real real(Real x)
{
	return x;
}

inline Real imag(Real)
{
	return 0.0;
}

inline Complex exp(const Complex& z)
{
	double cosy,siny;
	Complex w;
	
	sincos(z.im,&siny,&cosy);
	w.im=exp(z.re);
	w.re=w.im*cosy;
	w.im*=siny;
	return w;
}

inline Complex expm1(const Complex& z)
{
	double cosy,siny,sinyby2;
	Complex w;
	
	sincos(z.im,&siny,&cosy);
	w.im=exp(z.re);
	w.re=w.im*cosy-1.0;
	w.im*=siny;
	if(0.0 < cosy && w.re < 1.0) {
		sinyby2=sin(0.5*z.im);
		w.re=expm1(z.re+log1p(-2.0*sinyby2*sinyby2));
	}
	return w;
}

inline Complex expi(Real phase)
{
	double cosy,siny;
	sincos(phase,&siny,&cosy);
	return Complex(cosy,siny);
}

char *atos(const char *s);
Complex atoc(const char *s);
	
const int default_nperline=4;

static const char *outerror="Cannot write %s to output stream";

template<class T>	
inline void out_curve(oxstream& os, T *f, const char *text, unsigned int n,
					  int)
{
	os << n;
	for(unsigned int i=0; i < n; i++) os << f[i];
	if(!os) msg(WARNING,outerror,text);
}

template<class T>
inline void out_function(oxstream& os, T (*f)(unsigned int), const char *text, 
						 unsigned int n, int) 
{
	os << n;
	for(unsigned int i=0; i < n; i++) os << (*f)(i);
	if(!os) msg(WARNING,outerror,text);
}

template<class T>
inline void out_function(ostream& os, T (*f)(unsigned int), const char *text,
						 unsigned int n, int nperline)
{
	unsigned int i;
	os << "# " << text << newl;
	for(i=0; i < n-1;) {
		os << (*f)(i);
		if(++i % nperline) os << "\t"; else os << " \\\n";
	}
	os << (*f)(n-1) << newl;
	if(!os) msg(WARNING,outerror,text);
}

template<class T>	
inline void out_curve(ostream& os, T *f, const char *text, unsigned int n,
					  int nperline)
{
	unsigned int i;
	os << "# " << text << newl;
	for(i=0; i < n-1;) {
		os << f[i];
		if(++i % nperline) os << "\t"; else os << " \\\n";
	}
	os << f[n-1] << newl;
	if(!os) msg(WARNING,outerror,text);
}

template<class S>
inline void out_function(S& os, Real (*f)(unsigned int), const char *text,
						 unsigned int n)
{
	out_function(os,f,text,n,default_nperline);
}
	
template<class S, class T>
inline void out_curve(S& os, T *f, const char *text, unsigned int n)
{
	out_curve(os,f,text,n,default_nperline);
}

template<class S, class T>
inline void out_curve(S& os, T f, const char *text)
{
	out_curve(os,&f,text,(unsigned int) 1,default_nperline);
}

template<class S>
void out_real(S& os, Real *f, const char *textre, const char *,
			  unsigned int n, int nperline) 
{
	out_curve(os,f,textre,n,nperline);
}

template<class S>
void out_real(S& os, Real *f, const char *textre, const char *, unsigned int n)
{
	out_curve(os,f,textre,n,default_nperline);
}

extern Complex *out_base;
Real out_re(unsigned int i);
Real out_im(unsigned int i);

template<class S>
void out_real(S& os, Complex *f, const char *textre, const char *textim,
			  unsigned int n, int nperline) 
{
	out_base=f;
	out_function(os,out_re,textre,n,nperline);
	out_function(os,out_im,textim,n,nperline);
}

template<class S>
void out_real(S& os, Complex *f, const char *textre, const char *textim,
			  unsigned int n)
{
	out_base=f;
	out_function(os,out_re,textre,n,default_nperline);
	out_function(os,out_im,textim,n,default_nperline);
}

inline double drand()
{			  
	return ((double) rand())/RAND_MAX;
}

double drand_gauss();
Complex crand_gauss();

inline ixstream& operator >> (ixstream& s, Complex& y)
{
	s >> y.re >> y.im;
	return s;
}

inline oxstream& operator << (oxstream& s, const Complex& y)
{
	s << y.re << y.im;
	return s;
}

#endif

