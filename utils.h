#ifndef __utils_h__
#define __utils_h__ 1

#include "xstream.h"
#include <fstream.h>
#include <stddef.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>

#define INLINE inline

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
inline void set(T *to, const T * from, int n)
{
	memcpy(to,from,sizeof(*from)*n);
}

template<class T>
inline void swap(T& p, T& q)
{
	T temp; temp=p; p=q; q=temp;
}

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

char *upcase(const char *s, char *s2=NULL);
char *downcase(const char *s, char *s2=NULL);
char *dashify(const char *s, char *s2=NULL);
char *undashify(const char *s, char *s2=NULL);
char *convert(const char *s, char from, char to, char *s2);

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

int check_match(int match_type, char *object, char *s, int warn=1);
			   
#ifdef __GNUC__ 	
inline void vform(char *format, va_list& vargs, ostream& os=cout)
{
	os.vform(format,vargs);
	os.flush();
}
#else
// Lacking vsnprintf, formatting of arguments is available only for stdout.
inline void vform(char *format, va_list& vargs, ostream& os)
{
	os << format; 
	os.flush();
}
inline void vform(char *format, va_list& vargs)
{
	vprintf(format,vargs);
	fflush(stdout);
}
#endif

void msg(int severity, char *file, int line, char *format,...);

extern int abort_flag;
extern int beep_enabled;
extern int msg_override;
extern void (*inform)(char *);

enum ErrorCode {WARNING_,OVERRIDE_,ERROR_};

#define WARNING WARNING_,__FILE__,__LINE__
#define OVERRIDE OVERRIDE_,__FILE__,__LINE__
#define ERROR ERROR_,__FILE__,__LINE__

#define WARNING_GLOBAL WARNING_,"",0
#define OVERRIDE_GLOBAL OVERRIDE_,"",0
#define ERROR_GLOBAL ERROR_,"",0

enum ExitCode {FATAL=-1,CONTINUE,COMPLETE};
extern ExitCode exit_signal;

void mailuser(char *text);
void remove_dir(char *text);
int copy(char *oldname, char *newname);

char *machine(), *date(), *tempdir();

const int ncputime=3;
void cputime(double *cpu);

char *output_filename(char *basename, char *suffix);

inline Real divide0(Real x, Real y)
{
	return (y ? x/y : 0.0);
}

inline Real dominant(Real x)
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

Complex atoc(const char *s);
	
const int default_nperline=4;

template<class T>
inline void out_function(ostream& os, T (*f)(int), char *text, int n,
						 int nperline)
{
	int i;
	os << "# " << text << newl;
	for(i=0; i < n-1;) {
		os << (*f)(i);
		if(++i % nperline) os << "\t"; else os << " \\\n";
	}
	os << (*f)(n-1) << newl;
}

template<class T>	
inline void out_curve(ostream& os, T *f, char *text, int n, int nperline)
{
	int i;
	os << "# " << text << newl;
	for(i=0; i < n-1;) {
		os << f[i];
		if(++i % nperline) os << "\t"; else os << " \\\n";
	}
	os << f[n-1] << newl;
}

inline void out_function(ostream& os, Real (*f)(int), char *text, int n)
{
	out_function(os,f,text,n,default_nperline);
}
	
template<class T>
inline void out_curve(ostream& os, T *f, char *text, int n)
{
	out_curve(os,f,text,n,default_nperline);
}

template<class T>
inline void out_curve(ostream& os, T f, char *text)
{
	out_curve(os,&f,text,1,default_nperline);
}

void out_real(ostream& os, Complex *f, char *textre, char *, int n,
					 int nperline=default_nperline);
void out_real(ostream& os, Real *f, char *textre, char *textim, int n,
					 int nperline=default_nperline);

inline double drand()
{			  
	return ((double) rand())/RAND_MAX;
}

void crand_gauss(Real *w);
void crand_gauss(Complex *w);

inline double drand_gauss()
{
	double w;
	crand_gauss(&w);
	return w;
}

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

