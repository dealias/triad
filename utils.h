#ifndef __utils_h__
#define __utils_h__ 1

#include <stddef.h>
#include <stdarg.h>
#include <stdlib.h>
#include <fstream.h>
#include <math.h>

#if _MAC
#include ":arch:mac.h"
#else // !_MAC

#if __i386__
#include "arch/i386.h"
#else
#include "arch/unix.h"
#endif

#endif // _MAC

#include "new.h"
#include "precision.h"
#include "Complex.h"
#include "pow.h"

extern const double pi;
extern const double twopi;
extern const double twopi2;

extern char beep;

inline int isgn(const Real x)
{
	return ((x) == 0.0 ? 0 : ((x) > 0.0 ? 1 : -1));
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
inline void swap(T& p, T& q)
{
	T temp; temp=p; p=q; q=temp;
}

template<class T>
inline void sort2(T& p, T& q, int& sign)
{
	if(p > q) {swap(p,q); sign *= -1;}
}

template<class T>
inline T min(T a, T b)
{
	return (a < b) ? a : b;
}

template<class T>
inline T max(T a, T b)
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

void check_match(int match_type, char *object, char *s);
			   
#define OVERRIDE -1,__FILE__,__LINE__
#define WARNING 0,__FILE__,__LINE__
#define ERROR 1,__FILE__,__LINE__
#define ABORT 1,"",0

void msg(int fatal, char *file, int line, char *format,...);

#define INVALID_CALL msg(ERROR, "Invalid call")
#define INVALID_ARG msg(ERROR, "Invalid argument")

void mailuser(char *text);

const int ncputime=3;
void cputime(double *cpu);

char *output_filename(char *basename, char *suffix);

#if !_AIX
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

inline Real product(Real x, Real y)
{
	return x*y;
}

inline Complex product(const Complex& x, const Complex& y)
{
	return Complex(x.re*y.re,x.im*y.im);
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

inline Complex expim1(Real phase)
{
	double sinyby2=sin(0.5*phase);
	return Complex(-2.0*sinyby2*sinyby2,sin(phase));
}

Complex atoc(const char *s);

void out_curve(ostream& os, Real (*f)(int), char *text, int n);
void out_curve(ostream& os, Complex (*f)(int), char *text, int n);
void out_curve(ostream& os, int *f, char *text, int n);
void out_curve(ostream& os, Real *f, char *text, int n);
void out_curve(ostream& os, Complex *f, char *text, int n);

inline void out_curve(ostream& os, int f, char *text)
{
	out_curve(os,&f,text,1);
}

inline void out_curve(ostream& os, Real f, char *text)
{
	out_curve(os,&f,text,1);
}

inline void out_curve(ostream& os, Complex f, char *text)
{
	out_curve(os,&f,text,1);
}

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

#endif

