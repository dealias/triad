#ifndef __utils_h__
#define __utils_h__ 1

#include <stddef.h>
#include <stdarg.h>
#include <stdlib.h>
#include <fstream.h>
#include <math.h>
#include <stdio.h>

#if _MAC
#include ":arch:mac.h"
#else // !_MAC

#if __i386__
#include "arch/i386.h"
#else
#include "arch/unix.h"
#endif

#endif // _MAC

#include <new.h>
void *operator new(size_t size, void *ptr, int new_len);

#include "precision.h"
#include "Complex.h"
#include "pow.h"

extern const double pi;
extern const double twopi;
extern const double twopi2;

extern char beep;

const char newl='\n';

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
#define CONDITIONAL -1:1,__FILE__,__LINE__
#define ABORT 1,"",0

void msg(int fatal, char *file, int line, char *format,...);

#define INVALID_CALL msg(ERROR, "Invalid call")
#define INVALID_ARG msg(ERROR, "Invalid argument")

enum Exit_code {FATAL=-1,CONTINUE,COMPLETE};
extern Exit_code exit_signal;

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

inline Complex expim1(Real phase)
{
	double sinyby2=sin(0.5*phase);
	return Complex(-2.0*sinyby2*sinyby2,sin(phase));
}

Complex atoc(const char *s);
	
const int default_nperline=4;

#if _CRAY
// Cfront can't seem to handle a template here.
inline void out_function(ostream& os, Real (*f)(int), char *text, int n,
						 int nperline)
#else
template<class T>
void out_function(ostream& os, T (*f)(int), char *text, int n)
#endif	
{
	int i;
	os << "# " << text << newl;
	for(i=0; i < n-1;) {
		os << (*f)(i);
		if(++i % nperline) os << "\t"; else os << " \\\n";
	}
	os << (*f)(n-1) << newl;
}

#if _CRAY
// Cfront can't seem to handle a template here.
inline void out_curve(ostream& os, Real *f, char *text, int n, int nperline)
#else
template<class T>	
void out_curve(ostream& os, T *f, char *text, int n, int nperline)
#endif
{
	int i;
	os << "# " << text << newl;
	for(i=0; i < n-1;) {
		os << f[i];
		if(++i % nperline) os << "\t"; else os << " \\\n";
	}
	os << f[n-1] << newl;
}

#if 0
inline void out_function(ostream& os, Real (*f)(int), char *text, int n)
{
	out_function(os,f,text,n,default_nperline);
}
#endif
	
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

void convolve(Complex *H, Complex *F, Complex *G, unsigned int m, unsigned
			  int log2n);
void convolve0(Complex *H, Complex *F, Complex *g, unsigned int m, unsigned
			   int log2n);

#endif

