#ifndef __utils_h__
#define __utils_h__ 1

#include <stddef.h>
#include <stdarg.h>
#include <stdlib.h>
#include <fstream.h>
#include <math.h>

#if __mac
#include "arch/mac.h"
#else // !__mac

#if __i386__
#include "arch/i386.h"
#else
#include "arch/unix.h"
#endif

#endif // __mac

#include "new.h"
#include "precision.h"
#include "Complex.h"
#include "pow.h"

inline int isgn(Real x)
{
	return ((x) == 0.0 ? 0 : ((x) > 0.0 ? 1 : -1));
}

const double pi=PI;
const double twopi=2.0*pi;
const double twopi2=twopi*twopi;

const char beep='\a';

template<class T>
inline void swap(T& p, T& q) {
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

inline Real abs(const Real& x)
{
	return fabs(x);
}

inline Real abs2(const Real& x)
{
	return x*x;
}

inline Real abs2(const Complex& x)
{
	return x.real()*x.real()+x.imag()*x.imag();
}

inline Real product(const Real& x,const Real& y)
{
	return x*y;
}

inline Complex product(const Complex& x,const Complex& y)
{
	return Complex(x.real()*y.real(),x.imag()*y.imag());
}

inline Real conj(const Real& x)
{
	return x;
}

inline Real real(const Real& x)
{
	return x;
}

inline Real imag(const Real&)
{
	return 0.0;
}

inline Complex expi(const Real& phase)
{
	double cosy,siny;
	sincos(phase,&siny,&cosy);
	return Complex(cosy,siny);
}

inline Complex expim1(const Real& phase)
{
	double sinyby2=sin(0.5*phase);
	return Complex(-2.0*sinyby2*sinyby2,sin(phase));
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

inline Complex operator / (const Complex& x, const Complex& y)
{
	register double t1,t2,t3,t4;
	t3=y.re; t4=y.im; t1=t2=1.0/(t3*t3+t4*t4);
	t1*=t3; t2*=t4; t3=x.re; t4=x.im;
	return Complex(t4*t2+t3*t1,t4*t1-t3*t2);
}

inline Complex operator / (const Complex& x, Real y)
{
	return Complex(x.re/y,x.im/y);
}

inline Complex operator / (Real x, const Complex& y)
{
	register double factor;
	factor=1.0/(y.re*y.re+y.im*y.im);
	return Complex(x*y.re*factor,-x*y.im*factor);
}

inline Complex& Complex::operator /= (const Complex& y)
{
	register double t1,t2,t3,t4;
	t3=y.re; t4=y.im; t1=t2=1.0/(t3*t3+t4*t4);
	t1*=t3; t2*=t4; t3=re; t4=im;
	re*=t1;	re+=t4*t2;
	im*=t1;	im-=t3*t2;
	return *this;
}

inline Complex& Complex::operator /= (Real y)
{
	re/=y;
	im/=y;
	return *this;
}

Complex atoc(const char *s);

#ifndef __GNUC__
inline ostream& operator << (ostream& s, const Complex& y)
{
	s << "(" << y.real() << ", " << y.imag() << ")";
	return s;
}

inline istream& operator >> (istream& s, Complex& y)
{
	Real re,im;
	s >> "(" >> re >> "," >> im >> ")";
	y=Complex(re,im);
	return s;
}
#endif

void out_curve(ostream& os, Real (*f)(int), char *text, int n);
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

inline void crand_gauss(Real& w)
{
	double factor,r2,v1,v2;
	static int flag=0;
	static double save;
			  
	if (flag) {
		flag=0;
		w=save;
	} else {
		flag=1;
		do {
			v1=2.0*drand()-1.0;
			v2=2.0*drand()-1.0;
			r2=v1*v1+v2*v2;
		} while (r2 >= 1.0 || r2 == 0.0);
		factor=sqrt(-2.0*log(r2)/r2);
		w=v1*factor;
		save=v2*factor;
	}
}

inline void crand_gauss(Complex& w)
{
	double r2,v1,v2;
	do {
		v1=2.0*drand()-1.0;
		v2=2.0*drand()-1.0;
		r2=v1*v1+v2*v2;
	} while (r2 >= 1.0 || r2 == 0.0);
	w=Complex(v1,v2)*sqrt(-log(r2)/r2);
}

inline double drand_gauss()
{
	double w;
	crand_gauss(w);
	return w;
}

#endif

