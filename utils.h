#ifndef __utils_h__
#define __utils_h__ 1

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstddef>
#include <cstdarg>
#include <cstdlib>
#include <cmath>
#include <memory>
#include <cerrno>
#include <cstring>

#include "xstream.h"

#ifdef __APPLE__
#include "mac.h"
#else // !__APPLE__

#include "unix.h"

#endif // __APPLE__

#include <new>
void *operator new(size_t size, int);
void *operator new(size_t size, void *ptr, int new_len);
size_t memory();

#include "precision.h"
#include "Complex.h"

#include "pow.h"

#define __ArrayExtensions
#define __ExternalArrayExit
#define __ExternalDynVectorExit
#include "Array.h"

using std::istream;
using std::ostream;
using std::ifstream;
using std::ofstream;
using xdr::ixstream;
using xdr::oxstream;
using std::cin;
using std::cout;
using std::cerr;
using std::endl;
using std::ends;
using std::flush;
using std::setw;
using std::setfill;
using std::setprecision;

namespace triad {
extern const double pi;
extern const double twopi;
extern const double twopi2;
}

extern char beep;

inline ostream& newl(ostream& s) {s << '\n'; return s;}
inline oxstream& newl(oxstream& s) {return s;}

inline Real sgn1(Real x)
{
  return x >= 0.0 ? 1.0 : -1.0;
}

inline int isgn(Real x)
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
inline void set(T *to, const T *from, size_t n)
{
  memcpy((void *) to,(void *) from,sizeof(T)*n);
}

template<class T>
inline void set(Array::array1<T> to, const Array::array1<T> from, size_t n)
{
  memcpy((void *) to(),(void *) from(),sizeof(T)*n);
}

template<class T>
inline void swap(T& p, T& q)
{
  T temp=p; p=q; q=temp;
}

template<class T>
inline void sort2(T& p, T& q, int& sign)
{
  if(p > q) {swap(p,q); sign *= -1;}
}

template<class T, class S>
inline void sort2x2(T& p, T& q, S& a, S& b, int& sign)
{
  if(p > q || (p == q && a > b)) {swap(p,q); swap(a,b); sign *= -1;}
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

inline size_t min(size_t a, size_t b)
{
  return (a < b) ? a : b;
}

inline size_t max(size_t a, size_t b)
{
  return (a > b) ? a : b;
}

char *upcase(const char *s, char *s2=NULL);
char *downcase(const char *s, char *s2=NULL);
char *dashify(const char *s, char *s2=NULL);
char *undashify(const char *s, char *s2=NULL);
char *convert(const char *s, char from, char to, char *s2);

int RealCompare(const void *a, const void *b);

#if 0
extern "C" char *strdup(const char *);
extern "C" int strcasecmp (const char *s1, const char *s2);
extern "C" int strncasecmp (const char *s1, const char *s2, size_t n);
#endif

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

inline void vform(const char *format, va_list& vargs, ostream& os=cout)
{
  const int maxsize=256;
  char buf[maxsize];
  vsnprintf(buf,maxsize,format,vargs);
  os << buf;
  os.flush();
}

void msg(int severity, const char *file, int line, const char *format,...);

extern int abort_flag;
extern int beep_enabled;
extern int msg_override;
extern void (*inform)(const char *);

enum ErrorCode {WARNING_,OVERRIDE_,RETRY_,SLEEP_,ERROR_};

#define WARNING WARNING_,__FILE__,__LINE__
#define OVERRIDE OVERRIDE_,__FILE__,__LINE__
#define RETRY RETRY_,__FILE__,__LINE__
#define SLEEP SLEEP_,__FILE__,__LINE__
#define ERROR ERROR_,__FILE__,__LINE__

#define WARNING_GLOBAL WARNING_,"",0
#define OVERRIDE_GLOBAL OVERRIDE_,"",0
#define RETRY_GLOBAL RETRY_,"",0
#define SLEEP_GLOBAL SLEEP_,"",0
#define ERROR_GLOBAL ERROR_,"",0

inline void Array::ArrayExit(const char *x) {errno=0; msg(ERROR_GLOBAL,x);}
inline void DynVectorExit(const char *x) {errno=0; msg(ERROR_GLOBAL,x);}

enum ExitCode {FATAL=-1,CONTINUE,COMPLETE};
extern ExitCode exit_signal;

void mailuser(const char *text);
void remove_dir(const char *text);
int copy(const char *oldname, const char *newname);

const char *machine(), *date(), *tempdir();


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


#if !defined(_AIX) && !defined(__GNUC__)
inline Real abs(Real x)
{
  return fabs(x);
}
#endif

inline Real abs2(Real x)
{
  return x*x;
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

/*
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
*/

inline void LeastSquaresFit(unsigned n, Real *x, Real *y, Real &m, Real &b)
{
  Real sumx=0.0, sumy=0.0, sumxx=0.0, sumxy=0.0;
  for(unsigned i=0; i < n; i++) {
    sumx += x[i];
    sumy += y[i];
    sumxx += x[i]*x[i];
    sumxy += x[i]*y[i];
  }
  m=(n*sumxy-sumx*sumy)/(n*sumxx-sumx*sumx);
  b=(sumy-m*sumx)/n;
}

char *atos(const char *s);
Complex atoc(const char *s);
size_t atou(const char *s);

const int default_nperline=4;

const char *const outerror="Cannot write %s to output stream";

template<class T>
inline void out_curve(ofstream& os, T *f, const char *text, size_t n=1,
		      int nperline=default_nperline)
{
  size_t i;
  os << "# " << text << newl;
  if(n == 0) return;
  for(i=0; i < n-1;) {
    os << f[i];
    if(++i % nperline) os << "\t"; else os << " \\\n";
  }
  os << f[n-1] << newl;
  if(!os) msg(WARNING,outerror,text);
}

template<class T>
inline void out_curve(oxstream& os, T *f, const char *text, size_t n=1,
		      int=default_nperline)
{
  os << n;
  for(size_t i=0; i < n; i++) os << f[i];
  if(!os) msg(WARNING,outerror,text);
}

template<class T>
inline void out_curve(ofstream& os, T (*f)(size_t), const char *text,
		      size_t n, int nperline=default_nperline)
{
  size_t i;
  os << "# " << text << newl;
  if(n == 0) return;
  for(i=0; i < n-1;) {
    os << (*f)(i);
    if(++i % nperline) os << "\t"; else os << " \\\n";
  }
  os << (*f)(n-1) << newl;
  if(!os) msg(WARNING,outerror,text);
}

template<class T>
inline void out_curve(oxstream& os, T (*f)(size_t), const char *text,
		      size_t n, int=default_nperline)
{
  os << n;
  for(size_t i=0; i < n; i++) os << (*f)(i);
  if(!os) msg(WARNING,outerror,text);
}

template<class S, class T>
inline void out_curve(S& os, T f, const char *text,
		      int nperline=default_nperline)
{
  out_curve(os,&f,text,1,nperline);
}

template<class S, class T>
inline void out_curve(S& os, Array::array1<T> f, const char *text,
		      int nperline=default_nperline)
{
  out_curve(os,f(),text,1,nperline);
}

inline double drand()
{
  static const double factor=1.0/RAND_MAX;
  return random()*factor;
}

inline int rand_sign()
{
  return random() % 2 ? 1 : -1;
}

inline Complex rand_unityroot(int n)
{
  return expi(triad::twopi*(random() % n)/n);
}

inline Complex rand_phase()
{
  return expi(triad::twopi*drand());
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
