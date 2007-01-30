#ifndef __unix_h__
#define __unix_h__ 1

#ifndef __unix
#define __unix 1
#endif

const char* const dirsep="/";

#include <netdb.h>

#ifndef __GNUC__
extern "C" int getdomainname(char *name, size_t len);
extern "C" int vsnprintf(char *str, size_t size, const char  *format,
			 va_list ap);
extern "C" int putenv(char *string);
#endif

#if defined __i386__ && defined __GNUC__ && defined __inline_mathop
#define log2 __log2
#define pow2 __pow2
//#define expm1 __expm1
#define log1p __log1p
#define sincos __sincos
#define sgn __sgn
//#define hypot __hypot
#else
#include "extensions.h"
extern "C" double expm1(double);
extern "C" double log1p(double);
#endif

#ifndef PI
#define PI 3.14159265358979323846264338327950288
#endif

#ifdef sun

inline istream& ws(istream& s) {s >> skipws; return s;}

extern "C" double expm1(double);
extern "C" double log1p(double);

#elifdef _AIX
#undef hz

#elifdef mips
extern "C" double expm1(double);
extern "C" double log1p(double);

#elifdef _CRAY
#ifdef _CRAYMPP
inline double log1p(const double x) {return log(x)+1.0;} // Improve
inline double expm1(const double x) {return exp(x)-1.0;} // Improve
#else
#define _CRAYMVP 1
inline double log1p(const double x) {return logl(x)+1.0;} // Improve
inline double expm1(const double x) {return expl(x)-1.0;} // Improve
#define volatile
#endif // _CRAY

#endif // sun

#endif // __unix_h__

