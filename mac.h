#ifndef __mac_h__
#define __mac_h__ 1

char* const dirsep=":";

#include "extensions.h"

inline void setup_fpu()
{
}

inline double log1p(const double x) {return logl(x)+1.0;} // Improve
inline double expm1(const double x) {return expl(x)-1.0;} // Improve

#endif
