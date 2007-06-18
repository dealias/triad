#ifndef __options_h__
#define __options_h__ 1

#define DOUBLE_PRECISION 1
#define COMPLEX 0

#include "utils.h"

#if(COMPLEX)
typedef Complex Var;
inline Complex rand_gauss() {return crand_gauss();}
#else
typedef Real Var;
inline Real rand_gauss() {return drand_gauss();}
#endif

#if(NUCOMPLEX)
typedef Var Nu;
#else
typedef Real Nu;
#endif

#endif
