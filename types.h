#ifndef __types_h__
#define __types_h__ 1

#include "precision.h"
#include "Complex.h"

#define COMPLEX 0
#define MCREAL 1

#if(COMPLEX)
typedef Complex Var;
const Complex I(0.0,1.0);
#else
typedef Real Var;
const Real I=0.0;
#endif

#if(MCREAL)
typedef Real Mc;
typedef float McWeight;
#else
typedef Complex Mc;
typedef Complex McWeight;
#endif

typedef Real Nu;

#endif
