#ifndef __types_h__
#define __types_h__ 1

#include "precision.h"
#include "Complex.h"

#ifndef COMPLEX
#define COMPLEX 0
#endif

#if(COMPLEX)
typedef Complex Var;
const Complex I(0.0,1.0);
#else
typedef Real Var;
const Real I=0.0;
#endif

typedef Real Nu;
typedef Real Mc;

#endif
