#ifndef __options_h__
#define __options_h__ 1

#define COMPLEX 1
#define MCREAL 1
#define DOUBLE_PRECISION 1

#include "utils.h"
#include "precision.h"
#include "Complex.h"

#if(COMPLEX)
typedef Complex Var;
const Complex I(0.0,1.0);
#else
typedef Real Var;
const Real I=0.0;
#endif

typedef Var Nu;

#endif
