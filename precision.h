#ifndef __precision_h__
#define __precision_h__ 1

#include <float.h>

#ifndef DOUBLE_PRECISION
#define DOUBLE_PRECISION 1
#endif

#if(DOUBLE_PRECISION)
typedef double Real;
#define REAL_MIN DBL_MIN
#define REAL_MAX DBL_MAX
#define REAL_EPSILON DBL_EPSILON
#define REAL_DIG DBL_DIG
#else
typedef float Real;
#define REAL_MIN FLT_MIN
#define REAL_MAX FLT_MAX
#define REAL_EPSILON FLT_EPSILON
#define REAL_DIG FLT_DIG
#endif

#endif
