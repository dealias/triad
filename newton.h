#ifndef __newton_h__
#define __newton_h__ 1

#include "precision.h"

bool newton(Real &x1, Real x2, Real (*f)(Real x), Real (*dfdx)(Real x),
	    Real epsilon, bool verbose, unsigned int MaxIterations);

bool newton(Real &x1, Real x2, Real (*f)(Real x), Real (*dfdx)(Real x),
	    Real epsilon, bool verbose, unsigned int MaxIterations);

#endif
