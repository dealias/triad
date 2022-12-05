#ifndef __newton_h__
#define __newton_h__ 1

#include "precision.h"

bool newton(Real &x, Real (*f)(Real x), Real (*dfdx)(Real x),
	    bool verbose=false, size_t MaxIterations=100);

bool newton(Real &x1, Real x2, Real (*f)(Real x), Real (*dfdx)(Real x),
	    bool verbose=false, size_t MaxIterations=100);

#endif
