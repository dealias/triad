#ifndef __Mode_h__
#define __Mode_h__ 1

#if _CFRONT
#include "options.h" // Work around Cfront template instantiation problem
#endif

class Mode {
public:	
	virtual Real X() const=0;
	virtual Real Y() const {return 0.0;}
	virtual Real K2() const;
	virtual Real K() const=0;
	virtual Real Th() const {return 0.0;}
};

#endif
