#ifndef __Linearity_h__
#define __Linearity_h__ 1

#include "kernel.h"
#include "Polar.h"

#define LINEARITY(key) {new Entry<key,LinearityBase>(#key,LinearityTable);}

class LinearityBase {
public:	
	virtual char *Name()=0;
	virtual Real Denominator(Real k2) {return k2;}
	
	virtual Real Growth(const Polar& v) {return 0.0;}
	virtual Real Frequency(const Polar& v) {return 0.0;}
	
	Real Re(const Polar& v) {return -Growth(v);}
	Real Im(const Polar& v) {return Frequency(v);}
	void Evaluate(const Polar& v, Real &nu) {nu=Re();}
	void Evaluate(const Polar& v, Complex &nu) {nu=Re()+I*Im();}
};

extern LinearityBase *Linearity;
Compare_t LinearityCompare;
KeyCompare_t LinearityKeyCompare;

#endif
