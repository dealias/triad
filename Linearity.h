#ifndef __Linearity_h__
#define __Linearity_h__ 1

#include "Polar.h"

#define LINEARITY(key) {new Entry<key,LinearityBase>(#key,LinearityTable);}

class LinearityBase {
public:	
	virtual char *Name() {return "None";}
	virtual Real Denominator(const Mode &v) {return v.K2();}
	
	virtual Real Growth(const Mode&) {return 0.0;}
	virtual Real Frequency(const Mode&) {return 0.0;}
	
	Real Re(const Mode& v) {return -Growth(v);}
	Real Im(const Mode& v) {return Frequency(v);}
	void Evaluate(const Mode& v, Real &nu) {nu=Re(v);}
	void Evaluate(const Mode& v, Complex &nu) {nu=Complex(Re(v),Im(v));}
};

extern LinearityBase *Linearity;
Compare_t LinearityCompare;
KeyCompare_t LinearityKeyCompare;

#endif
