#ifndef __Bin_h__
#define __Bin_h__ 1

#include "utils.h"

template<class T>
class Bin {
public:
	Bin() {}
	T min,max,cen;
	T Delta() {return (max-min);}
	Real Area();
	Real K(),Th();
	Real Kx(),Ky();
};

template<class T>
inline ostream& operator << (ostream& os, const Bin<T>& y) {
	os << "[" << y.min << "\t" << y.cen << "\t" << y.max << "]" << endl;
	return os;
}

#endif
