#ifndef __Polar_h__
#define __Polar_h__ 1

class Polar {
public:	
	Real r,th;	// wavenumber components
	
	Polar(Real r0=0.0, Real th0=0.0) : r(r0), th(th0) {}

	Polar& operator = (const Polar& y) {
		r=y.r; th=y.th; return *this; 
	} 
	
	Real K() const {return r;}
	Real K2() const {return r*r;}
	Real Th() const {return th;}
	Real X() const {return r*cos(th);}
	Real Y() const {return r*sin(th);}
};

inline int operator == (const Polar& x, const Polar& y)
{
	return x.r == y.r && x.th == y.th;	
}

inline Polar operator + (const Polar& x, const Polar& y)
{
	return Polar(x.r+y.r, x.th+y.th);
}

inline Polar operator - (const Polar& x, const Polar& y)
{
	return Polar(x.r-y.r, x.th-y.th);
}

inline Polar operator * (const Polar& x, const Polar& y)
{
	return Polar(x.r*y.r, x.th*y.th);
}

inline Polar operator / (const Polar& x, const Polar& y)
{
	return Polar(x.r/y.r, x.th/y.th);
}

inline ostream& operator << (ostream& os, const Polar& y) {
	os << "(" << y.r << ", " << y.th << ")";
	return os;
}

#endif
