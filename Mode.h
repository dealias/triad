#ifndef __Mode_h__
#define __Mode_h__ 1

class Mode {
public:	
	virtual Real X() const=0;
	virtual Real Y() const=0;
	virtual Real K2() const=0;
	virtual Real K() const=0;
	virtual Real Th() const=0;
};

#endif