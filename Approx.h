#ifndef __Approx_h__
#define __Approx_h__ 1

#define APPROXIMATION(key) {new Entry<key,ApproximationBase>\
								(#key,ApproximationTable);}

class SR : public ApproximationBase {
public:
	char *Name() {return "Spectral Reduction";}
};

class None : public ApproximationBase {
public:
	char *Name() {return "None";}
};

#endif

