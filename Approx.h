#ifndef __Approx_h__
#define __Approx_h__ 1

#define APPROXIMATION(key) {new Entry<key,ApproximationBase>\
								(#key,ApproximationTable);}

class SR : public ApproximationBase {
public:
	char *Name() {return "Spectral Reduction";}
	void SetSrcRoutines(Source_t **LinearSrc, Source_t **NonlinearSrc,
						Source_t **ConstantSrc);
};

class None : public ApproximationBase {
public:
	char *Name() {return "None";}
	void SetSrcRoutines(Source_t **LinearSrc, Source_t **NonlinearSrc,
						Source_t **ConstantSrc);
};

#endif

