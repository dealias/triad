#ifndef __fftwlocal_h__
#define __fftwlocal_h__ 1

#define FFTW_1_0_COMPATIBILITY 0
#define FFTW_NOPROTO
#include <assert.h>
#include "fftw.h"
#include "rfftw.h"
#include "fstream.h"

extern "C" fftw_plan fftw_create_plan(int n, int dir, int flags);
extern "C" void fftw(fftw_plan plan, int howmany, Complex *in, int istride,
		     int idist, Complex *out, int ostride, int odist);
extern "C" void fftw_export_wisdom(void (*emitter)(char c, ofstream& s),
				   ofstream& s);
extern "C" char fftw_import_wisdom(char (*g)(ifstream& s), ifstream &s);

inline void put_wisdom(char c, ofstream& s)
{
  s.put(c);
}

inline char get_wisdom(ifstream& s)
{
  return s.get();
}

#endif
