#ifndef __fft_h__
#define __fft_h__ 1

#define FFTW_1_0_COMPATIBILITY 0
#define FFTW_NOPROTO
#include "fftw.h"
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

void fft_init(unsigned int log2n);

void fft(Complex *data, unsigned int log2n, int isign, Real scale=1.0,
		 int bitreverse=0); 
void rcfft(Complex *data, unsigned int log2n, int isign, Real scale=1.0, 
		   int bitreverse=0);
void crfft(Complex *data, unsigned int log2n, int isign, Real scale=1.0,
		   int bitreverse=0);

void mfft(Complex *data, unsigned int log2n, int isign, unsigned int nk,
		  unsigned int inc1=0, unsigned int inc2=1, Real scale=1.0,
		  int bitreverse=0);
void mrcfft(Complex *data, unsigned int log2n, int isign, unsigned int nk,
			unsigned int inc1=0, unsigned int inc2=1, Real scale=1.0,
			int bitreverse=0);
void mcrfft(Complex *data, unsigned int log2n, int isign, unsigned int nk,
			unsigned int inc1=0, unsigned int inc2=1, Real scale=1.0,
			int bitreverse=0);

void fft2d(Complex *data, unsigned int log2nx, unsigned int log2ny,
		   int isign, Real scale=1.0, int bitreverse=0);
void rcfft2d(Complex *data, unsigned int log2nx, unsigned int log2ny,
			 int isign, Real scale=1.0, int bitreverse=0);
void crfft2d(Complex *data, unsigned int log2nx, unsigned int log2ny,
			 int isign, Real scale=1.0, int bitreverse=0);
void crfft2d_sym(Complex *data, unsigned int log2nx, unsigned int log2ny,
				 int isign, Real scale=1.0, int bitreverse=0);
void rcfft2dT(Complex *data, unsigned int log2nx, unsigned int log2ny,
			  int isign, Real scale=1.0, int bitreverse=0);
void crfft2dT(Complex *data, unsigned int log2nx, unsigned int log2ny,
			  int isign, Real scale=1.0, int bitreverse=0);
void crfft2dT_sym(Complex *data, unsigned int log2nx, unsigned int log2ny,
				  int isign, Real scale=1.0, int bitreverse=0);

void fft3d(Complex *data, unsigned int log2nx, unsigned int log2ny,
		   unsigned int log2nz, int isign, Real scale=1.0, int bitreverse=0);
void rcfft3d(Complex *data, unsigned int log2nx, unsigned int log2ny,
			 unsigned int log2nz, int isign, Real scale=1.0, int bitreverse=0);
void crfft3d(Complex *data, unsigned int log2nx, unsigned int log2ny,
			 unsigned int log2nz, int isign, Real scale=1.0, int bitreverse=0);
void crfft3d_sym(Complex *data, unsigned int log2nx, unsigned int log2ny,
				 unsigned int log2nz, int isign, Real scale=1.0, 
				 int bitreverse=0);
	
void convolve(Complex *H, Complex *F, Complex *G, unsigned int m, unsigned
			  int log2n);
void convolve0(Complex *H, Complex *F, Complex *g, unsigned int m, unsigned
			   int log2n);

#endif
