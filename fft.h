#ifndef __fft_h__
#define __fft_h__ 1

#include "utils.h"

void fft_init(unsigned int log2n);

void fft(Complex *data, unsigned int log2n, int isign, int bitreverse=0);
void rcfft(Complex *data, unsigned int log2n, int isign, int bitreverse=0);
void crfft(Complex *data, unsigned int log2n, int isign, int bitreverse=0);

void mfft(Complex *data, unsigned int log2n, int isign, unsigned int nk,
		  unsigned int inc1=0, unsigned int inc2=1, int bitreverse=0);
void mrcfft(Complex *data, unsigned int log2n, int isign, unsigned int nk,
			unsigned int inc1=0, unsigned int inc2=1, int bitreverse=0);
void mcrfft(Complex *data, unsigned int log2n, int isign, unsigned int nk,
			unsigned int inc1=0, unsigned int inc2=1, int bitreverse=0);

void fft2d(Complex *data, unsigned int log2nx, unsigned int log2ny, int isign);
void rcfft2d(Complex *data, unsigned int log2nx, unsigned int log2ny,
			 int isign);
void crfft2d(Complex *data, unsigned int log2nx, unsigned int log2ny,
			 int isign);

void fft4(Complex *data, unsigned int log4n, int isign);

void convolve(Complex *H, Complex *F, Complex *G, unsigned int m, unsigned
			  int log2n);
void convolve0(Complex *H, Complex *F, Complex *g, unsigned int m, unsigned
			   int log2n);


#endif
