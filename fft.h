#ifndef __fft_h__
#define __fft_h__ 1

#include "utils.h"

void fft_init(unsigned int log2n);

void fft(Complex *data, unsigned int log2n, int isign, int bitreverse=0);
void rfft(Complex *data, unsigned int log2n, int bitreverse=0);
void rfft_inv(Complex *data, unsigned int log2n, int bitreverse=0);

void mfft(Complex *data, unsigned int log2n, int isign, unsigned int nk);
void fft4(Complex *data, unsigned int log4n, int isign);

#endif

