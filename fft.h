#ifndef __fft_h__
#define __fft_h__ 1

#include "utils.h"

void fft_br(Complex *data, unsigned int log2n);
void fft_brinv(Complex *data, unsigned int log2n);
void fft(Complex *data, unsigned int log2n, int isign, unsigned int nk);
void fft(Complex *data, unsigned int log2n, int isign);
void fft4(Complex *data, unsigned int log4n, int isign);

void rfft_br(Complex *data, unsigned int log2n);
void rfft_brinv(Complex *data, unsigned int log2n);
	
#endif

