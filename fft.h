#ifndef __fft_h__
#define __fft_h__ 1

void fft_init(unsigned int log2n);

void fft(Complex *data, unsigned int log2n, int isign, Real scale=1.0,
	 int bitreverse=0); 
void rcfft(Complex *data, unsigned int log2n, int isign, Real scale=1.0, 
	   int bitreverse=0);
void crfft(Complex *data, unsigned int log2n, int isign, Real scale=1.0,
	   int bitreverse=0);

void signscalefft(Complex *data, unsigned int n, int isign, unsigned int nk,
		  unsigned int inc1, unsigned int inc2, Real scale);
	
void mfft(Complex *data, unsigned int log2n, int isign, unsigned int nk,
	  unsigned int inc1=0, unsigned int inc2=1, Real scale=1.0,
	  int bitreverse=0);
void mrcfft(Complex *data, unsigned int log2n, int isign, unsigned int nk,
	    unsigned int inc1=0, unsigned int inc2=1, Real scale=1.0,
	    int bitreverse=0);
void mcrfft(Complex *data, unsigned int log2n, int isign, unsigned int nk,
	    unsigned int inc1=0, unsigned int inc2=1, Real scale=1.0,
	    int bitreverse=0);
void mrcfft0(Complex *data, unsigned int log2n, int isign, unsigned int nk,
	     unsigned int inc1=0, unsigned int inc2=1, Real scale=1.0,
	     int bitreverse=0);
void mcrfft0(Complex *data, unsigned int log2n, int isign, unsigned int nk,
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
