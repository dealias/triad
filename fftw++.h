#ifndef __fftwpp_h__
#define __fftwpp_h__ 1

#include <fstream>
#include <iostream>
#include <fftw3.h>
#include <complex>

#ifndef COMPLEX
typedef std::complex<double> Complex;
#endif

using std::ifstream;
using std::ofstream;
using std::cerr;
using std::endl;

inline void fftw_export_wisdom(void (*emitter)(char c, ofstream& s),
			       ofstream& s)
{
  fftw_export_wisdom((void (*) (char, void *)) emitter, (void *) &s);
}

inline int fftw_import_wisdom(int (*g)(ifstream& s), ifstream &s)
{
  return fftw_import_wisdom((int (*) (void *)) g, (void *) &s);
}

inline void PutWisdom(char c, ofstream& s) {s.put(c);}
inline int GetWisdom(ifstream& s) {return s.get();}


class fftw {
protected:
  unsigned int size;
  int sign;
  double norm;
  bool inplace;
  static const unsigned int effort=FFTW_PATIENT;
  fftw_plan plan;
  
  static bool Wise;
  static const char *WisdomName;
  static ifstream ifWisdom;
  static ofstream ofWisdom;
  
  // Shift the Fourier origin to (nx/2,0)
  void Shift(Complex *data, unsigned int nx, unsigned int nyp) {
    Complex *pstop=data+nx*nyp;
    int pinc=2*nyp;
    for(Complex *p=data+nyp; p < pstop; p += pinc) {
      //#pragma ivdep
      for(unsigned int j=0; j < nyp; j++) p[j]=-p[j];
    }
  }

  // Shift the Fourier origin to (nx/2,ny/2,0)
  void Shift(Complex *data, unsigned int nx, unsigned int ny,
	     unsigned int nz) {
    const unsigned int nzp=nz/2+1;
    const unsigned int nyzp=ny*nzp;
    const unsigned int pinc=2*nzp;
    Complex *p,*pstop;
    p=pstop=data;
    for(unsigned i=0; i < nx; i++) {
      if(i % 2) p -= nzp;
      else p += nzp;
      pstop += nyzp;
      for(; p < pstop; p += pinc) {
	//#pragma ivdep
	for(unsigned int k=0; k < nzp; k++) p[k]=-p[k];
      }
    }
  }
  
public:
  fftw(unsigned int size, int sign, unsigned int n=0) : size(size), sign(sign),
				      norm(1.0/(n ? n : size)) {}
  
  virtual void Plan(Complex *in, Complex *out)=0;
  
  void Setup(Complex *in, Complex *out) {
    if(!Wise) LoadWisdom();
    bool alloc=(in == NULL);
    if(alloc) in=new Complex[size];
    if(out == NULL) out=in;
    inplace=(out==in);
    
    Plan(in,out);
    
    if(alloc) delete [] in;
    SaveWisdom();
  }
  
  void LoadWisdom() {
    ifWisdom.open(WisdomName);
    fftw_import_wisdom(GetWisdom,ifWisdom);
    ifWisdom.close();
    Wise=true;
  }

  void SaveWisdom() {
    ofWisdom.open(WisdomName);
    fftw_export_wisdom(PutWisdom,ofWisdom);
    ofWisdom.close();
  }
  
  virtual void Execute(Complex *in, Complex *out) {
    fftw_execute_dft(plan,(fftw_complex *) in,(fftw_complex *) out);
  }
    
  void Setout(Complex *in, Complex *&out) {
    if(out == NULL) out=in;
    if(inplace ^ (out == in)) {
      cerr << "ERROR: fft constructor and call must be either both in place or both out of place" << endl; 
      exit(1);
    }
  }
  
  void fft(Complex *in, Complex *out=NULL) {
    Setout(in,out);
    Execute(in,out);
  }
    
  virtual void fftNormalized(Complex *in, Complex *out=NULL) {
    fft(in,out);
    for(unsigned int i=0; i < size; i++) out[i] *= norm;
  }
  
  void fftNormalized(Complex *in, Complex *out,
		     unsigned int nx, unsigned int m,
		     unsigned int stride, unsigned int dist) {
    if(stride == 1 && dist == nx) fftw::fftNormalized(in,out);
    else if(stride == nx && dist == 1) fftw::fftNormalized(in,out);
    else {
      fft(in,out);
      for(unsigned int k=0; k < m; k++) {
	for(unsigned int j=0; j < nx; j++) {
	  out[j*stride+k*dist] *= norm;
	}
      }
    }
  }

};

class fft : public fftw {
  unsigned int nx;
public:  
  fft(unsigned int nx, int sign, Complex *in=NULL, Complex *out=NULL) 
    : fftw(nx,sign), nx(nx) {Setup(in,out);} 
  
  void Plan(Complex *in, Complex *out) {
    plan=fftw_plan_dft_1d(nx,(fftw_complex *) in, (fftw_complex *) out,
		       sign,effort);
  }
};
  
class mfft : public fftw {
  unsigned int nx;
  unsigned int m;
  unsigned int stride;
  unsigned int dist;
public:  
  mfft(unsigned int nx, int sign, unsigned int m=1, unsigned int stride=1,
	 unsigned int dist=1, Complex *in=NULL, Complex *out=NULL) 
    : fftw((nx-1)*stride+(m-1)*dist+1,sign,nx), nx(nx), m(m),
      stride(stride), dist(dist) {Setup(in,out);} 
  
  void Plan(Complex *in, Complex *out) {
    int n[1]={nx};
    plan=fftw_plan_many_dft(1,n,m,
			    (fftw_complex *) in,NULL,stride,dist,
			    (fftw_complex *) out,NULL,stride,dist,
			    sign,effort);
  }
  
  void fftNormalized(Complex *in, Complex *out=NULL) {
    fftw::fftNormalized(in,out,nx,m,stride,dist);
  }
};
  
class fft2d : public fftw {
  unsigned int nx;
  unsigned int ny;
public:  
  fft2d(unsigned int nx, unsigned int ny, int sign, Complex *in=NULL,
	Complex *out=NULL)
    : fftw(nx*ny,sign), nx(nx), ny(ny) {Setup(in,out);} 
  
  void Plan(Complex *in, Complex *out) {
    plan=fftw_plan_dft_2d(nx,ny,(fftw_complex *) in, (fftw_complex *) out,
			  sign,effort);
  }
};

class fft3d : public fftw {
  unsigned int nx;
  unsigned int ny;
  unsigned int nz;
public:  
  fft3d(unsigned int nx, unsigned int ny, unsigned int nz,
	int sign, Complex *in=NULL, Complex *out=NULL)
    : fftw(nx*ny*nz,sign), nx(nx), ny(ny), nz(nz) {Setup(in,out);} 
  
  void Plan(Complex *in, Complex *out) {
    plan=fftw_plan_dft_3d(nx, ny, nz, (fftw_complex *) in,
			  (fftw_complex *) out, sign, effort);
  }
};

class rcfft : public fftw {
  unsigned int nx;
public:  
  rcfft(unsigned int nx, Complex *in=NULL, Complex *out=NULL) 
    : fftw(nx/2+1,-1,nx), nx(nx) {Setup(in,out);} 
  
  void Plan(Complex *in, Complex *out) {
    plan=fftw_plan_dft_r2c_1d(nx,(double *) in, (fftw_complex *) out, effort);
  }
  
  void Execute(Complex *in, Complex *out) {
    fftw_execute_dft_r2c(plan,(double *) in,(fftw_complex *) out);
  }
};
  
class crfft : public fftw {
  unsigned int nx;
public:  
  crfft(unsigned int nx, Complex *in=NULL, Complex *out=NULL) 
    : fftw(nx/2+1,1,nx), nx(nx) {Setup(in,out);} 
  
  void Plan(Complex *in, Complex *out) {
    plan=fftw_plan_dft_c2r_1d(nx,(fftw_complex *) in, (double *) out, effort);
  }
  
  void Execute(Complex *in, Complex *out) {
    fftw_execute_dft_c2r(plan,(fftw_complex *) in,(double *) out);
  }
};
  
class mrcfft : public fftw {
  unsigned int nx;
  unsigned int m;
  unsigned int stride;
  unsigned int dist;
public:  
  mrcfft(unsigned int nx, unsigned int m=1, unsigned int stride=1,
	 unsigned int dist=1, Complex *in=NULL, Complex *out=NULL) 
    : fftw(nx/2*stride+(m-1)*dist+1,-1,nx), nx(nx), m(m),
      stride(stride), dist(dist) {Setup(in,out);} 
  
  void Plan(Complex *in, Complex *out) {
    const int n[1]={nx};
    plan=fftw_plan_many_dft_r2c(1,n,m,
				(double *) in,NULL,stride,2*dist,
				(fftw_complex *) out,NULL,stride,dist,
				effort);
  }
  
  void Execute(Complex *in, Complex *out) {
    fftw_execute_dft_r2c(plan,(double *) in,(fftw_complex *) out);
  }
  
  void fftNormalized(Complex *in, Complex *out=NULL) {
    fftw::fftNormalized(in,out,(nx/2+1),m,stride,dist);
  }
};
  
class mcrfft : public fftw {
  unsigned int nx;
  unsigned int m;
  unsigned int stride;
  unsigned int dist;
public:
  mcrfft(unsigned int nx, unsigned int m=1, unsigned int stride=1,
	 unsigned int dist=1, Complex *in=NULL, Complex *out=NULL) 
    : fftw(nx/2*stride+(m-1)*dist+1,1,nx), nx(nx), m(m),
      stride(stride), dist(dist) {Setup(in,out);}
  
  void Plan(Complex *in, Complex *out) {
    const int n[1]={nx};
    plan=fftw_plan_many_dft_c2r(1,n,m,
				(fftw_complex *) in,NULL,stride,dist,
				(double *) out,NULL,stride,2*dist,
				effort);
  }
  
  void Execute(Complex *in, Complex *out) {
    fftw_execute_dft_c2r(plan,(fftw_complex *) in,(double *) out);
  }
  
  void fftNormalized(Complex *in, Complex *out=NULL) {
    fftw::fftNormalized(in,out,(nx/2+1),m,stride,dist);
  }
};
  
class rcfft2d : public fftw {
  unsigned int nx;
  unsigned int ny;
  bool shift;
public:  
  rcfft2d(unsigned int nx, unsigned int ny, Complex *in=NULL,
	  Complex *out=NULL) 
    : fftw(nx*(ny/2+1),-1,nx*ny), nx(nx), ny(ny), shift(true) {Setup(in,out);} 
  
  void Plan(Complex *in, Complex *out) {
    plan=fftw_plan_dft_r2c_2d(nx,ny,(double *) in, (fftw_complex *) out,
			      effort);
  }
  
  void fft(Complex *in, bool shift0=true, Complex *out=NULL) {
    shift=shift0;
    fftw::fft(in,out);
  }
    
  void Execute(Complex *in, Complex *out) {
    if(shift) Shift(in,nx,ny/2+1);
    fftw_execute_dft_r2c(plan,(double *) in,(fftw_complex *) out);
  }
};
  
class crfft2d : public fftw {
  unsigned int nx;
  unsigned int ny;
  bool shift;
public:  
  crfft2d(unsigned int nx, unsigned int ny, Complex *in=NULL,
	  Complex *out=NULL) 
    : fftw(nx*(ny/2+1),1,nx*ny), nx(nx), ny(ny), shift(true) {Setup(in,out);} 
  
  void Plan(Complex *in, Complex *out) {
    plan=fftw_plan_dft_c2r_2d(nx,ny,(fftw_complex *) in, (double *) out,
			      effort);
  }
  
  void fft(Complex *in, bool shift0=true, Complex *out=NULL) {
    shift=shift0;
    fftw::fft(in,out);
  }
    
  void Execute(Complex *in, Complex *out) {
    fftw_execute_dft_c2r(plan,(fftw_complex *) in,(double *) out);
    if(shift) Shift(out,nx,ny/2+1);
  }
};

class rcfft3d : public fftw {
  unsigned int nx;
  unsigned int ny;
  unsigned int nz;
  bool shift;
public:  
  rcfft3d(unsigned int nx, unsigned int ny, unsigned int nz, Complex *in=NULL,
	  Complex *out=NULL) : fftw(nx*ny*(nz/2+1),-1,nx*ny*nz),
			       nx(nx), ny(ny), nz(nz),
			       shift(true) {Setup(in,out);} 
  
  void Plan(Complex *in, Complex *out) {
    plan=fftw_plan_dft_r2c_3d(nx,ny,nz,(double *) in, (fftw_complex *) out,
			      effort);
  }
  
  void fft(Complex *in, bool shift0=true, Complex *out=NULL) {
    shift=shift0;
    fftw::fft(in,out);
  }
    
  void Execute(Complex *in, Complex *out) {
    if(shift) Shift(in,nx,ny,nz/2+1);
    fftw_execute_dft_r2c(plan,(double *) in,(fftw_complex *) out);
  }
};
  
class crfft3d : public fftw {
  unsigned int nx;
  unsigned int ny;
  unsigned int nz;
  bool shift;
public:  
  crfft3d(unsigned int nx, unsigned int ny, unsigned int nz, Complex *in=NULL,
	  Complex *out=NULL) : fftw(nx*ny*(nz/2+1),1,nx*ny*nz),
			       nx(nx), ny(ny), nz(nz),
			       shift(true) {Setup(in,out);} 
  
  void Plan(Complex *in, Complex *out) {
    plan=fftw_plan_dft_c2r_3d(nx,ny,nz,(fftw_complex *) in, (double *) out,
			      effort);
  }
  
  void fft(Complex *in, bool shift0=true, Complex *out=NULL) {
    shift=shift0;
    fftw::fft(in,out);
  }
    
  void Execute(Complex *in, Complex *out) {
    fftw_execute_dft_c2r(plan,(fftw_complex *) in,(double *) out);
    if(shift) Shift(out,nx,ny,nz/2+1);
  }
};
  
#endif
