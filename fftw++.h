#ifndef __fftwpp_h__
#define __fftwpp_h__ 1

#include <fstream>
#include <fftw3.h>

inline void fftw_export_wisdom(void (*emitter)(char c, ofstream& s),
			       ofstream& s)
{
  fftw_export_wisdom((void (*) (char, void *)) emitter, (void *) &s);
}

inline int fftw_import_wisdom(int (*g)(ifstream& s), ifstream &s)
{
  return fftw_import_wisdom((int (*) (void *)) g, (void *) &s);
}

inline void put_wisdom(char c, ofstream& s)
{
  s.put(c);
}

inline int get_wisdom(ifstream& s)
{
  return s.get();
}


class fftw {
protected:
  unsigned int size;
  int sign;
  double norm;
  bool inplace;
  static const unsigned int effort=FFTW_PATIENT;
  fftw_plan plan;
  
  static bool wise;
  static const char *wisdom_name;
  static ifstream ifwisdom;
  static ofstream ofwisdom;
  
public:
  fftw(unsigned int size, int sign, unsigned int n=0) : size(size), sign(sign),
				      norm(1.0/(n ? n : size)) {}
  
  virtual void Plan(Complex *in, Complex *out)=0;
  
  void Setup(Complex *in, Complex *out) {
    if(!wise) LoadWisdom();
    bool alloc=(in == NULL);
    if(alloc) in=new Complex[size];
    if(out == NULL) out=in;
    inplace=(out==in);
    
    Plan(in,out);
    
    if(alloc) delete [] in;
    SaveWisdom();
  }
  
  void LoadWisdom() {
    ifwisdom.open(wisdom_name);
    fftw_import_wisdom(get_wisdom,ifwisdom);
    ifwisdom.close();
    wise=true;
  }

  void SaveWisdom() {
    ofwisdom.open(wisdom_name);
    fftw_export_wisdom(put_wisdom,ofwisdom);
    ofwisdom.close();
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
    
  void fftNormalized(Complex *in, Complex *out=NULL) {
    Setout(in,out);
    Execute(in,out);
    for (unsigned int i=0; i < size; i++) out[i] *= norm;
  }
};

bool fftw::wise=false;
const char *fftw::wisdom_name="wisdom.txt";
ifstream fftw::ifwisdom;
ofstream fftw::ofwisdom;

class fft1d : public fftw {
  unsigned int nx;
public:  
  fft1d(unsigned int nx, int sign, Complex *in=NULL, Complex *out=NULL) 
    : fftw(nx,sign), nx(nx) {Setup(in,out);} 
  
  void Plan(Complex *in, Complex *out) {
    plan=fftw_plan_dft_1d(nx, (fftw_complex *) in, (fftw_complex *) out,
			  sign, effort);
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
    plan=fftw_plan_dft_2d(nx, ny, (fftw_complex *) in, (fftw_complex *) out,
			  sign, effort);
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

class rcfft1d : public fftw {
  unsigned int nx;
public:  
  rcfft1d(unsigned int nx, Complex *in=NULL, Complex *out=NULL) 
    : fftw(nx/2+1,-1,nx), nx(nx) {Setup(in,out);} 
  
  void Plan(Complex *in, Complex *out) {
    plan=fftw_plan_dft_r2c_1d(nx, (double *) in, (fftw_complex *) out, effort);
  }
  
  void Execute(Complex *in, Complex *out) {
    fftw_execute_dft_r2c(plan,(double *) in,(fftw_complex *) out);
  }
};
  
class crfft1d : public fftw {
  unsigned int nx;
public:  
  crfft1d(unsigned int nx, Complex *in=NULL, Complex *out=NULL) 
    : fftw(nx,nx/2+1,nx), nx(nx) {Setup(in,out);} 
  
  void Plan(Complex *in, Complex *out) {
    plan=fftw_plan_dft_c2r_1d(nx, (fftw_complex *) in, (double *) out, effort);
  }
  
  void Execute(Complex *in, Complex *out) {
    fftw_execute_dft_c2r(plan,(fftw_complex *) in,(double *) out);
  }
};
  
class rcfft2d : public fftw {
  unsigned int nx;
  unsigned int ny;
public:  
  rcfft2d(unsigned int nx, unsigned int ny, Complex *in=NULL,
	  Complex *out=NULL) 
    : fftw(nx*(ny/2+1),-1,nx*ny), nx(nx), ny(ny) {Setup(in,out);} 
  
  void Plan(Complex *in, Complex *out) {
    plan=fftw_plan_dft_r2c_2d(nx, ny, (double *) in, (fftw_complex *) out,
			      effort);
  }
  
  void Execute(Complex *in, Complex *out) {
    fftw_execute_dft_r2c(plan,(double *) in,(fftw_complex *) out);
  }
};
  
class crfft2d : public fftw {
  unsigned int nx;
  unsigned int ny;
public:  
  crfft2d(unsigned int nx, unsigned int ny, Complex *in=NULL,
	  Complex *out=NULL) 
    : fftw(nx*(ny/2+1),1,nx*ny), nx(nx), ny(ny) {Setup(in,out);} 
  
  void Plan(Complex *in, Complex *out) {
    plan=fftw_plan_dft_c2r_2d(nx, ny, (fftw_complex *) in, (double *) out,
			      effort);
  }
  
  void Execute(Complex *in, Complex *out) {
    fftw_execute_dft_c2r(plan,(fftw_complex *) in,(double *) out);
  }
};

class rcfft3d : public fftw {
  unsigned int nx;
  unsigned int ny;
  unsigned int nz;
public:  
  rcfft3d(unsigned int nx, unsigned int ny, unsigned int nz, Complex *in=NULL,
	  Complex *out=NULL) 
    : fftw(nx*ny*(nz/2+1),-1,nx*ny*nz), nx(nx), ny(ny) {Setup(in,out);} 
  
  void Plan(Complex *in, Complex *out) {
    plan=fftw_plan_dft_r2c_3d(nx, ny, nz, (double *) in, (fftw_complex *) out,
			      effort);
  }
  
  void Execute(Complex *in, Complex *out) {
    fftw_execute_dft_r2c(plan,(double *) in,(fftw_complex *) out);
  }
};
  
class crfft3d : public fftw {
  unsigned int nx;
  unsigned int ny;
  unsigned int nz;
public:  
  crfft3d(unsigned int nx, unsigned int ny, unsigned int nz, Complex *in=NULL,
	  Complex *out=NULL) 
    : fftw(nx*ny*(nz/2+1),1,nx*ny*nz), nx(nx), ny(ny) {Setup(in,out);} 
  
  void Plan(Complex *in, Complex *out) {
    plan=fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex *) in, (double *) out,
			      effort);
  }
  
  void Execute(Complex *in, Complex *out) {
    fftw_execute_dft_c2r(plan,(fftw_complex *) in,(double *) out);
  }
};
  
#endif
