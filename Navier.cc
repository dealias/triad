#include "options.h"
#include "NWave.h"
#include "Geometry.h"
#include "Partition.h"
#include "Basis.h"
#include "Polar.h"
#include "PolarBin.h"
#include "Cartesian.h"
#include "fft.h"

#include <sys/stat.h> // On the sun this must come after xstream.h

LinearityBase *Linearity;
GeometryBase *Geometry;

char *method="PS";
char *geometry="Cartesian";
char *integrator="PC";
char *linearity="BandLimited";

class NWaveVocabulary : public VocabularyBase {
public:
	char *Name() {return "N-Wave";}
	char *Abbrev() {return "nw";}
	NWaveVocabulary();
	Table<LinearityBase> *LinearityTable;
	Table<GeometryBase> *GeometryTable;
	GeometryBase *NewGeometry(char *& key) {
		return GeometryTable->Locate(key);
	}
	LinearityBase *NewLinearity(char *& key) {
		return LinearityTable->Locate(key);
	}
};

// Global variables
Real alpha=1.0;
Real beta=1.0;
Real E0=1.0;
Real U0=1.0;

Real kforce=0.0;
Real deltaf=1.0;
Real gammaf=0.0;
Real nu0=0.0;

Real force=0.0;
Real tauforce=0.0;

Real vd=0.0;

Real kL=STD_MAX;
Real kH=0.0;

Real nuL=0.0;
Real nuH=0.0;

int pL=-2;
int pH=1;

int randomIC=1;
int discrete=0;

int ngridx=0;
int ngridy=0;
int movie=0;
int truefield=0;
int weiss=0;

// (0 => evolve all modes, 1 => evolve only half of the modes).
int reality=1; // Reality condition flag 

Real mode_density;
Weight *weightBase;
WeightIndex WeightN;

const int maxbins=65536;

class Diamagnetic : public LinearityBase {
public:
	Real Frequency(const Mode& v) {return vd*v.Y()/Denominator(v);}
};

class BandLimited : public Diamagnetic {
public:
	char *Name() {return "Band-Limited";}
	
	Real Growth(const Mode& v) {
		Real k=v.K();
		Real gamma=0.0;
		if(k <= kL) gamma -= pow(k,pL)*nuL*pow(k,pL);
		if(abs(k-kforce) < 0.5*deltaf) gamma += gammaf/deltaf;
		if(k > kH) gamma -= pow(k,pH)*nuH*pow(k,pH);
		return gamma;
	}
};

class Waltz : public Diamagnetic {
public:
	char *Name() {return "Waltz";}
	
	Real Denominator(const Mode &v) {return 1.0+v.K2();}
	
	Real Growth(const Mode& v) {
		Real tempx=abs(v.X())/0.5-1.0;
		Real tempy=abs(v.Y())/0.5-1.0;
		return 0.06*(1.0-0.5*(tempx*tempx+tempy*tempy))-0.05;
	}
};

class ConstantFrequency : public BandLimited {
public:
	char *Name() {return "ConstantFrequency";}
	Real Frequency(const Mode&) {return vd;}
};

class FrequencyK : public BandLimited {
public:
	char *Name() {return "FrequencyK";}
	Real Frequency(const Mode& v) {return vd*v.K();}
};

class Convolution : public NWave {
public:
	Convolution() {}
	char *Name() {return "Convolution";}
	void NonLinearSrc(Var *source, Var *psi, double);
};

class PS : public NWave {
public:
	PS() {
		if(!reality) msg(ERROR,"Pseudospectral approximation needs reality=1");
	}
	char *Name() {return "Pseudospectral";}
	void NonLinearSrc(Var *source, Var *psi, double);
};

NWaveVocabulary::NWaveVocabulary()
{
	Vocabulary=this;
	
	VOCAB(reality,0,1);
	VOCAB(geometry,"","");
	VOCAB(linearity,"","");
	
	VOCAB(alpha,0.0,0.0);
	VOCAB(beta,0.0,0.0);
	VOCAB(E0,0.0,0.0);
	VOCAB(U0,0.0,0.0);
	
	VOCAB(kforce,0.0,STD_MAX);
	VOCAB(deltaf,0.0,STD_MAX);
	VOCAB(gammaf,0.0,STD_MAX);
	
	VOCAB(force,0.0,STD_MAX);
	VOCAB(tauforce,0.0,STD_MAX);
	
	VOCAB(vd,0.0,STD_MAX);
	
	VOCAB(kL,0.0,STD_MAX);
	VOCAB(kH,0.0,STD_MAX);
	
	VOCAB(nuL,0.0,STD_MAX);
	VOCAB(nuH,0.0,STD_MAX);
	
	VOCAB(pL,0,0);
	VOCAB(pH,0,0);
	
	VOCAB(Nx,1,INT_MAX);
	VOCAB(Ny,1,INT_MAX);
	
	VOCAB(Nr,1,INT_MAX);
	VOCAB(Nth,1,INT_MAX);
	VOCAB(krmin,0.0,DBL_MAX);
	VOCAB(krmax,0.0,DBL_MAX);
	VOCAB(kthmin,0.0,twopi);
	
	VOCAB(Nmoment,0,INT_MAX);
	VOCAB(randomIC,0,1);
	VOCAB(discrete,0,1);
	
	VOCAB(ngridx,0,INT_MAX);
	VOCAB(ngridy,0,INT_MAX);
	VOCAB(movie,0,1);
	VOCAB(truefield,0,1);
	VOCAB(weiss,0,1);
	
	GeometryTable=new Table<GeometryBase>("Geometry");
	LinearityTable=new Table<LinearityBase>("Linearity");

	METHOD(SR);
	METHOD(Convolution);
	METHOD(PS);
	
	INTEGRATOR(C_Euler);
	INTEGRATOR(I_PC);
	INTEGRATOR(C_PC);
	INTEGRATOR(E_PC);
	INTEGRATOR(CE_PC);
	INTEGRATOR(I_RK2);
	INTEGRATOR(C_RK2);
	INTEGRATOR(I_RK4);
	INTEGRATOR(C_RK4);
	INTEGRATOR(I_RK5);
	INTEGRATOR(C_RK5);
	
	PARTITION(Polar,Cartesian);
	BASIS(Cartesian);
	
	LINEARITY(BandLimited);
	LINEARITY(Waltz);
	LINEARITY(ConstantFrequency);
	LINEARITY(FrequencyK);
}

NWaveVocabulary NWave_Vocabulary;

void Convolution::NonLinearSrc(Var *source, Var *psi, double)
{
	int i;
	set(psibuffer,psi,Nmode);
	
	// Compute reflected psi's
#pragma ivdep		
	for(Var *p=psibuffer; p < psibufferR; p++) conjugate(*(p+Nmode),*p);
	
	for(i=0; i < Nmode; i++) source[i]=0.0;
		
	int q=0, nn=2*Nmode;
	for(int k=0; k < Nmode; k++) {
		Real kx=CartesianMode[k].X();
		Real ky=CartesianMode[k].Y();
		for(int p=0; p < nn; p++)	{
			Real px=CartesianMode[p].X();
			Real py=CartesianMode[p].Y();
			Cartesian mq = CartesianMode[k]-CartesianMode[p];
			if(q < nn-1 && CartesianMode[q+1] == mq) q++;
			else if(q > 1 && CartesianMode[q-1] == mq) q--;
			else for(q=0; q < nn && CartesianMode[q] != mq; q++);
			if(q < nn) source[k] += (kx*py-ky*px)*(px*px+py*py)*
				psibuffer[p]*psibuffer[q];
		}
	}

#pragma ivdep	
	for(i=0; i < Nmode; i++) source[i] *= kfactor[i];
	if(Nmoment) ComputeMoments(source,psi);
	ConstantForcing(source,t);
}


void PS::NonLinearSrc(Var *source, Var *psi, double)
{
#if COMPLEX
	const int bitreverse=1; // If 1, use faster bit-reversed FFT's.
	int i;
	
#pragma ivdep	
	for(i=0; i < Nmode; i++) {
		Real kx=CartesianMode[i].X();
		source[i].re=-psi[i].im*kx;
		source[i].im=psi[i].re*kx;
	}
	CartesianPad(psix,source);
 	crfft2dT_sym(psix,log2Nxb,log2Nyb,1,1.0,-bitreverse);

#pragma ivdep	
	for(i=0; i < Nmode; i++) source[i] *= knorm2[i];
	CartesianPad(vort,source);
	crfft2dT(vort,log2Nxb,log2Nyb,1,1.0,-bitreverse);

#pragma ivdep	
	for(i=0; i < Nmode; i++) {
		Real ky=CartesianMode[i].Y();
		source[i].re=-psi[i].im*ky;
		source[i].im=psi[i].re*ky;
	}
	CartesianPad(psiy,source);
	crfft2dT(psiy,log2Nxb,log2Nyb,1,1.0,-bitreverse);

#if 0 // Compute x-space velocity increments; requires bitreverse=0.
	Real *v2;
	// Strictly speaking, v2 should be divided by (Nxb*Nyb)^2 afterwards
	Real psix0=psix[0].re;
	Real psiy0=psiy[0].re;
	for(int j=0; j < Nyb; j++) {
		int jN=Nxb1*j;
#pragma ivdep	
		for(i=0; i < Nxb; i++) {
			Real vx1=psiy0-psiy[i+jN].re;
			Real vy1=psix[i+jN].re-psix0;
			*(v2++)=vx1*vx1+vy1*vy1;
			vx1=psiy0-psiy[i+jN].im;
			vy1=psix[i+jN].im-psix0;
			*(v2++)=vx1*vx1+vy1*vy1;
		}
	}
#endif	

#pragma ivdep	
	for(i=0; i < nfft; i++) {
		psiy[i].re *= vort[i].re;
		psiy[i].im *= vort[i].im;
	}

#pragma ivdep	
	for(i=0; i < Nmode; i++) source[i] *= knorm2[i];
	CartesianPad(vort,source);
	crfft2dT(vort,log2Nxb,log2Nyb,1,1.0,-bitreverse);

#pragma ivdep	
	for(i=0; i < nfft; i++) {
		psiy[i].re -= psix[i].re*vort[i].re;
		psiy[i].im -= psix[i].im*vort[i].im;
	}

	rcfft2dT(psiy,log2Nxb,log2Nyb,-1,1.0,bitreverse);
	CartesianUnPad(source,psiy);
	
#pragma ivdep	
	for(i=0; i < Nmode; i++) source[i] *= kfactor[i];
	if(Nmoment) ComputeMoments(source,psi);
	ConstantForcing(source,t);
#else	
	msg(ERROR,"Pseudospectral approximation requires COMPLEX=1");
#endif
}
	
void Basis<Cartesian>::Initialize()
{
	knorm2=new Real[Nmode];
	kfactor=new Real[Nmode];
	
	if(strcmp(Problem->Abbrev(),"PS") == 0) {
		cout << endl << "ALLOCATING FFT BUFFERS (" << Nxb << " x " << Nyp
			 << ")." << endl;
		psix=new Var[nfft];
		psiy=new Var[nfft];
		vort=new Var[nfft];
		Real scale=Nxb*Nyb;
		for(int k=0; k < Nmode; k++) {
			knorm2[k]=mode[k].K2();
			kfactor[k]=1.0/(scale*Normalization(k));
		}
	} else {
		psibuffer=new Var[n];
		psibufferR=(reality ? psibuffer+Nmode : psibuffer);
		for(int k=0; k < Nmode; k++) 
			kfactor[k]=1.0/Normalization(k);
	}
}

Real force_re(const Mode& v) 
{
	Real k=v.K();
	if(abs(k-kforce) < 0.5*deltaf) return force;
	else return 0.0;
}

static Real equilibrium(int i)
{
	return 0.5*mode_density*mode_density/(alpha+beta*Geometry->K2(i));
}

static ifstream ftin;
static ofstream fparam,fevt,fyvt,ft,favgy,fprolog;
static oxstream fpsi,fweiss;

Real continuum_factor=1.0;
static strstream *avgyre,*avgyim;
static int tcount=0;
static Complex *xcoeff, *ycoeff;
static Real *norm_factor;

void NWave::InitialConditions()
{
	int i,n;
	
	krmin2=krmin*krmin;
	mode_density=(strcmp(method,"SR") == 0 ? 1.0/krmin : 1.0);
		
	Geometry=NWave_Vocabulary.NewGeometry(geometry);
	if(!Geometry->Valid(Problem->Abbrev()))
		msg(ERROR,"Geometry \"%s\" is incompatible with method \"%s\"",
			Geometry->Name(),Problem->Abbrev());
	Linearity=NWave_Vocabulary.NewLinearity(linearity);
	Nmode=Geometry->Create();
	ny=Nmode*(1+Nmoment);
	Ntotal=Geometry->TotalNumber();
	NmodeR=Ntotal-Nmode;
	y=new Var[ny];
	nu=new Nu[Nmode];
	forcing=new Real[Nmode];
	
	for(i=0; i < Nmode; i++) {
		Real norm=1.0/sqrt(Geometry->Normalization(i));
		nu[i]=Geometry->Linear(i);
		y[i]=sqrt(2.0*equilibrium(i))*norm;
		forcing[i]=Geometry->Forcing(i)*norm;
	}
	
	int nindependent=Geometry->IndependentNumber();
	// Randomize the initial conditions.	
	if(randomIC) {
		for(i=0; i < nindependent; i++) {
			Var w;
			crand_gauss(&w);
			y[i] *= w;
		}
	}
	
	// If reality condition is not explicitly enforced and the number of
	// angular bins is even, force the initial conditions to respect reality. 
#pragma ivdep		
	for(i=nindependent; i < Nmode; i++) y[i]=conj(y[i-nindependent]);
	
	if(restart) {
		Real t0;
		ftin.open(Vocabulary->FileName(dirsep,"t"));
		while(ftin >> t0, ftin.good()) tcount++;
		ftin.close();
	}
	
	open_output(fparam,dirsep,"param",0);
	open_output(fevt,dirsep,"evt");
	open_output(fyvt,dirsep,"yvt");
	open_output(ft,dirsep,"t");
	open_output(fprolog,dirsep,"prolog");
	if(movie) {
		open_output(fpsi,dirsep,"psi");
		if(weiss) open_output(fweiss,dirsep,"weiss");
	}
	
	if(Nmoment) {
		avgyre=new strstream[Nmoment];
		avgyim=new strstream[Nmoment];
	}
	
	if(!restart) remove_dir(Vocabulary->FileName(dirsep,"avgy*"));
	for(n=0; n < Nmoment; n++) {
		if(!restart) {
			strstream buf;
			buf << "avgy" << n << ends;
			mkdir(Vocabulary->FileName(dirsep,buf.str()),0xFFFF);
			errno=0;
		}
		avgyre[n] << "y.re^" << n << ends;
		avgyim[n] << "y.im^" << n << ends;
	}
	
	Vocabulary->GraphicsDump(fparam);
	fparam.close();
	
	if(movie && strcmp(method,"SR") == 0) {
		if(discrete && truefield) {
			set_fft_parameters();
			psix=new Var[nfft];
			norm_factor=new Real[Nmode];
			for(int m=0; m < Nmode; m++) {
				norm_factor[m]=sqrt(Geometry->Normalization(m)/Linearity->
									Denominator(Geometry->ModeOf(m)));
			}
		} else {
			int m;
			Real k0=REAL_MAX;
			Real Dk0=REAL_MAX;
			for(i=0; i < Nmode; i++) {
				Real k=Geometry->K(i);
				if(k < k0) {
					k0=k;
					Dk0=Geometry->Area(i);
				}
			}
			if(!ngridx) ngridx=Nx;
			if(!ngridy) ngridy=Ny;
			Real L=twopi/(k0*Dk0);
			
			xcoeff=new Complex [ngridx*Nmode];
			for(m=0; m < Nmode; m++) {
				Real Dkinv=1.0/Geometry->Area(m);
				Real kx=Geometry->X(m);
				Real norm=sqrt(Geometry->Normalization(m)/
							   Linearity->Denominator(Geometry->ModeOf(m)));

				for(i=0; i < ngridx; i++) {
					Complex *p=xcoeff+i*Nmode;
					Real X=i*L/ngridx;
					p[m]=expi(kx*X*Dkinv)*norm;
				}				
			}
			
			ycoeff=new Complex [ngridy*Nmode];
			for(m=0; m < Nmode; m++) {
				Real ky=Geometry->Y(m);
				Real Dkinv=1.0/Geometry->Area(m);
				for(int j=0; j < ngridy; j++) {
					Complex *q=ycoeff+j*Nmode;
					Real Y=j*L/ngridy;
					q[m]=expi(ky*Y*Dkinv);
				}
			}
		}
	}
}

static Real K(int i) {return Geometry->K(i);}
static Real Th(int i) {return Geometry->Th(i);}
static Real Area(int i) {return Geometry->Area(i);}
static Real Normalization(int i) {return Geometry->Normalization(i);}
static Real nu_re(int i) {Complex nuC=nu[i]; return nuC.re;}
static Real nu_im(int i) {Complex nuC=nu[i]; return nuC.im;}

void NWave::Initialize()
{
	int i;
	
	out_function(fprolog,K,"K",Nmode);
	out_function(fprolog,Th,"Th",Nmode);
	out_function(fprolog,Area,"Area",Nmode);
	out_function(fprolog,nu_re,"nu.re",Nmode);
	out_function(fprolog,nu_im,"nu.im",Nmode);
	out_curve(fprolog,forcing,"f",Nmode);
	out_function(fprolog,equilibrium,"equil",Nmode);
	out_function(fprolog,Normalization,"normalization",Nmode);
	fprolog.flush();

	fevt << "#   t\t\t E\t\t Z\t\t P" << endl;

	// Initialize time integrals to zero.
	for(i=Nmode; i < ny; i++) y[i]=0.0;
}

void NWave::ComputeInvariants(Var *u, int Nmode, Real& E, Real& Z, Real& P)
{
	Real Ek,Zk,Pk,k2;
	E=Z=P=0.0;
	for(int i=0; i < Nmode; i++) {
		k2=Geometry->K2(i);
		Ek=Geometry->Normalization(i)*abs2(u[i])*Geometry->Area(i);
		Zk=k2*Ek;
		Pk=k2*Zk;
		E += Ek;
		Z += Zk;
		P += Pk;
	}
	
	Real factor=(reality ? 1.0: 0.5)*continuum_factor;
	E *= factor;
	Z *= factor;
	P *= factor;
}	

void out_field(oxstream& fout, Real *psir)
{
	fout << Nxb << Nyb << 1;
	int factor=2*(Nxb1-1);
	for(int j=Nyb-1; j >= 0; j--) {
		int jN=factor*(j/2)+j;
		for(int i=0; i < Nxb; i++)
			fout << (float) psir[2*i+jN];
	}
	fout.flush();
}

void NWave::Output(int)
{
	Real E,Z,P;
	int n;
	
	ComputeInvariants(y,Nmode,E,Z,P);
	
	fevt << t << "\t" << E << "\t" << Z << "\t" << P << endl;
	
	out_real(fyvt,y,"y.re","y.im",Nmode);
	fyvt.flush();
	
	Var *yavg=y+Nmode;
	for(n=0; n < Nmoment; n++) {
		strstream buf;
		buf << "avgy" << n << dirsep << "t" << tcount << ends;
		open_output(favgy,dirsep,buf.str(),0);
		out_curve(favgy,t,"t");
		out_real(favgy,yavg+Nmode*n,avgyre[n].str(),avgyim[n].str(),Nmode);
		favgy.close();
	}
	
	if(movie) {
		if(strcmp(method,"SR") == 0 && !(discrete && truefield)) {
			fpsi << ngridx << ngridy << 1;
			for(int j=ngridy-1; j >= 0; j--) {
				Complex *q=ycoeff+j*Nmode;
				for(int i=0; i < ngridx; i++) {
					Complex *p=xcoeff+i*Nmode;
					Real sum=0.0;
					for(int m=0; m < Nmode; m++) {
						Complex pq=p[m]*q[m];
						sum += pq.re*y[m].re-pq.im*y[m].im;
					}
					fpsi << (float) sum;
				}
			}
			fpsi.flush();
		} else {
			if(discrete) DiscretePad(psix,y,norm_factor);
			else CartesianPad(psix,y);
			
			crfft2dT(psix,log2Nxb,log2Nyb,1);
			Real ninv=1.0/(Nxb*Nyb);
			for(int i=0; i < nfft; i++) psix[i] *= ninv;
			
			out_field(fpsi,(Real *) psix);
		}
		if(!fpsi) msg(ERROR, "Cannot write to movie file psi");
	
		if(weiss && !discrete) {
			int i;
			for(i=0; i < Nmode; i++) {
				Real kx=CartesianMode[i].X();
				vort[i]=y[i]*kx*kx;
			}
			CartesianPad(psix,vort);
			crfft2dT(psix,log2Nxb,log2Nyb,1);
		
			for(i=0; i < Nmode; i++) {
				Real ky=CartesianMode[i].Y();
				vort[i]=y[i]*ky*ky;
			}
			CartesianPad(psiy,vort);
			crfft2dT(psiy,log2Nxb,log2Nyb,1);
			
			for(i=0; i < nfft; i++) psiy[i] *= psix[i];
				
			for(i=0; i < Nmode; i++) {
				Real kx=CartesianMode[i].X();
				Real ky=CartesianMode[i].Y();
				vort[i]=y[i]*kx*ky;
			}
			CartesianPad(psix,vort);
			crfft2dT(psix,log2Nxb,log2Nyb,1);
			
			Real ninv=1.0/(Nxb*Nyb);
			for(i=0; i < nfft; i++)
				psix[i]=(psix[i]*psix[i]-psiy[i])*ninv;
			
			out_field(fweiss,(Real *) psix);
		}
	}
	
	tcount++;
	ft << t << endl;
}

void ForcingAt(int i, Real &force)
{
	force=force_re(Geometry->ModeOf(i));
	return;
}
