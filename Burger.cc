#include <sys/stat.h>

#include "options.h"
#include "Cartesian1.h"
#include "Burger.h"
#include "fft.h"

LinearityBase *Linearity;

char *BurgerVocabulary::Name() {return "Burger";}
char *BurgerVocabulary::Abbrev() {return "burger";}

char *method="PS";
char *geometry="Cartesian1";
char *integrator="PC";
char *linearity="BandLimited";

int Nmode;   // number of explictly evolved modes
int NmodeR;  // number of reflected modes
int Ntotal; // total number of (evolved+reflected) modes

// (0 => evolve all modes, 1 => evolve only half of the modes).
int reality=1; // Reality condition flag 

// Global variables
Real alpha=1.0;
Real E0=1.0;

Real kforce=0.0;
Real deltaf=1.0;
Real gammaf=0.0;
Real nu0=0.0;

Real force=0.0;
Real tauforce=0.0;

Real kL=STD_MAX;
Real kH=0.0;

Real nuL=0.0;
Real nuH=0.0;

int pL=-2;
int pH=1;

int randomIC=1;

static Real mode_density;
Real scale;

Nu *nu,*nu_inv;
Real *nuR_inv,*nuI;
Real *forcing;

static Var random_factor=0.0;
static Real last_t=-REAL_MAX;
static Var *ux, *uconv;

void ConstantForcing(Var *source, double t)
{
	if(t-last_t > tauforce) {last_t=t; crand_gauss(&random_factor);}
#pragma ivdep
	for(int k=0; k < Nmode; k++) source[k] += forcing[k]*random_factor;
}

class BandLimited : public LinearityBase {
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

BurgerVocabulary::BurgerVocabulary()
{
	Vocabulary=this;
	
	VOCAB(reality,0,1);
	VOCAB(geometry,"","");
	VOCAB(linearity,"","");
	
	VOCAB(alpha,0.0,0.0);
	VOCAB(E0,0.0,0.0);
	
	VOCAB(kforce,0.0,STD_MAX);
	VOCAB(deltaf,0.0,STD_MAX);
	VOCAB(gammaf,0.0,STD_MAX);
	
	VOCAB(force,0.0,STD_MAX);
	VOCAB(tauforce,0.0,STD_MAX);
	
	VOCAB(kL,0.0,STD_MAX);
	VOCAB(kH,0.0,STD_MAX);
	
	VOCAB(nuL,0.0,STD_MAX);
	VOCAB(nuH,0.0,STD_MAX);
	
	VOCAB(pL,0,0);
	VOCAB(pH,0,0);
	
	VOCAB(Nx,1,INT_MAX);
	
	VOCAB(krmin,0.0,DBL_MAX);
	
	VOCAB(Nmoment,0,INT_MAX);
	VOCAB(randomIC,0,1);
	
	GeometryTable=new Table<GeometryBase>("Geometry");
	LinearityTable=new Table<LinearityBase>("Linearity");

	METHOD(PS);
	BASIS(Cartesian1);
	
	LINEARITY(BandLimited);
}

void PS::NonLinearSrc(Var *source, Var *u, double)
{
#if COMPLEX
	int i;
	
	for(i=0; i < Nmode; i++) {
		Real kx=Cartesian1Mode[i].X();
		source[i].re=-u[i].im*kx;
		source[i].im=u[i].re*kx;
	}
	Cartesian1Pad(ux,source);
	crfft(ux,log2Nxb,1);

	Cartesian1Pad(uconv,u);
	crfft(uconv,log2Nxb,1);
	
	for(i=0; i < nfft; i++) {
		uconv[i].re *= scale*ux[i].re;
		uconv[i].im *= scale*ux[i].im;
	}

	rcfft(uconv,log2Nxb,-1);
	Cartesian1UnPad(source,uconv);
	
//	if(Nmoment) ComputeMoments(source,u);
	ConstantForcing(source,t);
#else	
	msg(ERROR,"Pseudospectral approximation requires COMPLEX=1");
#endif
}

void Burger::LinearSrc(Var *source, Var *u, double)
{
#pragma ivdep
	for(int k=0; k < Nmode; k++) source[k] -= nu[k]*u[k];
}

BurgerVocabulary Burger_Vocabulary;

Real force_re(const Mode& v) 
{
	Real k=v.K();
	if(abs(k-kforce) < 0.5*deltaf) return force;
	else return 0.0;
}

static Real equilibrium(int)
{
	return 0.5*mode_density*mode_density/alpha;
}

static ifstream ftin;
static ofstream fparam,fevt,fyvt,ft,favgy,fprolog;

Real continuum_factor=1.0;
static strstream *avgyre,*avgyim;
static int tcount=0;

void Burger::InitialConditions()
{
	int i,n;
	
	krmin2=krmin*krmin;
	mode_density=(strcmp(method,"SR") == 0 ? 1.0/krmin : 1.0);
		
	Geometry=Burger_Vocabulary.NewGeometry(geometry);
	if(!Geometry->Valid(Problem->Abbrev()))
		msg(ERROR,"Geometry \"%s\" is incompatible with method \"%s\"",
			Geometry->Name(),Problem->Abbrev());
	Linearity=Burger_Vocabulary.NewLinearity(linearity);
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
}

static Real K(int i) {return Geometry->K(i);}
static Real Th(int i) {return Geometry->Th(i);}
static Real Area(int i) {return Geometry->Area(i);}
static Real Normalization(int i) {return Geometry->Normalization(i);}
static Real nu_re(int i) {Complex nuC=nu[i]; return nuC.re;}
static Real nu_im(int i) {Complex nuC=nu[i]; return nuC.im;}

void Burger::Initialize()
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

void compute_invariants(Var *y, int Nmode, Real& E)
{
	Real Ek,k2;
	E=0.0;
	for(int i=0; i < Nmode; i++) {
		k2=Geometry->K2(i);
		Ek=Geometry->Normalization(i)*abs2(y[i])*Geometry->Area(i);
		E += Ek;
	}
	
	Real factor=(reality ? 1.0: 0.5)*continuum_factor;
	E *= factor;
}	

void Burger::Output(int)
{
	Real E;
	int n;
	
	compute_invariants(y,Nmode,E);
	
	fevt << t << "\t" << E << endl;
	
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
	
	tcount++;
	ft << t << endl;
}

void ForcingAt(int i, Real &force)
{
	force=force_re(Geometry->ModeOf(i));
	return;
}

void display_invariants(Real E)
{
	cout << "Energy = " << E << newl;
}

void Burger::FinalOutput()
{
	Real E;
	int i;
	const int Nlimit=128;
	
	cout << newl << "FINAL VALUES:" << newl << endl;
	if(Nmode <= Nlimit || verbose > 1) {
		for(i=0; i < Nmode; i++) cout << "u[" << i << "] = " << y[i] << newl;
		cout << endl;
	}
	
	compute_invariants(y,Nmode,E);
	display_invariants(E);
	
	if(Nmoment > 2 && t) {
		cout << newl << "AVERAGED VALUES:" << newl << endl;
// We overwrite y+Nmode here, since it is no longer needed.
		Var *y2=y+3*Nmode;
		for(i=0; i < Nmode; i++) y2[i] = (real(y2[i])+imag(y2[i]))/t;
		
		if(Nmode <= Nlimit || verbose > 1) {
			for(i=0; i < Nmode; i++) {
				Real y2avg=real(y2[i]);
				cout << "|u|^2 [" << i << "] = " << y2avg << newl;
			}
			cout << endl;
		}
		
		for(i=0; i < Nmode; i++) y2[i] = sqrt(real(y2[i]));
		compute_invariants(y2,Nmode,E);
		display_invariants(E);
	}
}

void Basis<Cartesian1>::Initialize()
{
	if(strcmp(Problem->Abbrev(),"PS") == 0) {
		cout << endl << "ALLOCATING FFT BUFFERS (" << Nxb << ")." << endl;
		ux=new Var[nfft];
		uconv=new Var[nfft];
		scale=Nxb;
	}
}


