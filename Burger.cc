#include "options.h"
#include "NWave.h"
#include "Geometry.h"
#include "Basis.h"
#include "Cartesian1.h"
#include "fft.h"

#include <sys/stat.h> // On the sun this must come after xstream.h

LinearityBase *Linearity;
GeometryBase *Geometry;

char *method="PS";
char *geometry="Cartesian1";
char *integrator="PC";
char *linearity="BandLimited";

class BurgerVocabulary : public VocabularyBase {
public:
	char *Name() {return "Burger's Turbulence";}
	char *Abbrev() {return "burger";}
	BurgerVocabulary();
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

// (0 => evolve all modes, 1 => evolve only half of the modes).
int reality=1; // Reality condition flag 

Real mode_density;
static Real mscale;
static Var *ux,*uconv;

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

class PS : public NWave {
public:
	PS() {
		if(!reality) msg(ERROR,"Pseudospectral approximation needs reality=1");
	}
	char *Name() {return "Pseudospectral";}
	void NonLinearSrc(Var *source, Var *psi, double);
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
	
	BASIS(Cartesian1);
	
	LINEARITY(BandLimited);
}

BurgerVocabulary Burger;

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
		uconv[i].re *= mscale*ux[i].re;
		uconv[i].im *= mscale*ux[i].im;
	}

	rcfft(uconv,log2Nxb,-1);
	Cartesian1UnPad(source,uconv);
	
//	if(Nmoment) ComputeMoments(source,u);
	ConstantForcing(source,t);
#else	
	msg(ERROR,"Pseudospectral approximation requires COMPLEX=1");
#endif
}

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

void NWave::InitialConditions()
{
	int i,n;
	
	krmin2=krmin*krmin;
	mode_density=(strcmp(method,"SR") == 0 ? 1.0/krmin : 1.0);
		
	Geometry=Burger.NewGeometry(geometry);
	if(!Geometry->Valid(Problem->Abbrev()))
		msg(ERROR,"Geometry \"%s\" is incompatible with method \"%s\"",
			Geometry->Name(),Problem->Abbrev());
	Linearity=Burger.NewLinearity(linearity);
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

void NWave::ComputeInvariants(Var *y, int Nmode, Real& E, Real& Z, Real& P)
{
	Real Ek,Zk,Pk,k2;
	E=Z=P=0.0;
	for(int i=0; i < Nmode; i++) {
		k2=Geometry->K2(i);
		Ek=Geometry->Normalization(i)*abs2(y[i])*Geometry->Area(i);
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
	
	tcount++;
	ft << t << endl;
}

void ForcingAt(int i, Real &force)
{
	force=force_re(Geometry->ModeOf(i));
	return;
}

void Basis<Cartesian1>::Initialize()
{
	if(strcmp(Problem->Abbrev(),"PS") == 0) {
		cout << endl << "ALLOCATING FFT BUFFERS (" << Nxb << ")." << endl;
		ux=new Var[nfft];
		uconv=new Var[nfft];
		Real scale=1.0/Nxb;
		mscale=-scale;
	}
}


