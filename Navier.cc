#include "options.h"
#include "NWave.h"
#include "Polar.h"
#include "Cartesian.h"
#include "fft.h"

#include <sys/stat.h>

char *NWaveVocabulary::Name() {return "N-Wave";}
char *NWaveVocabulary::Abbrev() {return "nw";}

char *method="PS";
char *geometry="Cartesian";
char *integrator="PC";

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

Real kL=STD_MAX;
Real kH=0.0;

Real nuL=0.0;
Real nuH=0.0;

int pL=-2;
int pH=1;

int randomIC=1;
int ngridx=0;
int ngridy=0;
int movie=0;
int truefield=0;

Real Dx=1.0;

NWaveVocabulary::NWaveVocabulary()
{
	Vocabulary=this;
	
	VOCAB(reality,0,1);
	VOCAB(geometry,"","");
	
	VOCAB(alpha,0.0,0.0);
	VOCAB(beta,0.0,0.0);
	VOCAB(E0,0.0,0.0);
	VOCAB(U0,0.0,0.0);
	
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
	
	VOCAB(randomIC,0,1);
	VOCAB(Nmoment,0,INT_MAX);
	
	VOCAB(Nx,1,INT_MAX);
	VOCAB(Ny,1,INT_MAX);
	
	VOCAB(Nr,1,INT_MAX);
	VOCAB(Nth,1,INT_MAX);
	VOCAB(krmin,0.0,DBL_MAX);
	VOCAB(krmax,0.0,DBL_MAX);
	VOCAB(kthmin,0.0,twopi);
	
	VOCAB(ngridx,0,INT_MAX);
	VOCAB(ngridy,0,INT_MAX);
	VOCAB(movie,0,1);
	VOCAB(truefield,0,1);
	
	VOCAB(Dx,0.0,STD_MAX);
	
	GeometryTable=new Table<GeometryBase>("Geometry");

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
}

NWaveVocabulary NWave_Vocabulary;

Real force_re(const Polar& v) 
{
	Real k=v.K();
	if(abs(k-kforce) < 0.5*deltaf) return force;
	else return 0.0;
}

inline Real growth(const Polar& v) 
{
	Real k=v.K();
	Real gamma=0.0;
	if(k <= kL) gamma -= pow(k,pL)*nuL*pow(k,pL);
	if(abs(k-kforce) < 0.5*deltaf) gamma += gammaf/deltaf;
	if(k > kH) gamma -= pow(k,pH)*nuH*pow(k,pH);
	return gamma;
}

inline Real frequency(const Polar&)
{
	return 0.0;
}

Real linearity_re(const Polar& v)
{
	return -growth(v);
}

Real linearity_im(const Polar& v)
{
	return frequency(v);
}

static Real equilibrium(int i)
{
	Real k=Geometry->K(i);
	return 0.5/(alpha+beta*k*k);
}

static ifstream ftin;
static ofstream fparam,fevt,fyvt,ft,favgy,fprolog;
static oxstream fpsi;

Real continuum_factor=1.0;
typedef char Avgylabel[20];
static Avgylabel *avgyre,*avgyim;
static int tcount=0;
static Complex *xcoeff, *ycoeff;

void NWave::InitialConditions()
{
	int i,n;
	
	Geometry=NWave_Vocabulary.NewGeometry(geometry);
	if(!Geometry->Valid(Problem->Abbrev()))
		msg(ERROR,"Geometry \"%s\" is incompatible with method \"%s\"",
			Geometry->Name(),Problem->Abbrev());
	Npsi=Geometry->Create();
	ny=Npsi*(1+Nmoment);
	Ntotal=Geometry->TotalNumber();
	NpsiR=Ntotal-Npsi;
	y=new Var[ny];
	nu=new Nu[Npsi];
	forcing=new Real[Npsi];
	
	for(i=0; i < Npsi; i++) {
		Real norm=1.0/sqrt(Geometry->Normalization(i));
		nu[i]=Geometry->Linearity(i);
		y[i]=sqrt(2.0*equilibrium(i))*norm;
		forcing[i]=Geometry->Forcing(i)*norm;
	}
	
	// Randomize the initial conditions.	
	// If reality condition is not explicitly enforced and the number of
	// angular bins is even, force the initial conditions to respect reality. 
	if(randomIC) {
		int nindependent=Geometry->IndependentNumber();
		for(i=0; i < nindependent; i++) {
			Var w;
			crand_gauss(&w);
			y[i] *= w;
		}
#pragma ivdep		
		for(i=nindependent; i < Npsi; i++) y[i]=conj(y[i-nindependent]);
	}
	
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
	if(movie) open_output(fpsi,dirsep,"psi");
	
	if(Nmoment) {
		avgyre=new Avgylabel[Nmoment];
		avgyim=new Avgylabel[Nmoment];
	}
	
	if(!restart) remove_dir(Vocabulary->FileName(dirsep,"avgy*"));
	for(n=0; n < Nmoment; n++) {
		if(!restart) {
			static char tempbuffer[30];
			sprintf(tempbuffer,"avgy%d",n);
			mkdir(Vocabulary->FileName(dirsep,tempbuffer),0xFFFF);
			errno=0;
		}
		sprintf(avgyre[n],"y.re^%d",n);
		sprintf(avgyim[n],"y.im^%d",n);
	}
	
	Vocabulary->GraphicsDump(fparam);
	fparam.close();
	
	if(movie && strcmp(method,"SR") == 0) {
		if(discrete && truefield) {
			set_fft_parameters();
			psix=new Var[nfft];
		} else {
			int m;
			Real k0=REAL_MAX;
			Real Dk0=REAL_MAX;
			for(i=0; i < Npsi; i++) {
				Real k=Geometry->K(i);
				if(k < k0) {
					k0=k;
					Dk0=Geometry->Area(i);
				}
			}
			if(!ngridx) ngridx=Nx;
			if(!ngridy) ngridy=Ny;
			Real L=twopi/(k0*Dk0);
			
			xcoeff=new Complex [ngridx*Npsi];
			for(m=0; m < Npsi; m++) {
				Real DkDxinv=1.0/(Dx*Geometry->Area(m));
				Real kx=Geometry->X(m);
				Real norm=sqrt(Normalization(m)/Geometry->K2(m));
				for(i=0; i < ngridx; i++) {
					Complex *p=xcoeff+i*Npsi;
					Real X=i*L/ngridx;
					p[m]=expi(kx*X*DkDxinv)*norm;
				}				
			}
			
			ycoeff=new Complex [ngridy*Npsi];
			for(m=0; m < Npsi; m++) {
				Real ky=Geometry->Y(m);
				Real DkDxinv=1.0/(Dx*Geometry->Area(m));
				for(int j=0; j < ngridy; j++) {
					Complex *q=ycoeff+j*Npsi;
					Real Y=j*L/ngridy;
					q[m]=expi(ky*Y*DkDxinv);
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
	
	out_function(fprolog,K,"K",Npsi);
	out_function(fprolog,Th,"Th",Npsi);
	out_function(fprolog,Area,"Area",Npsi);
	out_function(fprolog,nu_re,"nu.re",Npsi);
	out_function(fprolog,nu_im,"nu.im",Npsi);
	out_curve(fprolog,forcing,"f",Npsi);
	out_function(fprolog,equilibrium,"equil",Npsi);
	out_function(fprolog,Normalization,"normalization",Npsi);
	fprolog.flush();

	fevt << "#   t\t\t E\t\t Z\t\t P" << endl;

	// Initialize time integrals to zero.
	for(i=Npsi; i < ny; i++) y[i]=0.0;
}

void compute_invariants(Var *y, int Npsi, Real& E, Real& Z, Real& P)
{
	Real Ek,Zk,Pk,k2;
	E=Z=P=0.0;
	for(int i=0; i < Npsi; i++) {
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
	
	compute_invariants(y,Npsi,E,Z,P);
	
	fevt << t << "\t" << E << "\t" << Z << "\t" << P << endl;
	
	out_real(fyvt,y,"y.re","y.im",Npsi);
	fyvt.flush();
	
	Var *yavg=y+Npsi;
	for(n=0; n < Nmoment; n++) {
		static char tempbuffer[30];
		sprintf(tempbuffer,"avgy%d%st%d",n,dirsep,tcount);
		open_output(favgy,dirsep,tempbuffer,0);
		out_curve(favgy,t,"t");
		out_real(favgy,yavg+Npsi*n,avgyre[n],avgyim[n],Npsi);
		favgy.close();
	}
	
	if(movie) {
		if(strcmp(method,"SR") == 0 && !(discrete && truefield)) {
			fpsi << ngridx << ngridy << 1;
			for(int j=ngridy-1; j >= 0; j--) {
				Complex *q=ycoeff+j*Npsi;
				for(int i=0; i < ngridx; i++) {
					Complex *p=xcoeff+i*Npsi;
					Real sum=0.0;
					for(int m=0; m < Npsi; m++) {
						Complex pq=p[m]*q[m];
						sum += pq.re*y[m].re-pq.im*y[m].im;
					}
					fpsi << (float) sum;
				}
			}
		} else {
			fpsi << Nxb << Nyb << 1;
			if(discrete) DiscretePad(psix,y);
			else CartesianPad(psix,y);
			
			crfft2dT(psix,log2Nxb,log2Nyb,1);
		
			Real *psir=(Real *) psix;
			Real ninv=1.0/(Nxb*Nyb);
			int factor=2*(Nxb1-1);
			for(int j=Nyb-1; j >= 0; j--) {
				int jN=factor*(j/2)+j;
				for(int i=0; i < Nxb; i++)
					fpsi << (float) (psir[2*i+jN]*ninv);
			}
		}
		fpsi.flush();
		if(!fpsi) msg(ERROR, "Cannot write to movie file psi");
	}
	
	tcount++;
	ft << t << endl;
}

void LinearityAt(int i,Real& nu)
{
	Polar v=Polar(Geometry->K(i),Geometry->Th(i));
	nu=linearity_re(v);
	return;
}

void LinearityAt(int i, Complex& nu)
{
	Polar v=Polar(Geometry->K(i),Geometry->Th(i));
	nu=linearity_re(v)+I*linearity_im(v);
	return;
}

void ForcingAt(int i, Real &force)
{
	force=force_re(Polar(Geometry->K(i),Geometry->Th(i)));
	return;
}
