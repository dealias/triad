#include "NWave.h"
#include "Polar.h"
#include "Cartesian.h"

#include <sys/stat.h>

char *NWave::Name() {return "N-Wave";}
char *NWave::Abbrev() {return "nw";}

//char *approximation="SR";
//char *geometry="Polar";

char *approximation="PS";
char *geometry="Cartesian";
char *integrator="PC";

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

Real kL=0.9*DBL_MAX;
Real kH=0.0;

Real nuL=0.0;
Real nuH=0.0;

int pL=-2;
int pH=1;

int randomIC=1;

NWave::NWave()
{
	Problem=this;
	GeometryProblem=this;
	
	VOCAB(reality,0,1);
	VOCAB(geometry,"","");
	
	VOCAB(alpha,0.0,0.0);
	VOCAB(beta,0.0,0.0);
	VOCAB(E0,0.0,0.0);
	VOCAB(U0,0.0,0.0);
	
	VOCAB(kforce,0.0,REAL_MAX);
	VOCAB(deltaf,0.0,REAL_MAX);
	VOCAB(gammaf,0.0,REAL_MAX);
	
	VOCAB(force,0.0,REAL_MAX);
	VOCAB(tauforce,0.0,REAL_MAX);
	
	VOCAB(kL,0.0,REAL_MAX);
	VOCAB(kH,0.0,REAL_MAX);
	
	VOCAB(nuL,0.0,REAL_MAX);
	VOCAB(nuH,0.0,REAL_MAX);
	
	VOCAB(pL,0,0);
	VOCAB(pH,0,0);
	
	VOCAB(randomIC,0,1);
	
	VOCAB(Nx,1,INT_MAX);
	VOCAB(Ny,1,INT_MAX);
	
	VOCAB(Nr,1,INT_MAX);
	VOCAB(Nth,1,INT_MAX);
	VOCAB(krmin,0.0,DBL_MAX);
	VOCAB(krmax,0.0,DBL_MAX);
	VOCAB(kthmin,0.0,twopi);
	
	GeometryTable=new Table<GeometryBase>("Geometry",GeometryCompare,
										  GeometryKeyCompare);
	APPROXIMATION(SR);
	APPROXIMATION(Convolution);
	APPROXIMATION(PS);
	
	INTEGRATOR(C_Euler);
	INTEGRATOR(I_PC);
	INTEGRATOR(C_PC);
	INTEGRATOR(E_PC);
	INTEGRATOR(CE_PC);
	INTEGRATOR(I_RK2);
	INTEGRATOR(C_RK2);
	INTEGRATOR(E_RK2);
	INTEGRATOR(C_RK4);
	INTEGRATOR(E_RK4);
	INTEGRATOR(C_RK5);
	INTEGRATOR(E_RK5);
	
	PARTITION(Polar,Cartesian);
	BASIS(Cartesian);
}

void SR::SetSrcRoutines(Source_t **LinearSrc, Source_t **NonlinearSrc,
						Source_t **ConstantSrc)
{
	*LinearSrc=StandardLinearity;
	*NonlinearSrc=PrimitiveNonlinearitySR;
	*ConstantSrc=ConstantForcing;
}

void Convolution::SetSrcRoutines(Source_t **LinearSrc, Source_t **NonlinearSrc,
								 Source_t **ConstantSrc)
{
	*LinearSrc=StandardLinearity;
	*NonlinearSrc=PrimitiveNonlinearity;
	*ConstantSrc=ConstantForcing;
}

void PS::SetSrcRoutines(Source_t **LinearSrc, Source_t **NonlinearSrc,
						Source_t **ConstantSrc)
{
	if(!reality) msg(ERROR,"Pseudospectral approximation needs reality=1");
	*LinearSrc=StandardLinearity;
	*NonlinearSrc=PrimitiveNonlinearityFFT;
	*ConstantSrc=ConstantForcing;
}

NWave NWaveProblem;

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
Real continuum_factor=1.0;
static char avgylabel[Nmoment][20];
static char tempbuffer[30];
static int tcount=0;

void NWave::InitialConditions()
{
	int i,n;
	
	Geometry=GeometryProblem->NewGeometry(geometry);
	if(!Geometry->ValidApproximation(Approximation->Abbrev()))
		msg(ERROR,"Geometry \"%s\" is incompatible with Approximation \"%s\"",
			Geometry->Name(),Approximation->Name());
	nyconserve=Npsi=Geometry->Create();
	Ntotal=Geometry->TotalNumber();
	NpsiR=Ntotal-Npsi;
	ny=(average ? Nmoment+1 : 1)*Npsi;
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
		ftin.open(Problem->FileName(dirsep,"t"));
		while(ftin >> t0, ftin.good()) tcount++;
		ftin.close();
	}
	
	open_output(fparam,dirsep,"param",0);
	open_output(fevt,dirsep,"evt");
	open_output(fyvt,dirsep,"yvt");
	open_output(ft,dirsep,"t");
	open_output(fprolog,dirsep,"prolog");
	
	for(n=0; n < Nmoment; n++) 
	
	if(average) for(n=0; n < Nmoment; n++) {
		sprintf(tempbuffer,"avgy%d",n+2);
		mkdir(Problem->FileName(dirsep,tempbuffer),0xFFFF);
		errno=0;
		sprintf(avgylabel[n],"y%%s^%d",n+2);
	}
	
	Problem->GraphicsDump(fparam);
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
	
	out_function(fprolog,K,"K",Npsi);
	out_function(fprolog,Th,"Th",Npsi);
	out_function(fprolog,Area,"Area",Npsi);
	out_function(fprolog,nu_re,"nu.re",Npsi);
	out_function(fprolog,nu_im,"nu.im",Npsi);
	out_function(fprolog,equilibrium,"equil",Npsi);
	out_function(fprolog,Normalization,"normalization",Npsi);

	fevt << "#   t\t\t E\t\t Z\t\t P" << endl;

	// Initialize time integrals to zero.
	if(average) for(i=Npsi; i < ny; i++) y[i]=0.0;
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
	
	out_real(fyvt,y,"y%s",Npsi);
	fyvt.flush();
	
	if(average) for(n=0; n < Nmoment; n++) {
		sprintf(tempbuffer,"avgy%d%st%d",n+2,dirsep,tcount);
		open_output(favgy,dirsep,tempbuffer,0);
		out_curve(favgy,t,"t");
		out_real(favgy,y+Npsi*(n+1),avgylabel[n],Npsi);
		favgy.close();
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
