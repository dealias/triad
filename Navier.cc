#include "NWave.h"
#include "Polar.h"

char *NWave::Name() {return "N-Wave";}
char *NWave::Abbrev() {return "nw";}

char *approximation="SR";
char *integrator="PC";
char *geometry="Polar";

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

Real kL=0.0;
Real kH=0.9*DBL_MAX;

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
	
	VOCAB(Nr,1,INT_MAX);
	VOCAB(Nth,1,INT_MAX);
	VOCAB(krmin,0.0,DBL_MAX);
	VOCAB(krmax,0.0,DBL_MAX);
	VOCAB(kthmin,0.0,twopi);
	
	LinearSrc=StandardLinearity;
	NonlinearSrc=PrimitiveNonlinearitySR;
	ConstantSrc=ConstantForcing;
	
	GeometryTable=new Table<GeometryBase>("Geometry",GeometryCompare,
										  GeometryKeyCompare);
	APPROXIMATION(SR);
	
	INTEGRATOR(C_Euler);
	INTEGRATOR(C_PC);
	INTEGRATOR(E_PC);
	INTEGRATOR(CE_PC);
	INTEGRATOR(C_RK2);
	INTEGRATOR(E_RK2);
	INTEGRATOR(C_RK4);
	INTEGRATOR(E_RK4);
	INTEGRATOR(C_RK5);
	
	GEOMETRY(Polar);
}

NWave NWaveProblem;

Real forcek(int i) 
{
	Real k=Geometry->K(i);
	if(abs(k-kforce) < 0.5*deltaf) return force;
	else return 0.0;
}

Real growth(const Polar& v) 
{
	Real k=v.r;
	Real gamma=0.0;
	if(k <= kL) gamma -= pow(k,pL)*nuL*pow(k,pL);
	if(abs(k-kforce) < 0.5*deltaf) gamma += gammaf/deltaf;
	if(k > kH) gamma -= pow(k,pH)*nuH*pow(k,pH);
	return gamma;
}

Real frequency(const Polar&)
{
	return 0.0;
}

void ComputeLinearity(int i, Real& nu)
{
	Real k=Geometry->K(i);
	Real th=Geometry->Th(i);
	nu=-growth(Polar(k,th));
}

void ComputeLinearity(int i, Complex& nu)
{
	Real k=Geometry->K(i);
	Real th=Geometry->Th(i);
	nu=-growth(Polar(k,th))+I*frequency(Polar(k,th));
}

inline Nu Linearity(int i)
{
	Nu nu;
	ComputeLinearity(i,nu);
	return nu;
}

static Real equilibrium(int i)
{
	Real k=Geometry->K(i);
	return 0.5/(alpha+beta*k*k);
}

static ofstream fparam,fevt,fekvt,favgy[Nmoment],fprolog;
Real *K,*Th,*Area,*y2;
Real continuum_factor=1.0/twopi2;
static Real *equil;
static Complex *nuC;

void NWave::InitialConditions()
{
	int i,n;
	
	Geometry=GeometryProblem->NewGeometry(geometry);
	nyconserve=Npsi=Geometry->Create();
	Ntotal=Geometry->TotalNumber();
	NpsiR=Ntotal-Npsi;
	ny=(average ? Nmoment+1 : 1)*Npsi;
	y=new Var[ny];
	y2=new Real[Npsi];
	equil=new Real[Npsi];
	nu=new Nu[Npsi];
	nuC=new Complex[Npsi];
	forcing=new Real[Npsi];
	
	K=new Real[Npsi];
	Th=new Real[Npsi];
	Area=new Real[Npsi];
	
	for(i=0; i < Npsi; i++) {
		nuC[i]=nu[i]=Geometry->BinAveragedLinearity(i);
		equil[i]=equilibrium(i);
		y[i]=sqrt(2.0*equil[i]);
		K[i]=Geometry->K(i);
		Th[i]=Geometry->Th(i);
		Area[i]=Geometry->Area(i);
		forcing[i]=forcek(i);
	}
	
	// Randomize the initial conditions.	
	// If reality condition is not explicitly enforced and the number of
	// angular bins is even, force the initial conditions to respect reality. 
	if(randomIC) {
		n=(reality || Npsi % 2) ? Npsi : Npsi/2;
		for(i=0; i < n; i++) {
			Var w;
			crand_gauss(&w);
			y[i] *= w;
		}
		for(i=n; i < Npsi; i++) y[i]=conj(y[i-n]);
	}
	
	open_output(fparam,dirsep,"param",0);
	open_output(fevt,dirsep,"evt");
	open_output(fekvt,dirsep,"ekvt");
	open_output(fprolog,dirsep,"prolog");
	
	if(average) for(n=0; n < Nmoment; n++) {
		char temp[10];
		sprintf(temp,"avgy%d",n+2);
		open_output(favgy[n],dirsep,temp);
	}
	
	Problem->GraphicsDump(fparam);
	fparam.close();
	
}

void NWave::Initialize()
{
	int i;
	
	out_curve(fprolog,K,"K",Npsi);
	out_curve(fprolog,Th,"Th",Npsi);
	out_curve(fprolog,Area,"Area",Npsi);
	out_curve(fprolog,nuC,"nu",Npsi);
	out_curve(fprolog,equil,"equil",Npsi);

	fevt << "#   t\t\t E\t\t Z\t\t P" << endl;

	// Initialize time integrals to zero.
	if(average) for(i=Npsi; i < ny; i++) y[i]=0.0;
	
	delete [] nuC;
	delete [] equil;
}

void NWave::Output(int)
{
	Real E,Z,P;
	int i,n;
	
	for(i=0; i < Npsi; i++) y2[i] = abs2(y[i]);
	compute_invariants(y2,Npsi,E,Z,P);
	
	fevt << t << "\t" << E << "\t" << Z << "\t" << P << endl << flush;
	
	out_curve(fekvt,t,"t");
	out_curve(fekvt,y,"y",Npsi);
	fekvt.flush();
	
	if(average) for(n=0; n < Nmoment; n++) {
		out_curve(favgy[n],t,"t");
		out_curve(favgy[n],y+Npsi*(n+1),"y^n",Npsi);
		favgy[n].flush();
	}
}
