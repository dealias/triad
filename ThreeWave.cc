#include "NWave.h"

char *NWave::Name() {return "Three-Wave";}
char *NWave::Abbrev() {return "w3";}

char *approximation="none";
char *integrator="PC";

const Nu Zero=0.0;

// Three-wave variables
int randomIC=0;
Mc Mkpq[]={1.0,1.0,-2.0};
static Nu nu0[]={Zero,Zero,Zero};
static Var IC[]={sqrt(1.5),0.0,sqrt(1.5)};
Real K0[]={sqrt(3.0),sqrt(9.0),sqrt(6.0)};
Real Area0[]={1.0,1.0,1.0};
Real *K=K0,*Area=Area0,*y2;
Real continuum_factor=1.0;

NWave::NWave()
{
	Problem=this;
	
	VOCAB(randomIC,0,1);
	
	VOCAB_ARRAY(Mkpq);
	VOCAB_ARRAY(IC);
	VOCAB_ARRAY(nu0);
	
	reality=0;
	LinearSrc=StandardLinearity;
	NonlinearSrc=PrimitiveNonlinearity;
	
	APPROXIMATION(None);
	
	INTEGRATOR(C_Euler);
	INTEGRATOR(C_PC);
	INTEGRATOR(E_PC);
	INTEGRATOR(CE_PC);
	INTEGRATOR(C_RK2);
	INTEGRATOR(E_RK2);
	INTEGRATOR(C_RK4);
	INTEGRATOR(E_RK4);
	INTEGRATOR(C_RK5);
}

NWave ThreeWaveProblem;

ofstream fout;

void NWave::InitialConditions()
{
	int k,j;
	
	nyconserve=Npsi=3; 
	psibuffer=new Var[Npsi];
	ny=(average ? Nmoment+1 : 1)*Npsi;
	y=new Var[ny];
	y2=new Real[Npsi];
	
	nu=nu0;

	Npair=3;
	pair=new Pair[Npair];
	pair[0].Store(&psibuffer[1],&psibuffer[2]);
	pair[1].Store(&psibuffer[2],&psibuffer[0]);
	pair[2].Store(&psibuffer[0],&psibuffer[1]);

	j=0;
	for(k=0; k < Npsi; k++) {
		triad[j++].Store(pair[k].Index(),Mkpq[k]);
		triad[j++].Store(NULL,0.0);
	}
	triadBase=triad.Base();
	
	for(k=0; k < Npsi; k++) y[k]=IC[k];
	if(randomIC) {
		for(int i=0; i < Npsi; i++) {
			Var w;
			crand_gauss(w);
			y[i] *= w;
		}
	}
	
	open_output(fout,".",downcase(undashify(Integrator->Abbrev())));
}

void NWave::Initialize()
{
	int i;

	// Initialize time integrals to zero.
	if(average) for(i=Npsi; i < ny; i++) y[i]=0.0;
}

void NWave::Output(int)
{
	int i;
	
	fout << t << "\t";
	for(i=0; i < Npsi-1; i++) fout << y[i] << "\t";
	fout << y[i] << endl;
}
