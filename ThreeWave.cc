#include "NWave.h"

char *NWave::Name() {return "Three-Wave";}
char *NWave::Abbrev() {return "w3";}

char *approximation="none";
char *integrator="PC";

// Three-wave variables
int randomIC=0;
Mc Mkpq[]={1.0,1.0,-2.0};
static Nu nu0[]={0.0,0.0,0.0};
static Var IC[]={sqrt(1.5),0.0,sqrt(1.5)};
Real Area0[]={1.0,1.0,1.0};
Real *K,*Area=Area0,*y2;
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
	int k;

	nyconserve=Npsi=Ntotal=3;
	NpsiR=Ntotal-Npsi;
	psibuffer=psibufferR=new Var[Npsi];
	psibufferStop=psibuffer+Npsi;
	ny=(average ? Nmoment+1 : 1)*Npsi;
	y=new Var[ny];
	y2=new Real[Npsi];
	K=new Real [Npsi];
	
	K[0]=sqrt(3.0);
	K[1]=sqrt(9.0);
	K[2]=sqrt(6.0);
	
	nu=nu0;

	pqbuffer=new Var[Npsi*(Npsi+1)/2];
	pqIndex=new Var*[Npsi];
	qStart=new int[Npsi];
	triadLimits=new TriadLimits[Npsi];
	
	Var *pq=pqbuffer;
	for(int p=0; p < Npsi; p++) {
			pqIndex[p]=pq-p;
			pq += Npsi-p;
	}
	for(p=0; p < n; p++) qStart[p]=p;
	
	triad[0].Store(pqbuffer+4,Mkpq[0]);
	triad[1].Store(pqbuffer+2,Mkpq[1]);
	triad[2].Store(pqbuffer+1,Mkpq[2]);

	triadLimits[0].start=triad.Base();
	for(k=0; k < Npsi-1; k++) {
		triadLimits[k+1].start=triadLimits[k].stop=triad.Base()+k+1;
	}
	triadLimits[Npsi-1].stop=triad.Base()+Npsi;
	
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
