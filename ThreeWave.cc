#include "options.h"
#include "NWave.h"

class ThreeWaveVocabulary : public VocabularyBase {
public:
	char *Name() {return "Three-Wave";}
	char *Abbrev() {return "w3";}
	ThreeWaveVocabulary();
	Table<LinearityBase> *LinearityTable;
	LinearityBase *NewLinearity(char *& key) {
		return LinearityTable->Locate(key);
	}
};

char *method="SR";
char *integrator="PC";

LinearityBase *Linearity=new LinearityBase[1];

// Global variables
int randomIC=0;
Real tauforce=0.0;

// Local variables
static Mc Mk[]={1.0,1.0,-2.0};
static Nu nu0[]={0.0,0.0,0.0};
static Var IC[]={sqrt(1.5),0.0,sqrt(1.5)};
static Real *K;

ThreeWaveVocabulary::ThreeWaveVocabulary()
{
	Vocabulary=this;
	
	VOCAB(randomIC,0,1);
	VOCAB(Nmoment,0,INT_MAX);
	
	VOCAB_ARRAY(Mk);
	VOCAB_ARRAY(IC);
	VOCAB_ARRAY(nu0);
	
	METHOD(SR);
	
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
}

ThreeWaveVocabulary ThreeWave;

ofstream fout;

void NWave::InitialConditions()
{
	int k,p;

	Nmode=Ntotal=3;
	ny=Nmode*(1+Nmoment);
	NmodeR=Ntotal-Nmode;
	psibuffer=psibufferR=new Var[Nmode];
	psibufferStop=psibuffer+Nmode;
	y=new Var[ny];
	K=new Real [Nmode];
	forcing=new Real[Nmode];
		
	K[0]=sqrt(3.0);
	K[1]=sqrt(9.0);
	K[2]=sqrt(6.0);
	
	nu=nu0;

	pqbuffer=new Var[Nmode*(Nmode+1)/2];
	pqIndex=new Var*[Nmode];
	qStart=new int[Nmode];
	triadLimits=new TriadLimits[Nmode];
	
	Var *pq=pqbuffer;
	for(p=0; p < Nmode; p++) {
		pqIndex[p]=pq-p;
		pq += Nmode-p;
	}
	for(p=0; p < Nmode; p++) qStart[p]=p;
	
	triad[0].Store(pqbuffer+4,Mk[0]);
	triad[1].Store(pqbuffer+2,Mk[1]);
	triad[2].Store(pqbuffer+1,Mk[2]);

	triadLimits[0].start=triad.Base();
	for(k=0; k < Nmode-1; k++) {
		triadLimits[k+1].start=triadLimits[k].stop=triad.Base()+k+1;
	}
	triadLimits[Nmode-1].stop=triad.Base()+Nmode;
	
	for(k=0; k < Nmode; k++) {
		y[k]=IC[k];
		forcing[k]=0.0;
	}
	if(randomIC) {
		for(int i=0; i < Nmode; i++) {
			Var w;
			crand_gauss(&w);
			y[i] *= w;
		}
	}
	
	open_output(fout,dirsep,downcase(undashify(Integrator->Abbrev())));
}

void NWave::Initialize()
{
	int i;

	// Initialize time integrals to zero.
	for(i=Nmode; i < ny; i++) y[i]=0.0;
}

void NWave::ComputeInvariants(Var *y, int Nmode, Real& E, Real& Z, Real& P)
{
	Real Ek,Zk,Pk,k2;
	E=Z=P=0.0;
	for(int i=0; i < Nmode; i++) {
		k2=K[i]*K[i];
		Ek=abs2(y[i]);
		Zk=k2*Ek;
		Pk=k2*Zk;
		E += Ek;
		Z += Zk;
		P += Pk;
	}
	
	Real factor=0.5;
	E *= factor;
	Z *= factor;
	P *= factor;
}	

void NWave::Output(int)
{
	int i;
	
	fout << t << "\t";
	for(i=0; i < Nmode-1; i++) fout << y[i] << "\t";
	fout << y[i] << endl;
}

void ForcingAt(int, Real &)
{
}
