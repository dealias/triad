#include "options.h"
#include "kernel.h"

char *problem="Kepler";
char *integrator="C_PC";

// Vocabulary
double r0=1.0;
double v0=0.0;
double th0=0.0;
double L=1.0;
double K=1.5;
double m=1.0;

// Local variables;
double A0,E0,beta;

enum {R, V, Th};

class KeplerVocabulary : public VocabularyBase {
public:
	char *Name() {return "Kepler";}
	char *Abbrev() {return "kepler";}
	KeplerVocabulary();
};

class Kepler : public ProblemBase {
public:
	void InitialConditions();
	void Output(int it);
	Source_t NonLinearSrc;
};

class C_PC : public PC {
public:
	char *Name() {return "Conservative Predictor-Corrector";}
	int Corrector(double, double&, int, int);
};

KeplerVocabulary::KeplerVocabulary()
{
	Vocabulary=this;
	problem=Name();
	
	VOCAB(r0,0.0,0.0);
	VOCAB(v0,0.0,0.0);
	VOCAB(K,0.0,0.0);
	VOCAB(L,0.0,0.0);
	VOCAB(m,0.0,0.0);
	
	PROBLEM(Kepler);
	INTEGRATOR(C_PC);
}

KeplerVocabulary Kepler_Vocabulary;

ofstream fout;

void Kepler::InitialConditions()
{
	ny=3;
	Nmoment=0;
	y=new Var[ny];
	
	y[R]=r0;
	y[V]=v0;
	y[Th]=th0;
	
	A0=L*L/(m*r0)-K;	// Assumes th0=0.
	beta=K/A0;
	
	dynamic=0;
	open_output(fout,dirsep,downcase(undashify(Integrator->Abbrev())));
}

void Kepler::Output(int it)
{
	int i;
	double arg,A,E;
	
	arg=y[V]*cos(y[Th])-L*sin(y[Th])/(m*y[R]);
	if(arg != 0.0) A=-K*y[V]/arg;
	else A=A0;
	
	E=0.5*m*y[V]*y[V]+0.5*L*L/(m*y[R]*y[R])-K/y[R];
	if(it == 0) E0=E;
	
	fout << t << "\t";
	for(i=0; i < ny-1; i++) fout << y[i] << "\t";
	fout << y[i] << "\t" << E-E0 << "\t" << A-A0 << "\t" << endl; 
	
}

void Kepler::NonLinearSrc(Var *source, Var *y, double)
{
	source[R]=y[V];
	source[V]=(L*L/(m*y[R])-K)/(m*y[R]*y[R]);
	source[Th]=L/(m*y[R]*y[R]);
}

int C_PC::Corrector(double dt, double&, int, int) 
{
	Real DE,arg,lambda,f,diff,lastdiff,th;

	DE=0.5*dt*K*(y0[V]/(y0[R]*y0[R])+y[V]/(y[R]*y[R]));
	y[R]=-K/(-K/y0[R] + DE);
	
	arg=y0[V]*y0[V]+L*L/(m*m*y0[R]*y0[R])-2.0/m*DE-L*L/(m*m*y[R]*y[R]);
	if(arg >= 0) y[V]=sgn(y[V])*sqrt(arg);
	else return 0;

	lambda=L/(m*y[R]);
	
	int i=0;
	diff=DBL_MAX;
	do {
		th=y[Th];
		f=y[V]*(cos(y[Th])+beta)-lambda*sin(y[Th]);
		y[Th] -= f/(-y[V]*sin(y[Th])-lambda*cos(y[Th]));
		lastdiff=diff;
		diff=fabs(y[Th]-th);
		if(++i == 100) msg(ERROR,"Theta iteration did not converge");
	} while(diff != 0.0 && diff < lastdiff);
	return 1;
}
