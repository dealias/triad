#define COMPLEX 0
#include "kernel.h"
#include "Param.h"
#include "Integrator.h"
#include "Approx.h"

char *approximation="none";
char *integrator="C_PC";

// Vocabulary
double r0=1.0,v0=0.0,th0=0.0,L=1.0,K=1.5,m=1.0;

// Variables;
double A0,E0,beta;

enum {R, V, Th};

void KeplerSource(Var *, Var *, double);

class Kepler : public ProblemBase {
public:
	Kepler();
	char *Name() {return "Kepler";}
	char *Abbrev() {return "kepler";}
	void InitialConditions();
	void Output(int it, int final);
};

class C_PC : public PC {
public:
	char *Name() {return "Conservative Predictor-Corrector";}
	Corrector_t Corrector;
};

Kepler::Kepler()
{
	Problem=this;
	ny=3;
	y=new Var[ny];
	
	VOCAB(r0,0.0,0.0);
	VOCAB(v0,0.0,0.0);
	VOCAB(K,0.0,0.0);
	VOCAB(L,0.0,0.0);
	VOCAB(m,0.0,0.0);
	
	NonlinearSrc=KeplerSource;
	
	APPROXIMATION(None);
	INTEGRATOR(C_PC);
}

Kepler KeplerProblem;

ofstream fout;

void Kepler::InitialConditions()
{
	y[R]=r0;
	y[V]=v0;
	y[Th]=th0;
	
	A0=L*L/(m*r0)-K;	// Assumes th0=0.
	beta=K/A0;
	
	dynamic=0;
	open_output(fout,".",downcase(undashify(Integrator->Abbrev())));
}

void Kepler::Output(int it, int final)
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

void KeplerSource(Var *source, Var *y, double t)
{
	source[R]=y[V];
	source[V]=(L*L/(m*y[R])-K)/(m*y[R]*y[R]);
	source[Th]=L/(m*y[R]*y[R]);
}

inline int C_PC::Corrector(Var *y0, double dt, double& errmax, int start,
						   int stop)
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
