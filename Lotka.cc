// Note: COMPLEX must be set to 0 in types.h

#include "kernel.h"
#include "Param.h"
#include "Integrator.h"
#include "Approx.h"

char *approximation="none";
char *integrator="C_PC";

// Vocabulary
double x0=1.0,y0=0.4,mu=1.5;

// Variables;
double E,E0;

enum {X, Y};

void LotkaSource(Var *, Var *, double);

class Lotka : public ProblemBase {
public:
	Lotka();
	char *Name() {return "Lotka-Volterra";}
	char *Abbrev() {return "lotka";}
	void InitialConditions();
	void Output(int it);
};

class C_PC : public PC {
public:
	char *Name() {return "Conservative Predictor-Corrector";}
	char *Abbrev() {return "C-PC";}
	int Corrector(double, double&, int, int);
};

Lotka::Lotka()
{
	Problem=this;
	ny=2;
	Nmoment=0;
	y=new Var[ny];
	
	VOCAB(x0,0.0,0.0);
	VOCAB(y0,0.0,0.0);
	VOCAB(mu,0.0,0.0);
	
	APPROXIMATION(None);
	INTEGRATOR(C_PC);
}

void None::SetSrcRoutines(Source_t **, Source_t **NonlinearSrc, Source_t **)
{
	*NonlinearSrc=LotkaSource;
}

Lotka LotkaProblem;

ofstream fout;

Real Mu[2];

void Lotka::InitialConditions()
{
	y[X]=x0;
	y[Y]=y0;
	
	Mu[0]=1.0;
	Mu[1]=mu;
	
	dynamic=0;
	open_output(fout,dirsep,downcase(undashify(Integrator->Abbrev())));
}

void Lotka::Output(int it)
{
	int i;
	
	E=0.0;
	for(int j=0; j < ny; j++) E += Mu[j]*(y[j]-log(y[j]));
	
	if(it == 0) E0=E;
	
	fout << t << "\t";
	for(i=0; i < ny-1; i++) fout << y[i] << "\t";
	fout << y[i] << "\t" << E-E0 << endl; 
}

void LotkaSource(Var *source, Var *y, double)
{
	source[X]=-mu*y[X]*(1.0-y[Y]);
	source[Y]=y[Y]*(1.0-y[X]);
}

int C_PC::Corrector(double dt, double&, int, int)
{
	Real xi[2],lastdiff;
	int j;
	
	for(j=0; j < ny; j++) {
		if(y0[j] <= 0.0) return 0;
		xi[j]=Mu[j]*(y0[j]-log(y0[j]));
	}
	
	Real DE=0.5*dt*mu*((y0[X]-1.0)*(y0[Y]-1.0)+(y[X]-1.0)*(y[Y]-1.0));

	xi[0] += DE;
	xi[1] -= DE;
	
	if(xi[0] < 1.0 || xi[1] < mu) return 0;
	
	for(j=0; j < ny; j++) {
		int i=0;
		Real diff=DBL_MAX;
		do {
			Real old=y[j];
			if(old <= 0.0) return 0;
			Real f=Mu[j]*(y[j]-log(old))-xi[j];
			y[j] -= f/(Mu[j]-Mu[j]/old);
			lastdiff=diff;
			diff=fabs(y[j]-old);
			if(++i == 100) msg(ERROR,"Iteration did not converge");
		} while(diff != 0.0 && diff < lastdiff);
	}
	
	return 1;
}
