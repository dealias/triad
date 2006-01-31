// Sample file that illustrates the use of Triad

#include "options.h"
#include "kernel.h"
#include "Exponential.h"
#include "Conservative.h"

using namespace Array;

const double ProblemVersion=1.0;

const char *method="Ode";
const char *integrator="PC";

typedef array1<Nu>::opt Nuvector;

// Global variables
Var A;
Var B;

// Local vocabulary declarations and default values
static Nu nu0=1.0;

class OdeVocabulary : public VocabularyBase {
public:
  const char *Name() {return "Ode";}
  const char *Abbrev() {return "Ode";}
  OdeVocabulary();
};

OdeVocabulary Ode_Vocabulary;

class Ode : public ProblemBase {
  Nuvector nu;
public:
  void InitialConditions();
  void Initialize();
  void Output(int it);
  void FinalOutput();
	
  Nu LinearCoeff(unsigned int i) {
    return nu[i];
  }

  void NonLinearSource(const vector2& Src, const vector2& Y, double t);
  void LinearSource(const vector2& Src, const vector2& Y, double t);
  
  void ConservativeSource(const vector2& Src, const vector2& Y, double t) {
    NonLinearSource(Src,Y,t);
    LinearSource(Src,Y,t);
  }
  void NonConservativeSource(const vector2& Src, const vector2& Y, double t) {}
  void ExponentialSource(const vector2& Src, const vector2& Y, double t) {
    NonLinearSource(Src,Y,t);
    NonConservativeSource(Src,Y,t);
  }
  void Source(const vector2& Src, const vector2& Y, double t) {
    ConservativeSource(Src,Y,t);
    NonConservativeSource(Src,Y,t);
  }
  
  Ode() {
    ExponentialIntegrators(Ode_Vocabulary.IntegratorTable,this);
    ConservativeIntegrators(Ode_Vocabulary.IntegratorTable,this);
  }
  
  void IndexLimits(unsigned int& start, unsigned int& stop,
		   unsigned int& startT, unsigned int& stopT,
		   unsigned int& startM, unsigned int& stopM) {
    start=0;
    stop=1;
    startT=0;
    stopT=0;
    startM=0;
    stopM=0;
  }
  
};

Ode *OdeProblem;

//#include "Ode.h"

class TestIntegrator : public IntegratorBase {
  IntegratorBase *Integrator, *Integrator1, *Integrator2;
public:  
  const char *Name() {return "Test";}
  void Unswap() {Integrator->Unswap();}
  void Allocator() {
    IntegratorBase::Allocator();
    Integrator1=new PC;
    Integrator1->Allocator(*Problem);
    Integrator2=new PC;
    Integrator2->Allocator(*Problem);
    Integrator=Integrator1;
  }
  Solve_RC Solve(IntegratorBase *I) {
    I->Sync(Integrator);
    Integrator=I;
    I->SetTime(t,dt);  // copy errmax as well
    Solve_RC rc=I->Solve();
    return rc;
  }
  
  Solve_RC Solve() {
    Solve_RC rc=Solve(Integrator1);
    if(rc == SUCCESSFUL || rc == ADJUST) rc=Solve(Integrator2);
    return rc;
  }
};

OdeVocabulary::OdeVocabulary()
{
  Vocabulary=this;
	
  VOCAB(nu0,(Real) 0.0,(Real) 0.0,"linear damping");
  VOCAB(A,(Var) 0.0,(Var) 0.0,"coefficient A");
  VOCAB(B,(Var) 0.0,(Var) 0.0,"coefficient B");
	
  METHOD(Ode);
	
  // Specialized integrators:
  
  INTEGRATOR(TestIntegrator);
  INTEGRATOR(RK2p);
  INTEGRATOR(RK5p);
#if 0  
  INTEGRATOR(E_Euler);
  INTEGRATOR(I_Euler);
  INTEGRATOR(RB1);
  INTEGRATOR(I_PC);
  
  INTEGRATOR(Implicit);
#endif  
}

ofstream fout;

void Ode::InitialConditions()
{
  OdeProblem=this;
  NY[0]=1;
  Allocator();
  
  Allocate(nu,ny);
	
  vector y=Y[0];
  y[0]=1.0;
  nu[0]=nu0;
	
  open_output(fout,dirsep,downcase(undashify(Integrator->Abbrev())));
}

void Ode::Initialize()
{
  fout << "# " << Integrator->Name() << endl;
}

void Ode::Output(int)
{
  fout << t << "\t" << abs(y[0]) << "\t" << endl;
}

unsigned int count=0;
void Ode::NonLinearSource(const vector2& Src, const vector2& Y, double)
{
  count++;
  vector source=Src[0];
  vector y=Y[0];
//  source[0]=cos(y[0]);
//  source[0]=cos(t)*y[0];
//  source[0]=-A*y[0]-B*y[0]*y[0];
    source[0]=-A*y[0]-B*exp(y[0].re);
}

void Ode::LinearSource(const vector2& Src, const vector2& Y, double)
{
  vector source=Src[0];
  vector y=Y[0];
  source[0] -= nu[0]*y[0];
}

void Ode::FinalOutput() {
  cout << endl << "Number of source evaluations: " << count << endl;
}
