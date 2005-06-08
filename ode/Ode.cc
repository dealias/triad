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
  
  class C_PC;
  class C_RK2;
  class C_RK4;
  class C_RK5;
  
  class E_PC;
  class E_RK2;
  class E_RK3;
  class E_RK4;
  
  Ode() {
    Table<IntegratorBase> *t=Ode_Vocabulary.IntegratorTable;
    new entry<C_PC,IntegratorBase,Ode>("C_PC",t,this);
    new entry<C_RK2,IntegratorBase,Ode>("C_RK2",t,this);
    new entry<C_RK4,IntegratorBase,Ode>("C_RK4",t,this);
    new entry<C_RK5,IntegratorBase,Ode>("C_RK5",t,this);
    
    new entry<E_PC,IntegratorBase,Ode>("E_PC",t,this);
    new entry<E_RK2,IntegratorBase,Ode>("E_RK2",t,this);
    new entry<E_RK3,IntegratorBase,Ode>("E_RK3",t,this);
    new entry<E_RK4,IntegratorBase,Ode>("E_RK4",t,this);
  }
  
  void IndexLimits(unsigned int& start, unsigned int& stop,
		   unsigned int& startN, unsigned int& stopN) {
    start=0;
    stop=1;
    startN=0;
    stopN=0;
  }
  
};

class Ode::C_PC : public ::C_PC<Ode> {
public:
  C_PC(Ode *parent) : ::C_PC<Ode>(parent) {}
};

class Ode::C_RK2 : public ::C_RK2<Ode> {
public:
  C_RK2(Ode *parent) : ::C_RK2<Ode>(parent) {}
};

class Ode::C_RK4 : public ::C_RK4<Ode> {
public:
  C_RK4(Ode *parent) : ::C_RK4<Ode>(parent) {}
};

class Ode::C_RK5 : public ::C_RK5<Ode> {
public:
  C_RK5(Ode *parent) : ::C_RK5<Ode>(parent) {}
};

class Ode::E_PC : public ::E_PC<Ode> {
public:
  E_PC(Ode *parent) : ::E_PC<Ode>(parent) {}
};

class Ode::E_RK2 : public ::E_RK2<Ode> {
public:
  E_RK2(Ode *parent) : ::E_RK2<Ode>(parent) {}
};

class Ode::E_RK3 : public ::E_RK3<Ode> {
public:
  E_RK3(Ode *parent) : ::E_RK3<Ode>(parent) {}
};

class Ode::E_RK4 : public ::E_RK4<Ode> {
public:
  E_RK4(Ode *parent) : ::E_RK4<Ode>(parent) {}
};

Ode *OdeProblem;

//#include "Ode.h"

class TestIntegrator : public IntegratorBase {
  IntegratorBase *Integrator, *Integrator1, *Integrator2;
public:  
  const char *Name() {return "Test";}
//  void Unswap() {Integrator->Unswap();}
  void Allocator() {
    IntegratorBase::Allocator();
    Integrator1=new PC;
    Integrator1->Allocator(*Problem);
    Integrator2=new PC;
    Integrator2->Allocator(*Problem);
  }
  Solve_RC Solve(IntegratorBase *Integrator0) {
    Integrator=Integrator0;
    Integrator->SetTime(t,dt);  // copy errmax as well
    Solve_RC rc=Integrator->Solve();
//    Integrator->Unswap();
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

void Ode::NonLinearSource(const vector2& Src, const vector2& Y, double)
{
  vector source=Src[0];
  vector y=Y[0];
//  source[0]=cos(y[0]);
//  source[0]=cos(t)*y[0];
  source[0]=-A*y[0]-B*y[0]*y[0];
}

void Ode::LinearSource(const vector2& Src, const vector2& Y, double)
{
  vector source=Src[0];
  vector y=Y[0];
  source[0] -= nu[0]*y[0];
}
