// Sample file that illustrates the use of Triad

#include "options.h"
#include "kernel.h"

using namespace Array;

const char *method="Ode";
const char *integrator="PC";

typedef array1<Nu>::opt Nuvector;

// Global variables
Var A;
Var B;
Nuvector nu;

// Local vocabulary declarations and default values
static Nu nu0=1.0;

class OdeVocabulary : public VocabularyBase {
public:
  const char *Name() {return "Ode";}
  const char *Abbrev() {return "Ode";}
  OdeVocabulary();
};

class Ode : public ProblemBase {
public:
  Ode();
	
  void InitialConditions();
  void Initialize();
  void Output(int it);
	
  void NonLinearSource(const vector2& Src, const vector2& Y, double t);
  void LinearSource(const vector2& Src, const vector2& Y, double t);
  
  void Source(const vector2& Src, const vector2& Y, double t) {
    NonLinearSource(Src,Y,t);
    LinearSource(Src,Y,t);
  }
};

Ode *OdeProblem;

#include "Ode.h"

Ode::Ode()
{
  OdeProblem=this;
}

OdeVocabulary::OdeVocabulary()
{
  Vocabulary=this;
	
  VOCAB(nu0,(Var) 0.0,(Var) 0.0,"linear damping");
  VOCAB(A,(Var) 0.0,(Var) 0.0,"coefficient A");
  VOCAB(B,(Var) 0.0,(Var) 0.0,"coefficient B");
	
  METHOD(Ode);
	
  // Specialized integrators:
  
  INTEGRATOR(I_Euler);
  INTEGRATOR(E_Euler);
  INTEGRATOR(RB1);
  INTEGRATOR(I_PC);
  INTEGRATOR(E_PC);
  INTEGRATOR(LE_PC);
  INTEGRATOR(C_PC);
  INTEGRATOR(Implicit);
}

OdeVocabulary Ode_Vocabulary;

ofstream fout;

void Ode::InitialConditions()
{
  NY[0]=1;
  Allocate();
  
  Allocate1(nu,ny);
	
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

