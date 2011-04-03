#ifndef __MultiIntegrator_h__
#define __MultiIntegrator_h__ 1

#include "kernel.h"
#include "Array.h"

using namespace Array;

extern unsigned Ngrids;
extern const char *subintegrator;

class MultiProblem : public ProblemBase {
 protected:
  DynVector<unsigned> gridindex;
 public:
  enum Field {Nfields};
  unsigned Ngrids, grid, nfields;
  vector3 mY;
  vector y0; // saved data
  virtual unsigned getNfields(unsigned g) {return (unsigned) nfields;};

  void InitialConditions() {}
  virtual void InitialConditions(unsigned Ngrids);
  virtual unsigned getnfields(unsigned g)=0;

  virtual void Project(unsigned toG)=0;
  virtual void Prolong(unsigned toG)=0;
  virtual int Rescale()=0;
};

void MultiProblem::InitialConditions(unsigned Ngrids0) 
{
  Ngrids=Ngrids0;
  Allocate(mY,Ngrids);
}

class MultiIntegrator : public IntegratorBase {
 public:
  vector y0;
 private:
 protected:
  bool new_y0;
  Array1<RK *> Integrator;
  Array1<DynVector<unsigned> > nY; //NB: ProblemBase has a variable NY.
  MultiProblem *MProblem;
  unsigned Nfields;
  vector3 mY;
  unsigned grid;
 public:
  const char *Name() {return "MultiIntegrator";}
  virtual void Allocator(ProblemBase& problem,size_t Align=0);

  void setGrid(unsigned g) {
    grid=g;
    MProblem->grid=g;
  }

  void TimestepDependence() {
    for (unsigned g=0; g < Ngrids ; g++) {
      setGrid(g);
      Integrator[g]->SetTime(t,dt);
      Integrator[g]->TimestepDependence();
    }
  }

  Solve_RC Solve();
};

extern MultiProblem *MProblem;

void MultiIntegrator::Allocator(ProblemBase& problem, size_t Align)
{
  const vector2& YP=problem.YVector();
  DynVector<unsigned int>* NYP=problem.Sizes();
  align=Align;
  
  ny=0;
  NY.SetDynVector(*NYP);
  unsigned int nfields=YP.Size();
  for(unsigned int i=0; i < nfields; i++)
    ny += NY[i];
  
  Dimension(Y,YP);
  Dimension(y,ny,(Var *) (YP[0]));
  Allocate(y0,ny,align);
  
  Ngrids=::Ngrids;
  if(Ngrids < 2) msg(ERROR,"Need more grids");
  
  MProblem=::MProblem;
  SetProblem(problem);
  Allocate(Integrator,Ngrids);
  
  Dimension(mY,MProblem->mY);
  Set(MProblem->y0,y0);
  Dimension(MProblem->y0,ny);
  Allocate(nY,Ngrids);

  for (unsigned g=0; g < Ngrids; g++) {
    setGrid(g);
    RK *integrator=dynamic_cast<RK *>(Vocabulary->NewIntegrator(subintegrator));
    if(!integrator) msg(ERROR,"subintegrator must be an RK integrator");

    Integrator[g]=integrator;
    Integrator[g]->SetProblem(problem);
    Integrator[g]->SetParam(*this);
    Nfields=MProblem->getNfields(g);
    mY[g].Allocate(Nfields);
    for (unsigned F=0; F < Nfields; ++F) {
      nY[g][F]=Problem->Size(Nfields*g+F);
      Dimension(mY[g][F],nY[g][F],Problem->YVector()[Nfields*g+F]);
    }
    Integrator[g]->Allocator(mY[g],&nY[g],Problem->ErrorMask(),align);
  }
  
  new_y0=true;
    
  // Assumes that sub-integrators are all the same order
  // TODO: this could just take the lowest order and be fine.
  order=Integrator[0]->Order();
  pgrow=(order > 0) ? 0.5/order : 0;
  pshrink=(order > 1) ? 0.5/(order-1) : pgrow;


  // this should also give an option for rescaling.
  // Add rescaling options to MultiIntegrator vocab or something?
}

Solve_RC MultiIntegrator::Solve() {
  Solve_RC flag;
  
  errmax=0.0;
  
  if(new_y0)
    for(unsigned i=0; i < ny; ++i) y0[i]=y[i]; // Save the initial data.
  else
    for(unsigned i=0; i < ny; ++i) y[i]=y0[i]; // Reload the initial data.
  
  TimestepDependence();
  unsigned lastgrid=Ngrids-1;
  
  for (unsigned g=0; g <= lastgrid; g++) {
    setGrid(g);
    Integrator[g]->initialize0();
    Integrator[g]->iSource();
    Integrator[g]->Predictor(0,Integrator[g]->Ny());
      
    if(Integrator[g]->Corrector(0,Integrator[g]->Ny())) {
      double err=Integrator[g]->Errmax();
      if(err > errmax) errmax=err;
      flag=dynamic ? CheckError() : SUCCESSFUL;
      new_y0=(flag != UNSUCCESSFUL);
    } else {
      flag=NONINVERTIBLE;
      new_y0=false;
      break;
    }
        
    if(new_y0) {
      if(g < lastgrid)
        MProblem->Project(g+1); 
    }
  }


  if(MProblem->Rescale() > 0) {
    flag=NONINVERTIBLE;
    new_y0=false;
  }

  if (new_y0) {
    for (unsigned g=lastgrid; g > 0; g--)
      MProblem->Prolong(g-1);
    for (unsigned g=0; g < Ngrids; g++) {
      setGrid(g);
      MProblem->Stochastic(mY[g],t,dt);
    }
  }
  
  return flag;
}

#endif
