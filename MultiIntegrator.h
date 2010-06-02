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
  unsigned grid;
  // Grid class goes here?
};

class MultiIntegratorBase : public IntegratorBase {
 protected:
  Array1<RK *> Integrator;
  Array1<DynVector<unsigned> > nY;
  MultiProblem *MProblem;
  vector3 Y;
  vector2 YsaveUK; // for recovery from failed time steps
  bool Ysaved;
  unsigned grid;
 public:
  virtual void Allocator(ProblemBase& problem,size_t);

  // These should be provided by the problem class:
  virtual void Project(unsigned)=0;
  virtual void Prolong(unsigned)=0;
  virtual void Rescale()=0; // necessary?

  void Grid(unsigned int g) {
    grid=g;
    //GoyProblem->grid=g;
    MProblem->grid=g;
  }

};

void MultiIntegratorBase::Allocator(ProblemBase& problem,size_t)
{
  if(Ngrids < 2) msg(ERROR,"Need more grids");
  Allocate(Integrator,Ngrids);
  SetProblem(problem);

  Ysaved=0;
  Allocate(Y,Ngrids);
  Allocate(YsaveUK,Ngrids);
  
  Allocate(nY,Ngrids);

  // Assumes that sub-integrators are all the same order
  order=Integrator[0]->Order();
  pgrow=(order > 0) ? 0.5/order : 0;
  pshrink=(order > 1) ? 0.5/(order-1) : pgrow;

  for (unsigned int i=0; i< Ngrids; i++) {
    Grid(i);
    RK *integrator;
    integrator=
      dynamic_cast<RK *>(Vocabulary->NewIntegrator(subintegrator));
    if(!integrator) msg(ERROR,"subintegrator must be an RK integrator");
    
    Integrator[i]=integrator;
    Integrator[i]->SetProblem(problem);
    Integrator[i]->SetParam(*this);
  }
  
  // this should also give an option for rescaling.
  // I guess adding rescaling options to MultiIntegratorBase vocab or something?
  //for (unsigned int i=1; i< Ngrids; i++) Project(i);




}


class GridBase {
  //grid stuff, in the house!
};

#endif
