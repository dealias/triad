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
  vector3 mY; // FIXME: should this be public?

  virtual void InitialConditions(unsigned Ngrids);

  // FIXME: Grid base class here?
  virtual void Project(unsigned toG)=0;
  virtual void Prolong(unsigned toG)=0;
  //  virtual void Rescale()=0;
};

void MultiProblem::InitialConditions(unsigned Ngrids0) 
{
  Ngrids=Ngrids0;
  Allocate(mY,Ngrids);
}

class MultiIntegrator : public IntegratorBase {
 protected:
  Array1<RK *> Integrator;
  Array1<DynVector<unsigned> > nY;
  //Note that ProblemBase already has a variable NY.
  MultiProblem *MProblem;
  unsigned Nfields;
  vector3 mY;
  vector2 YsaveUK; // for recovery from failed time steps
  bool Ysaved;
  unsigned grid;
 public:
  virtual void Allocator(ProblemBase& problem,size_t);

  // These should be provided by the problem class:
  //  virtual void Project(unsigned)=0;
  //  virtual void Prolong(unsigned)=0;
  

  const char *Name() {return "MultiIntegrator";}

  void Grid(unsigned int g) {
    grid=g;
    //GoyProblem->grid=g;
    MProblem->grid=g;
  }

  void TimestepDependence() {
    for (unsigned int i=0; i < Ngrids ; i++) {
      Grid(i);
      Integrator[i]->SetTime(t,dt);
      Integrator[i]->TimestepDependence();
    }
  }

  Solve_RC CheckError() {
    bool adjust=true;
    for (unsigned int i=0; i < Ngrids ; i++) {
      Solve_RC rc=Integrator[i]->CheckError();
      if(rc == UNSUCCESSFUL)
	return UNSUCCESSFUL;
      adjust &= (rc == ADJUST);
    }
    return adjust ? ADJUST : SUCCESSFUL;
  }
};

void MultiIntegrator::Allocator(ProblemBase& problem,size_t)
{
  if(Ngrids < 2) msg(ERROR,"Need more grids");
  Allocate(Integrator,Ngrids);
  SetProblem(problem);
  Nfields=MProblem->nfields;
  
  Ysaved=0;
  // FIXME: Y should poing to MProblem->Y
  // FIXME: if Y gets swapped at any point, we're fucked.
  Dimension(mY,MProblem->mY);
  Allocate(YsaveUK,Ngrids); // FIXME: should this deal with more than just UK?
  
  Allocate(nY,Ngrids);

  // Assumes that sub-integrators are all the same order
  // TODO: this could just take the lowest order and be fine.
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
    mY[i].Allocate(Nfields);
    for (unsigned F=0; F < Nfields; ++F) {
      nY[i][F]=Problem->Size(Nfields*grid+F);
      Dimension(mY[i][F],nY[i][F],Problem->YVector()[Nfields*grid+F]);
    }
    Integrator[i]->Allocator(mY[i],&nY[i],Problem->ErrorMask());
  }
  
  // this should also give an option for rescaling.
  // I guess adding rescaling options to MultiIntegrator vocab or something?

  // FIXME: Project should be moved out of the Integrator
  // into the problem class, since it contains a bunch of grids
  /*
  for (unsigned int i=1; i< Ngrids; i++) 
    MProblem->Project(i);
  */
}

class Grid {
  //grid stuff, in the house!
};

#endif
