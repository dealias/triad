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
  unsigned saveF; // index of fields to save

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
 private:
  unsigned saveF;
 protected:
  Array1<RK *> Integrator;
  Array1<DynVector<unsigned> > nY; //NB: ProblemBase has a variable NY.
  MultiProblem *MProblem;
  unsigned Nfields;
  vector3 mY;
  vector2 Ysave; // for recovery from failed time steps
  bool Ysaved;
  unsigned grid;
 public:
  const char *Name() {return "MultiIntegrator";}
  virtual void Allocator(ProblemBase& problem,size_t);

  void Grid(unsigned g) {
    grid=g;
    MProblem->grid=g;
  }

  void TimestepDependence() {
    for (unsigned i=0; i < Ngrids ; i++) {
      Grid(i);
      Integrator[i]->SetTime(t,dt);
      Integrator[i]->TimestepDependence();
    }
  }

  Solve_RC CheckError() {
    bool adjust=true;
    for (unsigned i=0; i < Ngrids ; i++) {
      Solve_RC rc=Integrator[i]->CheckError();
      if(rc == UNSUCCESSFUL)
	return UNSUCCESSFUL;
      adjust &= (rc == ADJUST);
    }
    return adjust ? ADJUST : SUCCESSFUL;
  }
  Solve_RC Solve();
};

void MultiIntegrator::Allocator(ProblemBase& problem,size_t)
{
  if(Ngrids < 2) msg(ERROR,"Need more grids");
  Allocate(Integrator,Ngrids);
  SetProblem(problem);
  cout << MProblem << endl;
  cout << &problem << endl;
  //MProblem=&problem; // FIXME: this doesn't work.
  exit(1);
  
  Nfields=MProblem->nfields;
  saveF=MProblem->saveF;

  Ysaved=0;
  // FIXME: if Y gets swapped at any point, we're fucked.
  Dimension(mY,MProblem->mY);
  Allocate(Ysave,Ngrids); // FIXME: should we save more than just one field?
  
  Allocate(nY,Ngrids);

      

  // Assumes that sub-integrators are all the same order
  // TODO: this could just take the lowest order and be fine.
  order=Integrator[0]->Order();
  pgrow=(order > 0) ? 0.5/order : 0;
  pshrink=(order > 1) ? 0.5/(order-1) : pgrow;

  for (unsigned i=0; i< Ngrids; i++) {
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
    Allocate(Ysave[i],nY[i][saveF]);
  }


  
  // this should also give an option for rescaling.
  // I guess adding rescaling options to MultiIntegrator vocab or something?

  // FIXME: Project should be moved out of the Integrator
  // into the problem class, since it contains a bunch of grids

  for (unsigned i=1; i< Ngrids; i++) 
    MProblem->Project(i);

}

Solve_RC MultiIntegrator::Solve() {
  Solve_RC rc;
  errmax=0.0;
  
  // initialize integrators
  for (unsigned i=0; i < Ngrids; i++) {
    Integrator[i]->SetTime(t,dt);
    Integrator[i]->initialize0();
    // TODO: put first source calculation here?
  }
  if(dynamic) TimestepDependence();
  unsigned laststage=Integrator[0]->NStages()-1;
  unsigned lastgrid=Ngrids-1;
  
  bool new_y0=1;
  Grid(0);
  
  // save a copy of mY[toG] for recovering from failed time steps
  if (!Ysaved) {
    for (unsigned j=0; j <= lastgrid; j++) {
      unsigned stop=nY[j][saveF];
      for(unsigned i=0; i < stop; i++)
	Ysave[j][i]= mY[j][MProblem->saveF][i];
    }
    Ysaved=1;
  }
  for (unsigned j=0; j <= lastgrid; j++) {
    if (new_y0) {
      Grid(j);
      Integrator[j]->iSource();
      for (unsigned i=0; i < laststage; ++i) {
	Integrator[j]->PStage(i);
	Integrator[j]->Source(i);
      }
      
      Integrator[j]->Corrector(0,Integrator[j]->Ny());
      double err=Integrator[j]->Errmax();
      if(err > errmax) errmax=err;
      if (errmax < tolmax2) {
	if(j < lastgrid) {
	  MProblem->Project(j+1); 
	  Integrator[j+1]->initialize0();
	}
      } else 
	new_y0=0;
    }
  }
  for (unsigned i=0; i < Ngrids; i++)  {
    double err=Integrator[i]->Errmax();
    if(err > errmax) errmax=err;
  }
  
  rc=(dynamic ? CheckError() : SUCCESSFUL);
  if (new_y0) 
    new_y0=rc != UNSUCCESSFUL;
  
  if (new_y0) {
    for (unsigned j=lastgrid; j > 0; j--)
      MProblem->Prolong(j-1);
    Ysaved=0;
  } else {
    // revert completed runs
    for (unsigned j=0; j <= lastgrid; j++) {
      unsigned stop=nY[j][saveF];
      for(unsigned i=0; i < stop; i++)
	mY[j][MProblem->saveF][i]=Ysave[j][i]; 
    // FIXME: some grids may not have been projected onto, so they
    // don't need to be reverted
    }
  }
  
  for (unsigned i=0; i < Ngrids; i++) {
    Integrator[i]->setnew_y0(new_y0);
    if(new_y0) {
      Grid(i);
      MProblem->Stochastic(mY[i],t,dt);
    }
  }
  return rc;
}

class Grid {
  //grid stuff, in the house!
};

#endif
