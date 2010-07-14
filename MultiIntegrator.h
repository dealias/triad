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
  unsigned getNfields() {return (unsigned) nfields;};

  virtual void InitialConditions(unsigned Ngrids);
  virtual unsigned getnfields(unsigned g)=0;

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
  array1<unsigned> nsave;
  bool Ysaved;
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

  Solve_RC CheckError() {
    bool adjust=true;
    for (unsigned g=0; g < Ngrids ; g++) {
      Solve_RC rc=Integrator[g]->CheckError();
      if(rc == UNSUCCESSFUL)
	return UNSUCCESSFUL;
      adjust &= (rc == ADJUST);
    }
    return adjust ? ADJUST : SUCCESSFUL;
  }
  Solve_RC Solve();
};

extern MultiProblem *MProblem;

void MultiIntegrator::Allocator(ProblemBase& problem,size_t Align)
{
  align=Align;
  Ngrids=::Ngrids;
  if(Ngrids < 2) msg(ERROR,"Need more grids");
  Allocate(Integrator,Ngrids);
  SetProblem(problem);
  MProblem=::MProblem;
  saveF=MProblem->saveF;

  Ysaved=false;
  // FIXME: if Y gets swapped at any point, we're fucked.
  Dimension(mY,MProblem->mY);
  Allocate(Ysave,Ngrids); // should we save more than just one field?
  Allocate(nsave,Ngrids); 
  Allocate(nY,Ngrids);

  // Assumes that sub-integrators are all the same order
  // TODO: this could just take the lowest order and be fine.
  order=Integrator[0]->Order();

  for (unsigned g=0; g < Ngrids; g++) {
    setGrid(g);
    RK *integrator;
    integrator=
      dynamic_cast<RK *>(Vocabulary->NewIntegrator(subintegrator));
    if(!integrator) msg(ERROR,"subintegrator must be an RK integrator");

    Integrator[g]=integrator;
    Integrator[g]->SetProblem(problem);
    Integrator[g]->SetParam(*this);
    Nfields=MProblem->getNfields();
    mY[g].Allocate(Nfields);
    for (unsigned F=0; F < Nfields; ++F) {
      nY[g][F]=Problem->Size(Nfields*grid+F);
      Dimension(mY[g][F],nY[g][F],Problem->YVector()[Nfields*grid+F]);
    }
    Integrator[g]->Allocator(mY[g],&nY[g],Problem->ErrorMask(),align);
    nsave[g]=nY[g][saveF];
    Allocate(Ysave[g],nsave[g]);
  }
  
  order=Integrator[0]->Order();
  pgrow=(order > 0) ? 0.5/order : 0;
  pshrink=(order > 1) ? 0.5/(order-1) : pgrow;


  // this should also give an option for rescaling.
  // Add rescaling options to MultiIntegrator vocab or something?

  for (unsigned g=1; g< Ngrids; g++) 
    MProblem->Project(g);
}

Solve_RC MultiIntegrator::Solve() {
  Solve_RC rc;
  errmax=0.0;
  
  // initialize integrators
  for (unsigned g=0; g < Ngrids; g++) {
    Integrator[g]->SetTime(t,dt);
    Integrator[g]->initialize0();
    // TODO: put first source calculation here?
  }
  if(dynamic) TimestepDependence();
  unsigned laststage=Integrator[0]->NStages()-1;
  unsigned lastgrid=Ngrids-1;
  
  bool new_y0=1;
  setGrid(0);
  
  // save a copy of mY[toG] for recovering from failed time steps
  if(dynamic) {
    if (!Ysaved) {
      for (unsigned g=0; g <= lastgrid; g++) {
	unsigned stop=nsave[g];
	for(unsigned i=0; i < stop; i++)
	  Ysave[g][i]= mY[g][MProblem->saveF][i];
      }
      Ysaved=1;
    }
  }

  for (unsigned j=0; j <= lastgrid; j++) {
    if (new_y0) {
      setGrid(j);
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
  for (unsigned g=0; g < Ngrids; g++)  {
    double err=Integrator[g]->Errmax();
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
  
  for (unsigned g=0; g < Ngrids; g++) {
    Integrator[g]->setnew_y0(new_y0);
    if(new_y0) {
      setGrid(g);
      MProblem->Stochastic(mY[g],t,dt);
    }
  }
  return rc;
}

#endif
