#include "options.h"
#include "kernel.h"

inline void IntegratorBase::ChangeTimestep(double dtnew)
{
  // New time step must be <= sample.
  if(sample > 0.0 && dtnew > sample) dtnew=sample;
	
  if(abs(dt-dtnew) <= tprecision*abs(dt)) return; // Don't adjust time step.
	
  if(verbose > 1) cout << newl << "Time step changed from " << dt <<
		    " to " << dtnew << " at t=" << t << "." << endl;
  if(dtnew == 0.0) msg(ERROR,"Zero time step encountered");
  dt=dtnew;
  TimestepDependence();
}

static clock_t realtime,lasttime=0;
static const int nperline=10;

void IntegratorBase::Integrate(double& t0, double tmax, 
			       double& dt0, double sample0,
			       int& iteration, unsigned long& nout)
  // Don't dump or microprocess if sample is negative.
{
  double dtold=0.0, dtorig=0.0;
  int it,itx,cont;
  int final=1;
  
  t=t0;
  sample=sample0;
  
  const int forwards=(tmax >= t);
  const double sign=(forwards ? 1.0 : -1.0);
	
  double tstop=((sample > 0.0) ? 0.0 : tmax);
  dt=min(dt0/Microfactor(),dtmax);
  dt *= sign;
  if(dt == 0.0) msg(ERROR,"Zero time step encountered");
	
  TimestepDependence();
  microprocess=(sample >= 0.0) ? Problem->Microprocess() : 0;
	
  // Main integration loop
  for(it=0; it < itmax && (forwards ? t < tmax : t > tmax); it++) {
		
    if(verbose) {
      if((it % nperline) == 0) cout << newl;
      cout << "[" << it << flush;
    }
		
    unsigned long count=nout;
    if(sample == 0.0) dump(t0=t,it,0,tmax);
    else if(sample > 0) {
      if(abs(tstop-t) <= tprecision*abs(tstop)) tstop=t;
      if (forwards ? t >= tstop : t <= tstop) {
        count++;
	tstop=sign*count*sample;
	if((forwards ? tstop > tmax : tstop < tmax) ||
	   abs(tmax-tstop) <= tprecision*abs(tmax)) tmax=tstop=t;
	else dump(t0=t,it,0,tmax);
	if(dtorig) {ChangeTimestep(dtorig); dtorig=0.0;}
      }
    }
		
    statistics(t,dt,it);
    nout=count;
		
    for(itx=0; itx < microsteps; itx++) {
      if(microprocess) Problem->Microprocess();
      if(polltime) {
	realtime=time(NULL);
	if(realtime-lasttime > polltime) {
	  poll();
	  lasttime=realtime;
	}
      }
			
      if(forwards ? t+dt > tstop : t+dt < tstop) {
	if(abs(tstop-t) <= tprecision*abs(tstop)) t=tstop;
	if(t >= tstop) break;
	dtorig=dt;
	ChangeTimestep(tstop-t);
	itx=microsteps-1; // This is the final iteration.
      }
      cont=1;
      do {
	switch(Solve()) {
	case ADJUST:
	  t += dt;
	  ChangeTimestep(sign*max(min(sign*dt*growfactor,dtmax),dtmin));
	  cont=0;
	  break;
	case SUCCESSFUL:
	  t += dt;
	  if(dtold) {ChangeTimestep(dtold); dtold=0.0;}  
	  cont=0;
	  break;
	case NONINVERTIBLE:
	  if(!dtold) dtold=dtorig ? dtorig : dt;
	  invert_cnt++;
	  ChangeTimestep(dt*stepnoninverse);
	  break;
	case UNSUCCESSFUL:
	  if(sign*dt <= dtmin) 
	    msg(ERROR,"Minimum timestep restriction violated");
	  ChangeTimestep(sign*max(sign*dt*shrinkfactor,dtmin));
	}
      } while(cont);
      iteration++;
    }
    
    if(Yout[0] != Y[0]) set(Yout[0],Y[0],ny);
    
    if(verbose) cout << "] ";
  }
	
  if(verbose) cout << endl;
  if(dtorig) ChangeTimestep(dtorig);
  t0=t;
  if(sample >= 0.0) dump(t0,it,final,tmax);
  dt0=dt*sign;
  statistics(t0,dt0,it);
}	

void IntegratorBase::SetProblem(ProblemBase& problem)
{
  Problem=&problem;
}

void IntegratorBase::Alloc0(vector2& Y,
			    const vector& y)
{
  DynVector<unsigned int> *NY=Problem->Index();
  unsigned int nfields=NY->Size();
  Allocate1(Y,nfields);
  Var *p=y;
  for(unsigned int i=0; i < nfields; i++) {
    unsigned int n=(*NY)[i];
    Dimension1(Y[i],n,p);
    p += n;
  }
}

void IntegratorBase::Alloc(vector2& Y,
			   vector& y)
{
  Allocate1(y,ny);
  Alloc0(Y,y);
}

void IntegratorBase::Alloc(vector2& Y)
{
  Var *y;
  Allocate1(y,ny);
  Alloc0(Y,y);
}

void IntegratorBase::Allocate()
{
  check_compatibility(DEBUG);
  
  DynVector<unsigned int> *NY=Problem->Index();
  unsigned int nfields=NY->Size();
  ny=0;
  for(unsigned int i=0; i < nfields; i++) ny += (*NY)[i];
  Alloc(Src,source);
  
  Dimension1(y,Problem->Vector());
  Dimension1(Y,Problem->Vector2());
  Dimension1(Yout,Problem->Vector2());
}

Solve_RC Euler::Solve()
{
  Source(Src,Y,t);
  Problem->Transform(Y,t,dt,YI);
  for(unsigned int j=0; j < ny; j++) y[j] += dt*source[j];
  Problem->BackTransform(Y,t+dt,dt,YI);
  Problem->Stochastic(Y,t,dt);
  return SUCCESSFUL;
}

Solve_RC SYM1::Solve()
{
  Source(Src,Y,t);
  Problem->Transform(Y,t,dt,YI);
  for(unsigned int j=0; j < ny; j += 2) y[j] += dt*source[j];
  Problem->BackTransform(Y,t+dt,dt,YI);
  Source(Src,Y,t);
  Problem->Transform(Y,t,dt,YI);
  for(unsigned int j=1; j < ny; j += 2) y[j] += dt*source[j];
  Problem->BackTransform(Y,t+dt,dt,YI);
  Problem->Stochastic(Y,t,dt);
  return SUCCESSFUL;
}

Solve_RC PC::Solve()
{
  Solve_RC flag;
  
  errmax=0.0;
  errmask=Problem->ErrorMask();
	
  if(new_y0) {
    swaparray(Y0,Y);
    Set1(y,Y[0]);
    Set1(y0,Y0[0]);
    Source(Src0,Y0,t);
  }
  Problem->Transform(Y0,t,dt,YI);
	
  Predictor();
  
  if(Corrector()) {
    flag=(dynamic ? CheckError() : SUCCESSFUL);
    new_y0=(flag != UNSUCCESSFUL);
  } else {
    flag=NONINVERTIBLE;
    new_y0=0;
  }
  
  if(new_y0) {
    Problem->BackTransform(Y,t+dt,dt,YI);
    Problem->Stochastic(Y,t,dt);
  } else if(YI.Size()) {
    swaparray(Y0,YI);
    Set1(y0,Y0[0]);
  }
  
  return flag;
}

void PC::Predictor()
{
  for(unsigned int j=0; j < ny; j++) y[j]=y0[j]+dt*source0[j];
  Problem->BackTransform(Y,t+dt,dt,YI);
}

int PC::Corrector()
{
  Source(Src,Y,t+dt);
  if(dynamic) {
    for(unsigned int j=0; j < ny; j++) {
      y[j]=y0[j]+halfdt*(source0[j]+source[j]);
      if(!errmask || errmask[j]) CalcError(y0[j],y[j],
					   y0[j]+dt*source0[j],y[j]);
    }
    ExtrapolateTimestep();
  } else for(unsigned int j=0; j < ny; j++) {
    y[j]=y0[j]+halfdt*(source0[j]+source[j]);
  }
	
  return 1;
}

void SYM2::Predictor()
{
  for(unsigned int j=0; j < ny; j++) {
    y[j]=y0[j]+halfdt*source0[j];
    if(++j < ny) y[j]=y0[j];
  }
  Problem->BackTransform(Y,t+dt,dt,YI);
  Source(Src,Y,t+halfdt);
  for(unsigned int j=1; j < ny; j += 2) {
    y[j] += dt*source[j];
  }
  Problem->BackTransform(Y,t+dt,dt,YI);
}

int SYM2::Corrector()
{
  Source(Src,Y,t+dt);
  if(dynamic) {
    for(unsigned int j=0; j < ny; j += 2) {
      y[j] += halfdt*source[j];
      if(!errmask || errmask[j]) CalcError(y0[j],y[j],
					   y0[j]+dt*source0[j],y[j]);
    }
    ExtrapolateTimestep();
  } else for(unsigned int j=0; j < ny; j += 2) {
    y[j] += halfdt*source[j];
  }
	
  return 1;
}

// *** Need to add Transform and BackTransform here...

Solve_RC AdamsBashforth::Solve()
{
  swaparray(Y0,Y);
  vector y0=Y0[0];
  vector y=Y[0];
  
  switch(init) {
  case 0:
    {
    leftshiftarray(Src,Src1,Src0);
    vector source0=Src0[0];
    vector source1=Src1[0];
    vector source=Src[0];
    Source(Src0,Y0,t);
    for(unsigned int j=0; j < ny; j++) 
      y[j]=y0[j]+a0*source0[j]+a1*source1[j]+a2*source[j];
    Source(Src,Y,t);
    for(unsigned int j=0; j < ny; j++) 
      y[j]=y0[j]+b0*source[j]+b1*source0[j]+b2*source1[j];
    break;
    }
  case 1:
    swaparray(Src1,Src);
  case 2:
    {
    vector source0=Src0[0];
    vector source=Src[0];
    // Initialize with 2nd-order predictor-corrector
    Source(Src0,Y0,t);
    for(unsigned int j=0; j < ny; j++) y[j]=y0[j]+dt*source0[j];
    Source(Src,Y,t);
    double halfdt=0.5*dt;
    for(unsigned int j=0; j < ny; j++) 
      y[j]=y0[j]+halfdt*(source0[j]+source[j]);
    init--;
    break;
    }
  }
  
  Problem->Stochastic(Y,t,dt);
  return SUCCESSFUL;
}

// *** Need to add Transform and BackTransform and here...
Solve_RC Midpoint::Solve()
{
  int niterations=10;
  swaparray(Y0,Y);
  vector y0=Y0[0];
  vector y=Y[0];
  
  if(verbose > 1) cout << endl;
  Source(Src,Y0,t);
  
  for(int i=0; i < niterations; i++) {
    for(unsigned int j=0; j < ny; j++) y[j]=y0[j]+halfdt*source[j];
    if(verbose > 1) cout << y[0] << endl;
    Source(Src,Y,t);
  }
  for(unsigned int j=0; j < ny; j++) y[j]=y0[j]+dt*source[j];
  if(verbose > 1) cout << y[0] << endl;
	
  Problem->Stochastic(Y,t,dt);
  return SUCCESSFUL;
}

void LeapFrog::Predictor()
{
  if(new_y0) {oldhalfdt=lasthalfdt;}
  else yp=yp0;
  double dtprime=halfdt+oldhalfdt;
  Problem->Transform(YP,t-oldhalfdt,dtprime,YP0);
  for(unsigned int j=0; j < ny; j++) yp[j] += dtprime*source0[j];
  Problem->BackTransform(YP,t+halfdt,dtprime,YP0);
  lasthalfdt=halfdt;
}
	
int LeapFrog::Corrector()
{
  Source(Src,YP,t+halfdt);
  if(dynamic) {
    for(unsigned int j=0; j < ny; j++) {
      y[j]=y0[j]+dt*source[j];
      if(!errmask || errmask[j]) 
	CalcError(y0[j],y[j],y0[j]+dt*source0[j],y[j]);
    }
    ExtrapolateTimestep();
  } else for(unsigned int j=0; j < ny; j++) y[j]=y0[j]+dt*source[j];

  return 1;
}

void RK2::Predictor()
{
  for(unsigned int j=0; j < ny; j++) y[j]=y0[j]+halfdt*source0[j];
  Problem->BackTransform(Y,t+halfdt,halfdt,YI);
}

int RK2::Corrector()
{
  Source(Src,Y,t+halfdt);
  if(dynamic) {
    for(unsigned int j=0; j < ny; j++) {
      y[j]=y0[j]+dt*source[j];
      if(!errmask || errmask[j])
	CalcError(y0[j],y[j],y0[j]+dt*source0[j],y[j]);
    }
    ExtrapolateTimestep();
  } else for(unsigned int j=0; j < ny; j++) y[j]=y0[j]+dt*source[j];

  return 1;
}

void RK4::TimestepDependence()
{
  halfdt=0.5*dt;
  sixthdt=dt/6.0;
}

void RK4::Predictor()
{
  for(unsigned int j=0; j < ny; j++) y[j]=y0[j]+halfdt*source0[j];
  Problem->BackTransform(Y,t+halfdt,halfdt,YI);
  Source(Src1,Y,t+halfdt);
  if(YI.Size()) {swaparray(YI,Y); Set1(y,Y[0]);}
  for(unsigned int j=0; j < ny; j++) y[j]=y0[j]+halfdt*source1[j];
  Problem->BackTransform(Y,t+halfdt,halfdt,YI);
  Source(Src2,Y,t+halfdt);
  if(YI.Size()) {swaparray(YI,Y); Set1(y,Y[0]);}
  for(unsigned int j=0; j < ny; j++) y[j]=y0[j]+dt*source2[j];
  Problem->BackTransform(Y,t+dt,dt,YI);
}

int RK4::Corrector()
{
  Source(Src,Y,t+dt);
  if(YI.Size()) {swaparray(YI,Y); Set1(y,Y[0]);}
  if(dynamic) {
    for(unsigned int j=0; j < ny; j++) {
      y[j]=y0[j]+sixthdt*(source0[j]+2.0*(source1[j]+source2[j])+source[j]);
      if(!errmask || errmask[j])
	CalcError(y0[j],y[j],y0[j]+dt*source2[j],y[j]);
    }
    ExtrapolateTimestep();
  } else for(unsigned int j=0; j < ny; j++) 
    y[j]=y0[j]+sixthdt*(source0[j]+2.0*(source1[j]+source2[j])+source[j]);
  return 1;
}

void RK5::TimestepDependence()
{
  a1=0.2*dt; a2=0.3*dt; a3=0.6*dt; a4=dt; a5=0.875*dt;
  b10=0.2*dt;
  b20=3.0/40.0*dt; b21=9.0/40.0*dt;
  b30=0.3*dt; b31=-0.9*dt; b32=1.2*dt;
  b40=-11.0/54.0*dt; b41=2.5*dt; b42=-70.0/27.0*dt; b43=35.0/27.0*dt;
  b50=1631.0/55296.0*dt; b51=175.0/512.0*dt; b52=575.0/13824.0*dt;
  b53=44275.0/110592.0*dt; b54=253.0/4096.0*dt;
  c0=37.0/378.0*dt; c2=250.0/621.0*dt; c3=125.0/594.0*dt; c5=512.0/1771.0*dt;
  d0=2825.0/27648.0*dt; d2=18575.0/48384.0*dt; d3=13525.0/55296.0*dt;
  d4=277.0/14336.0*dt; d5=0.25*dt;
}

void RK5::Predictor()
{
  //#pragma ivdep		
  for(unsigned int j=0; j < ny; j++) y[j]=y0[j]+b10*source0[j];
  Problem->BackTransform(Y,t+a1,a1,YI);
  Source(Src,Y,t+a1);
  if(YI.Size()) {swaparray(YI,Y); Set1(y,Y[0]);}
  //#pragma ivdep		
  for(unsigned int j=0; j < ny; j++) y[j]=y0[j]+b20*source0[j]+b21*source[j];
  Problem->BackTransform(Y,t+a2,a2,YI);
  Source(Src2,Y,t+a2);
  if(YI.Size()) {swaparray(YI,Y); Set1(y,Y[0]);}
  //#pragma ivdep		
  for(unsigned int j=0; j < ny; j++) y[j]=y0[j]+b30*source0[j]+b31*source[j]+
				b32*source2[j];
  Problem->BackTransform(Y,t+a3,a3,YI);
  Source(Src3,Y,t+a3);
  if(YI.Size()) {swaparray(YI,Y); Set1(y,Y[0]);}
  //#pragma ivdep		
  for(unsigned int j=0; j < ny; j++) y[j]=y0[j]+b40*source0[j]+b41*source[j]+
				b42*source2[j]+b43*source3[j];
  Problem->BackTransform(Y,t+a4,a4,YI);
  Source(Src4,Y,t+a4);
  if(YI.Size()) {swaparray(YI,Y); Set1(y,Y[0]);}
  //#pragma ivdep		
  for(unsigned int j=0; j < ny; j++) y[j]=y0[j]+b50*source0[j]+b51*source[j]+
				b52*source2[j]+b53*source3[j]+
				b54*source4[j];
  Problem->BackTransform(Y,t+a5,a5,YI);
}

int RK5::Corrector()
{
  Source(Src,Y,t+a5);
  if(YI.Size()) {swaparray(YI,Y); Set1(y,Y[0]);}
  if(dynamic) {
    //#pragma ivdep		
    for(unsigned int j=0; j < ny; j++) {
      y[j]=y0[j]+c0*source0[j]+c2*source2[j]+c3*source3[j]+c5*source[j];
      if(!errmask || errmask[j]) {
	Var pred=y0[j]+d0*source0[j]+d2*source2[j]+
	  d3*source3[j]+d4*source4[j]+d5*source[j];
	CalcError(y0[j],y[j],pred,y[j]);
      }
    }
    ExtrapolateTimestep();
  } else {
    //#pragma ivdep		
    for(unsigned int j=0; j < ny; j++) {
      y[j]=y0[j]+c0*source0[j]+c2*source2[j]+c3*source3[j]+c5*source[j];
    }
  }
  return 1;
}
