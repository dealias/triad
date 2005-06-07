#include "options.h"
#include "kernel.h"

using namespace Array;

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
	  ExtrapolateTimestep();
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
	  ExtrapolateTimestep();
	  if(sign*dt <= dtmin) 
	    msg(ERROR,"Minimum timestep restriction violated");
	  ChangeTimestep(sign*max(sign*dt*shrinkfactor,dtmin));
	}
      } while(cont);
      iteration++;
    }
    
    Unswap();
    
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

void IntegratorBase::Alloc(vector2& Y0, vector& y)
{
  Allocate(y,ny);
  unsigned int nfields=Y.Size();
  Allocate(Y0,nfields);
  Var *p=y;
  for(unsigned int i=0; i < nfields; i++) {
    unsigned int n=Problem->Size(i);
    Dimension(Y0[i],n,p);
    p += n;
  }
}

void IntegratorBase::Alloc0(vector2& Y, vector& y)
{
  Alloc(Y,y);
  for(unsigned int i=0; i < ny; i++) y[i]=0.0;
}

void IntegratorBase::Allocator(const vector2& Y0,
			       const ivector& errmask0)
{
  check_compatibility(DEBUG);
  
  unsigned int nfields=Y0.Size();
  ny=0;
  for(unsigned int i=0; i < nfields; i++)
    ny += Y0[i].Size();
  
  Dimension(Y,Y0);
  Alloc0(Src,source);
  
  Dimension(errmask,errmask0);
  Dimension(Yout,Y0);
  Dimension(y,ny,Y[0]());
  
  pgrow=(order > 0) ? 0.5/order : 0;
  pshrink=(order > 1) ? 0.5/(order-1) : pgrow;
}

Solve_RC Euler::Solve()
{
  Source(Src,Y,t);
  Problem->Transform(Y,t,dt,YI);
  Predictor(0,ny);
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
	
  if(new_y0) {
//swaparray(Y0,Y); // TODO: Reinstate this optimization
    for(unsigned i=0; i < Y0.Size(); ++i) {
      Y0[i]=Y[i];
    }
    Set(y,Y[0]);
    Set(y0,Y0[0]);
    Source(Src0,Y0,t);
  }
  Problem->Transform(Y0,t,dt,YI);
	
  Predictor(0,ny);
  
  if(Corrector(0,ny)) {
    flag=(dynamic ? CheckError() : SUCCESSFUL);
    new_y0=(flag != UNSUCCESSFUL);
  } else {
    flag=NONINVERTIBLE;
    new_y0=0;
  }
  
  if(new_y0) {
    Problem->BackTransform(Y,t+dt,dt,YI);
    Problem->Stochastic(Y,t,dt);
  } else if(Active(YI)) { // TODO: Reinstate this optimization
//    swaparray(Y0,YI);
//    Set(y0,Y0[0]);
    for(unsigned i=0; i < Y0.Size(); ++i) {
      YI[i]=Y0[i];
    }
  }
  
  return flag;
}

void PC::Predictor(unsigned int start, unsigned int stop)
{
  for(unsigned int j=start; j < stop; j++) y[j]=y0[j]+dt*source0[j];
  Problem->BackTransform(Y,t+dt,dt,YI);
}

int PC::Corrector(unsigned int start, unsigned int stop)
{
  CSource(Src,Y,t+dt);
  if(dynamic) {
    for(unsigned int j=start; j < stop; j++) {
      y[j]=y0[j]+halfdt*(source0[j]+source[j]);
      if(!Active(errmask) || errmask[j])
	CalcError(y0[j],y[j],y0[j]+dt*source0[j],y[j]);
    }
  } else for(unsigned int j=start; j < stop; j++) {
    y[j]=y0[j]+halfdt*(source0[j]+source[j]);
  }
	
  return 1;
}

void SYM2::Predictor(unsigned int start, unsigned int stop)
{
  for(unsigned int j=start; j < stop; j++) {
    y[j]=y0[j]+halfdt*source0[j];
    if(++j < stop) y[j]=y0[j];
  }
  Problem->BackTransform(Y,t+dt,dt,YI);
  PSource(Src,Y,t+halfdt);
  for(unsigned int j=1; j < stop; j += 2) {
    y[j] += dt*source[j];
  }
  Problem->BackTransform(Y,t+dt,dt,YI);
}

int SYM2::Corrector(unsigned int start, unsigned int stop)
{
  CSource(Src,Y,t+dt);
  if(dynamic) {
    for(unsigned int j=start; j < stop; j += 2) {
      y[j] += halfdt*source[j];
      if(!Active(errmask) || errmask[j])
	CalcError(y0[j],y[j],y0[j]+dt*source0[j],y[j]);
    }
  } else for(unsigned int j=start; j < stop; j += 2) {
    y[j] += halfdt*source[j];
  }
	
  return 1;
}

Solve_RC AdamsBashforth::Solve()
{
  Solve_RC flag=UNSUCCESSFUL;
  errmax=0.0;
  
  swaparray(Y0,Y);
  Set(y,Y[0]);
  Set(y0,Y0[0]);
  
  switch(init) {
  case 0:
    {
    leftshiftarray(Src,Src1,Src0);
    Set(source0,Src0[0]);
    Set(source1,Src1[0]);
    Set(source,Src[0]);
    Source(Src0,Y0,t);
    Problem->Transform(Y0,t,dt,YI);
    for(unsigned int j=0; j < ny; j++) 
      y[j]=y0[j]+a0*source0[j]+a1*source1[j]+a2*source[j];
    Problem->BackTransform(Y,t+dt,dt,YI);
    Source(Src,Y,t);
    if(dynamic) {
      for(unsigned int j=0; j < ny; j++) {
	Var temp=y0[j]+b0*source[j]+b1*source0[j]+b2*source1[j];
	CalcError(y0[j],temp,y[j],temp);
	y[j]=temp;
      }
      flag=CheckError();
    } else {
      for(unsigned int j=0; j < ny; j++) {
	y[j]=y0[j]+b0*source[j]+b1*source0[j]+b2*source1[j];
      }
      flag=SUCCESSFUL;
    }
    if (flag != UNSUCCESSFUL) Problem->BackTransform(Y,t+dt,dt,YI);
    break;
    }
  case 1:
    swaparray(Src1,Src);
  case 2:
    {
    Set(source0,Src0[0]);
    Set(source,Src[0]);
    // Initialize with 2nd-order predictor-corrector **IMPROVE: USE RK4?**
    Source(Src0,Y0,t);
    Problem->Transform(Y0,t,dt,YI);
    for(unsigned int j=0; j < ny; j++) y[j]=y0[j]+dt*source0[j];
    Problem->BackTransform(Y,t+dt,dt,YI);
    Source(Src,Y,t);
    double halfdt=0.5*dt;
    for(unsigned int j=0; j < ny; j++) 
      y[j]=y0[j]+halfdt*(source0[j]+source[j]);
    Problem->BackTransform(Y,t+dt,dt,YI);
    flag=SUCCESSFUL;
    init--;
    break;
    }
  }
  
  Problem->Stochastic(Y,t,dt);
  return flag;
}

Solve_RC Midpoint::Solve()
{
  int niterations=10;
  swaparray(Y0,Y);
  Set(y,Y[0]);
  Set(y0,Y0[0]);
  
  if(verbose > 1) cout << endl;
  Source(Src,Y0,t);
  
  Problem->Transform(Y0,t,dt,YI);
  for(int i=0; i < niterations; i++) {
    for(unsigned int j=0; j < ny; j++) y[j]=y0[j]+halfdt*source[j];
    if(verbose > 1) cout << y[0] << endl;
    Problem->BackTransform(Y,t+dt,dt,YI);
    Source(Src,Y,t);
  }
  for(unsigned int j=0; j < ny; j++) y[j]=y0[j]+dt*source[j];
  if(verbose > 1) cout << y[0] << endl;
  Problem->BackTransform(Y,t+dt,dt,YI);
	
  Problem->Stochastic(Y,t,dt);
  return SUCCESSFUL;
}

void LeapFrog::Predictor(unsigned int start, unsigned int stop)
{
  if(new_y0) {oldhalfdt=lasthalfdt;}
  else yp=yp0;
  double dtprime=halfdt+oldhalfdt;
  Problem->Transform(YP,t-oldhalfdt,dtprime,YP0);
  for(unsigned int j=start; j < stop; j++) yp[j] += dtprime*source0[j];
  Problem->BackTransform(YP,t+halfdt,dtprime,YP0);
  lasthalfdt=halfdt;
}
	
int LeapFrog::Corrector(unsigned int start, unsigned int stop)
{
  CSource(Src,YP,t+halfdt);
  if(dynamic) {
    for(unsigned int j=start; j < stop; j++) {
      y[j]=y0[j]+dt*source[j];
      if(!Active(errmask) || errmask[j]) 
	CalcError(y0[j],y[j],y0[j]+dt*source0[j],y[j]);
    }
  } else for(unsigned int j=start; j < stop; j++)
    y[j]=y0[j]+dt*source[j];

  return 1;
}

void RK2::Predictor(unsigned int start, unsigned int stop)
{
  for(unsigned int j=start; j < stop; j++)
    y[j]=y0[j]+halfdt*source0[j];
  Problem->BackTransform(Y,t+halfdt,halfdt,YI);
}

int RK2::Corrector(unsigned int start, unsigned int stop)
{
  CSource(Src,Y,t+halfdt);
  if(dynamic) {
    for(unsigned int j=start; j < stop; j++) {
      y[j]=y0[j]+dt*source[j];
      if(!Active(errmask) || errmask[j])
	CalcError(y0[j],y[j],y0[j]+dt*source0[j],y[j]);
    }
  } else for(unsigned int j=start; j < stop; j++)
    y[j]=y0[j]+dt*source[j];

  return 1;
}

void RK3::TimestepDependence()
{
  PC::TimestepDependence();
  sixthdt=dt/6.0;
}

void RK3::Predictor(unsigned int start, unsigned int stop)
{
  for(unsigned int j=start; j < stop; j++)
    y[j]=y0[j]+halfdt*source0[j];
  PSource(Src1,Y,t+halfdt);
  for(unsigned int j=start; j < stop; j++)	
    y[j]=y0[j]+dt*(2.0*source1[j]-source0[j]);
}

int RK3::Corrector(unsigned int start, unsigned int stop)
{
  CSource(Src,Y,t+dt);
  if(dynamic) {
    for(unsigned int j=start; j < stop; j++) {
      Var pred=y0[j]+dt*source1[j];
      y[j]=y0[j]+sixthdt*(source0[j]+4.0*source1[j]+source[j]);
      if(!Active(errmask) || errmask[j])
	CalcError(y0[j],y[j],pred,y[j]);
    }
  } else {
    for(unsigned int j=start; j < stop; j++)
      y[j]=y0[j]+sixthdt*(source0[j]+4.0*source1[j]+source[j]);
  }
  return 1;
}

void RK4::Predictor(unsigned int start, unsigned int stop)
{
  for(unsigned int j=start; j < stop; j++) y[j]=y0[j]+halfdt*source0[j];
  Problem->BackTransform(Y,t+halfdt,halfdt,YI);
  PSource(Src1,Y,t+halfdt);
  if(Active(YI)) {swaparray(YI,Y); Set(y,Y[0]);}
  for(unsigned int j=start; j < stop; j++) y[j]=y0[j]+halfdt*source1[j];
  Problem->BackTransform(Y,t+halfdt,halfdt,YI);
  PSource(Src2,Y,t+halfdt);
  if(Active(YI)) {swaparray(YI,Y); Set(y,Y[0]);}
  for(unsigned int j=start; j < stop; j++) y[j]=y0[j]+dt*source2[j];
  Problem->BackTransform(Y,t+dt,dt,YI);
}

int RK4::Corrector(unsigned int start, unsigned int stop)
{
  CSource(Src,Y,t+dt);
  if(Active(YI)) {swaparray(YI,Y); Set(y,Y[0]);}
  if(dynamic) {
    for(unsigned int j=start; j < stop; j++) {
      y[j]=y0[j]+sixthdt*(source0[j]+2.0*(source1[j]+source2[j])+source[j]);
      if(!Active(errmask) || errmask[j])
	CalcError(y0[j],y[j],y0[j]+dt*source2[j],y[j]);
    }
  } else for(unsigned int j=start; j < stop; j++) 
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

void RK5::Predictor(unsigned int start, unsigned int stop)
{
  //#pragma ivdep		
  for(unsigned int j=start; j < stop; j++)
    y[j]=y0[j]+b10*source0[j];
  Problem->BackTransform(Y,t+a1,a1,YI);
  PSource(Src,Y,t+a1);
  if(Active(YI)) {swaparray(YI,Y); Set(y,Y[0]);}
  //#pragma ivdep		
  for(unsigned int j=start; j < stop; j++)
    y[j]=y0[j]+b20*source0[j]+b21*source[j];
  Problem->BackTransform(Y,t+a2,a2,YI);
  PSource(Src2,Y,t+a2);
  if(Active(YI)) {swaparray(YI,Y); Set(y,Y[0]);}
  //#pragma ivdep		
  for(unsigned int j=start; j < stop; j++)
    y[j]=y0[j]+b30*source0[j]+b31*source[j]+b32*source2[j];
  Problem->BackTransform(Y,t+a3,a3,YI);
  PSource(Src3,Y,t+a3);
  if(Active(YI)) {swaparray(YI,Y); Set(y,Y[0]);}
  //#pragma ivdep		
  for(unsigned int j=start; j < stop; j++)
    y[j]=y0[j]+b40*source0[j]+b41*source[j]+b42*source2[j]+b43*source3[j];
  Problem->BackTransform(Y,t+a4,a4,YI);
  PSource(Src4,Y,t+a4);
  if(Active(YI)) {swaparray(YI,Y); Set(y,Y[0]);}
  //#pragma ivdep		
  for(unsigned int j=start; j < stop; j++) 
    y[j]=y0[j]+b50*source0[j]+b51*source[j]+b52*source2[j]+b53*source3[j]+
      b54*source4[j];
  Problem->BackTransform(Y,t+a5,a5,YI);
}

int RK5::Corrector(unsigned int start, unsigned int stop)
{
  CSource(Src,Y,t+a5);
  if(Active(YI)) {swaparray(YI,Y); Set(y,Y[0]);}
  if(dynamic) {
    //#pragma ivdep		
    for(unsigned int j=start; j < stop; j++) {
      y[j]=y0[j]+c0*source0[j]+c2*source2[j]+c3*source3[j]+c5*source[j];
      if(!Active(errmask) || errmask[j]) {
	Var pred=y0[j]+d0*source0[j]+d2*source2[j]+
	  d3*source3[j]+d4*source4[j]+d5*source[j];
	CalcError(y0[j],y[j],pred,y[j]);
      }
    }
  } else {
    //#pragma ivdep		
    for(unsigned int j=start; j < stop; j++) {
      y[j]=y0[j]+c0*source0[j]+c2*source2[j]+c3*source3[j]+c5*source[j];
    }
  }
  return 1;
}
