#include "Exponential.h"

using namespace Array;

// Exponential integrators

typedef array1<Nu>::opt Nuvector;

void E_PC::Predictor()
{
  for(unsigned int j=0; j < Nspecial; j++)	
    y[j]=coeff0[j]*y0[j]+coeff1[j]*source0[j];
  for(unsigned int j=Nspecial; j < ny; j++)
    y[j]=y0[j]+dt*source0[j];
}

int E_PC::Corrector()
{
  Exponential::Source(Src,Y,t+dt);
  if(dynamic) {
    for(unsigned int j=0; j < Nspecial; j++) {
      Var temp=0.5*(source0[j]+source[j]);
      y[j]=coeff0[j]*y0[j]+coeff1[j]*temp;
      if(!Active(errmask) || errmask[j])
	CalcError(y0[j]*dtinv,temp,source0[j],temp);
    }
    for(unsigned int j=Nspecial; j < ny; j++) {
      Var pred=y[j];
      y[j]=y0[j]+halfdt*(source0[j]+source[j]);
      if(!Active(errmask) || errmask[j])
	CalcError(y0[j],y[j],pred,y[j]);
    }
    ExtrapolateTimestep();
  } else {
    for(unsigned int j=0; j < Nspecial; j++)
      y[j]=coeff0[j]*y0[j]+coeff1[j]*0.5*(source0[j]+source[j]);
    for(unsigned int j=Nspecial; j < ny; j++)
      y[j]=y0[j]+halfdt*(source0[j]+source[j]);
  }
  
  return 1;
}

void E_RK2::Predictor()
{
  for(unsigned int j=0; j < Nspecial; j++)	
    y[j]=coeff0[j]*y0[j]+coeff1[j]*source0[j];
  
  for(unsigned int j=Nspecial; j < ny; j++)
    y[j]=y0[j]+dt*source0[j];
}

int E_RK2::Corrector()
{
  Exponential::Source(Src,Y,t+dt);
  if(dynamic) {
    for(unsigned int j=0; j < Nspecial; j++) {
      Var pred=y[j];
      y[j] += coeff2[j]*(source[j]-source0[j]);
      if(!Active(errmask) || errmask[j])
	CalcError(y0[j],y[j],pred,y[j]);
    }
    for(unsigned int j=Nspecial; j < ny; j++) {
      Var pred=y[j];
      y[j]=y0[j]+halfdt*(source0[j]+source[j]);
      if(!Active(errmask) || errmask[j])
	CalcError(y0[j],y[j],pred,y[j]);
    }
    ExtrapolateTimestep();
  } else {
    for(unsigned int j=0; j < Nspecial; j++)
      y[j] += coeff2[j]*(source[j]-source0[j]);
    for(unsigned int j=Nspecial; j < ny; j++)
      y[j]=y0[j]+halfdt*(source0[j]+source[j]);
  }
  
  return 1;
}

void E_RK3::Predictor()
{
  for(unsigned int j=0; j < Nspecial; j++)	
    y[j]=coeff0h[j]*y0[j]+coeff1h[j]*source0[j];
  Exponential::Source(Src1,Y,t+halfdt);
  for(unsigned int j=0; j < Nspecial; j++)	
    y[j]=coeff0[j]*y0[j]+coeff1[j]*(2.0*source1[j]-source0[j]);
  
  for(unsigned int j=Nspecial; j < ny; j++)
    y[j]=y0[j]+dt*source0[j];
}

int E_RK3::Corrector()
{
  Exponential::Source(Src,Y,t+dt);
  if(dynamic) {
    for(unsigned int j=0; j < Nspecial; j++)	
      y[j]=coeff0[j]*y0[j]+coeff1[j]*source0[j];
    Exponential::Source(Src2,Y,t+dt);
    
    for(unsigned int j=0; j < Nspecial; j++) {
      Var pred=y[j]+coeff2[j]*(source2[j]-source0[j]);
      y[j]=coeff0[j]*y0[j]+coeffA[j]*source0[j]+2.0*coeffB[j]*source1[j]
	+coeffC[j]*source[j];
      if(!Active(errmask) || errmask[j]) {
	CalcError(y0[j],y[j],pred,y[j]);
      }
    }
    for(unsigned int j=Nspecial; j < ny; j++) {
      Var pred=y[j];
      y[j]=y0[j]+halfdt*(source0[j]+source[j]);
      if(!Active(errmask) || errmask[j])
	CalcError(y0[j],y[j],pred,y[j]);
    }
    ExtrapolateTimestep();
  } else {
    for(unsigned int j=0; j < Nspecial; j++)
      y[j]=coeff0[j]*y0[j]+coeffA[j]*source0[j]+2.0*coeffB[j]*source1[j]
	+coeffC[j]*source[j];
    for(unsigned int j=Nspecial; j < ny; j++)
      y[j]=y0[j]+halfdt*(source0[j]+source[j]);
  }
  
  return 1;
}

void E_RK4::Predictor()
{
//  msg(ERROR,"This integrator is broken.");
  for(unsigned int j=0; j < Nspecial; j++)	
    y[j]=coeff0h[j]*y0[j]+coeff1h[j]*source0[j];
  Exponential::Source(Src1,Y,t+halfdt);
  for(unsigned int j=0; j < Nspecial; j++)	
    source[j]=coeff0h[j]*y0[j]+coeff1h[j]*source1[j];
  Exponential::Source(Src2,Src,t+halfdt);
  for(unsigned int j=0; j < Nspecial; j++)	
    y[j]=coeff0h[j]*y[j]+coeff1h[j]*(2.0*source2[j]-source0[j]);
  
  for(unsigned int j=Nspecial; j < ny; j++)
    y[j]=y0[j]+dt*source0[j];
}

int E_RK4::Corrector()
{
  Exponential::Source(Src,Y,t+dt);
  
  if(dynamic) {
    for(unsigned int j=0; j < Nspecial; j++)	
      y[j]=coeff0[j]*y0[j]+coeff1[j]*(2.0*source1[j]-source0[j]);
    Exponential::Source(Src3,Y,t+dt);

    for(unsigned int j=0; j < Nspecial; j++) {
      y[j]=coeff0[j]*y0[j]+coeffA[j]*source0[j]
	+coeffB[j]*(source1[j]+source2[j])+coeffC[j]*source[j];
      if(!Active(errmask) || errmask[j]) {
	Var pred=coeff0[j]*y0[j]+coeffA[j]*source0[j]+2.0*coeffB[j]*source1[j]
	  +coeffC[j]*source3[j];
	CalcError(y0[j],y[j],pred,y[j]);
      }
    }
    for(unsigned int j=Nspecial; j < ny; j++) {
      Var pred=y[j];
      y[j]=y0[j]+halfdt*(source0[j]+source[j]);
      if(!Active(errmask) || errmask[j])
	CalcError(y0[j],y[j],pred,y[j]);
    }
    ExtrapolateTimestep();
  } else {
    for(unsigned int j=0; j < Nspecial; j++)
      y[j]=coeff0[j]*y0[j]+coeffA[j]*source0[j]
	+coeffB[j]*(source1[j]+source2[j])+coeffC[j]*source[j];
    for(unsigned int j=Nspecial; j < ny; j++)
      y[j]=y0[j]+halfdt*(source0[j]+source[j]);
  }
  
  return 1;
}
