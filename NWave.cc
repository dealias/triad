#include "options.h"
#include "NWave.h"
#include "Triad.h"

int Nmode;
int NmodeR;
int Ntotal;

Var *psibuffer,*psibufferR,*psibufferStop;
Var *psix,*psiy,*vort;

int Ntriad;
DynVector<Triad> triad;
TriadLimits *triadLimits;

Var *pqbuffer;
Var **pqIndex;
int *qStart;

Real *knorm2,*kfactor;

Nu *nu,*nu_inv;
Real *nuR_inv,*nuI;
Real *forcing;
extern Real tauforce;

static Var random_factor=0.0;
static Real last_t=-REAL_MAX;

void NWave::ConstantForcing(Var *source, double t)
{
	if(t-last_t > tauforce) {last_t=t; crand_gauss(&random_factor);}
#pragma ivdep
	for(int k=0; k < Nmode; k++) source[k] += forcing[k]*random_factor;
}

void NWave::ComputeMoments(Var *source, Var *psi)
{
	Var *p=source, *q=source+Nmode;
#pragma ivdep
	double randre=random_factor.re, randim=random_factor.im;
	for(int i=0; i < Nmode; i++) {
#if COMPLEX
		double psire=psi[i].re;
		double psiim=psi[i].im;
		double sumre=p[i].re*psire;
		double sumim=p[i].im*psiim;
		q[i].re=sumre+sumim;                // S_k*psi_k
		q[i].im=randre*psire+randim*psiim;  // w_k*psi_k
#else
		q[i]=realproduct(p[i],psi[i]);        // S_k*psi_k
#endif		
	}
	q += Nmode;
	if(Nmoment > 1) {
#pragma ivdep
		set(q,psi,Nmode);
		q += Nmode;
		Complex *kstop=psi+Nmode;
#pragma ivdep
		for(int n=2; n < Nmoment; n++) { // psi_k^n
#pragma ivdep
			for(Complex *k=psi; k < kstop; k++, q++) *q=product(*(q-Nmode),*k);
		}
	}
}

void SR::NonLinearSrc(Var *source, Var *psi, double)
{
	set(psibuffer,psi,Nmode);
	
	// Compute reflected psi's
#pragma ivdep		
	for(Var *k=psibuffer; k < psibufferR; k++) conjugate(*(k+Nmode),*k);
	
#if _CRAYMVP
#pragma _CRI taskloop private(ip) value(psibuffer,Ntotal,pqIndex,qStart)
#endif	
	for(int ip=0; ip < Ntotal; ip++) {
		Var *pq=pqIndex[ip]; 
#if COMPLEX
		Real psipre=psibuffer[ip].re;
		Real psipim=psibuffer[ip].im;
#pragma ivdep		
		for(int iq=qStart[ip]; iq < Ntotal; iq++) {
			pq[iq].re=psipre*psibuffer[iq].re-psipim*psibuffer[iq].im;
			pq[iq].im=psipre*psibuffer[iq].im+psipim*psibuffer[iq].re;
		}
#else
		Var psip=psibuffer[ip];
#pragma ivdep		
		for(int iq=qStart[ip]; iq < Ntotal; iq++) pq[iq]=psip*psibuffer[iq];
#endif		
	}
	
#if _CRAYMVP
#pragma _CRI taskloop private(i) value(source,triadLimits,Nmode)
#endif	
	for(int i=0; i < Nmode; i++) {
		Triad *t=triadLimits[i].start, *tstop=triadLimits[i].stop;
#if !_CRAYMVP && COMPLEX && MCREAL
		Real sum0re=0.0, sum0im=0.0, sum1re=0.0, sum1im=0.0,
			sum2re=0.0, sum2im=0.0, sum3re=0.0, sum3im=0.0;
		Triad *tstop0=t+((tstop-t) & 3);
		for(; t < tstop0; t++) {
			sum0re += t->Mkpq * t[0].pq->re;
			sum0im -= t->Mkpq * t[0].pq->im;
		}
		for(; t < tstop; t += 4) {
			sum0re += t[0].Mkpq * t[0].pq->re;
			sum0im -= t[0].Mkpq * t[0].pq->im;
			sum1re += t[1].Mkpq * t[1].pq->re;
			sum1im -= t[1].Mkpq * t[1].pq->im;
			sum2re += t[2].Mkpq * t[2].pq->re;
			sum2im -= t[2].Mkpq * t[2].pq->im;
			sum3re += t[3].Mkpq * t[3].pq->re;
			sum3im -= t[3].Mkpq * t[3].pq->im;
		}
		source[i].re=sum0re+sum1re+sum2re+sum3re;
		source[i].im=sum0im+sum1im+sum2im+sum3im;
#else		
		Var sum;
		for(sum=0.0; t < tstop; t++) sum += t->Mkpq * conj(*(t->pq));
		source[i]=sum;
#endif		
	}
	if(Nmoment) ComputeMoments(source,psi);
	ConstantForcing(source,t);
}

void NWave::LinearSrc(Var *source, Var *psi, double)
{
#pragma ivdep
	for(int k=0; k < Nmode; k++) source[k] -= nu[k]*psi[k];
}

void NWave::ExponentialLinearity(Var *source, Var *, double)
{
#pragma ivdep
	for(int k=0; k < Nmode; k++) source[k] *= nu_inv[k];
}

void NWave::ConservativeExponentialLinearity(Real *source, Real *, double)
{
#pragma ivdep
	for(int k=0; k < Nmode; k++) source[k] *= nuR_inv[k];
}

void NWave::ConservativeExponentialLinearity(Complex *source, Complex *psi,
											 double)
{
#pragma ivdep
	for(int k=0; k < Nmode; k++) {
		source[k].re += imag(nu[k])*imag(psi[k]);
		source[k].im -= imag(nu[k])*real(psi[k]);
		source[k] *= nuR_inv[k];
	}
}

Solve_RC C_Euler::Solve(double t, double dt)
{
	int j,iter,cont,jfix=-1;
	Real *ry0=(Real *) y0, *ry=(Real *) y;
	Real *rsource=(Real *)source, *rlastdiff=(Real *)lastdiff;
	double tau,mu,temp;
	int ny0=nyprimary*sizeof(Var)/sizeof(Real);
	
	tau=dt;
	Source(source,y0,t);
	
	mu=temp=0.0;
	for(j=0; j < ny0; j++) {
		if(ry0[j]) {
			temp=tau*rsource[j]/ry0[j];
			if(fabs(temp) > fabs(mu)) {mu=temp;	if(fabs(mu) >= 0.5) jfix=j;}
		}
		else if(rsource[j]) jfix=j;
	}
	mu += 2.0;
	
	if(jfix == -1)	// Evolve forwards.
		for(j=0; j < ny0; j++)
			ry0[j]=sgn(ry0[j])*sqrt(ry0[j]*(ry0[j]+mu*tau*rsource[j]));
	else {			// Iterate backwards.
		set(y,y0,ny);
		iter=0;
		for(j=0; j < ny0; j++) rlastdiff[j]=0.0;
		ry0[jfix]=ry[jfix]+tau*rsource[jfix];
		temp=(ry[jfix]/ry0[jfix]+1.0)*rsource[jfix];
		do {
			cont=0;
			Source(source,y0,t);
			mu=temp/rsource[jfix];
			for(j=0; j < ny0; j++) {
				Real yj=ry0[j];
				if(j !=jfix) {
					Real discr=ry[j]*ry[j]+mu*tau*rsource[j]*ry0[j];
					if(discr >= 0.0) ry0[j]=sgn(ry[j])*sqrt(discr);
					else {
						set(y0,y,ny);
						return NONINVERTIBLE;
					}
				}
				Real diff=abs(ry0[j]-yj);
				if(diff != 0.0 && 
				   (diff < rlastdiff[j] || rlastdiff[j] == 0.0)) cont=1;
				rlastdiff[j]=diff;
			}
			if(++iter == 1000) msg(ERROR,"Iteration did not converge");
		} while (cont);
	}
	for(j=nyprimary; j < ny; j++) y0[j] += dt*source[j];
	return SUCCESSFUL;
}

inline int CorrectC_PC::Correct(const Real y0, const Real y1, Real& y,
								const Real source0, const Real source,
								const double dt)
{
	Real discr=y0*y0+dt*(y0*source0+y1*source);
	if(discr < 0.0) return 0;
	y=sgn(y1)*sqrt(discr);
	return 1;
}

inline int CorrectC_PC::Correct(const Complex y0, const Complex y1, Complex& y,
								const Complex source0, const Complex source,
								const double dt)
{
	if(!Correct(y0.re,y1.re,y.re,source0.re,source.re,dt)) return 0;
	if(!Correct(y0.im,y1.im,y.im,source0.im,source.im,dt)) return 0;
	return 1;
}

int C_PC::Corrector(double dt, int dynamic, int start, int stop)
{
	int j;
	for(j=start; j < stop; j++)
		if(!Correct(y0[j],y1[j],y[j],source0[j],source[j],dt)) return 0;
	if(dynamic) 
		for(j=start; j < stop; j++)
			if(!errmask || errmask[j]) 
				CalcError(y0[j],y[j],y1[j],y[j]);
	return 1;
}

void E_PC::Allocate(int n)
{
	PC::Allocate(n);
	expinv=new Nu[nyprimary];
	onemexpinv=new Nu[nyprimary];
	
	nu_inv=new Nu[nyprimary];

	for(int j=0; j < nyprimary; j++) {
		if(nu[j] != 0.0) nu_inv[j]=1.0/nu[j];
		else nu_inv[j]=1.0;
	}
}

void E_PC::TimestepDependence(double dt)
{
	for(int j=0; j < nyprimary; j++) {
		onemexpinv[j]=-expm1(-nu[j]*dt);
		expinv[j]=1.0-onemexpinv[j];
		if(nu[j] == 0.0) onemexpinv[j]=dt;
	}
	dtinv=1.0/dt;
}

void E_PC::Predictor(double t, double, int start, int stop)
{
	for(int j=start; j < stop; j++)	
		y1[j]=expinv[j]*y0[j]+onemexpinv[j]*source0[j];
	Source(source,y1,t+dt);
}

int E_PC::Corrector(double, int dynamic, int start, int stop)
{
	int j;
	for(j=start; j < stop; j++) {
		source[j]=0.5*(source0[j]+source[j]);
		y[j]=expinv[j]*y0[j]+onemexpinv[j]*source[j];
	}
	if(dynamic)
		for(j=start; j < stop; j++)
			if(!errmask || errmask[j]) 
				CalcError(y0[j]*dtinv,source[j],source0[j],source[j]);
	return 1;
}

void I_PC::Allocate(int n)
{
	PC::Allocate(n);
	expinv=new Nu[nyprimary];
}

void I_PC::TimestepDependence(double dt)
{
	for(int j=0; j < nyprimary; j++) expinv[j]=exp(-nu[j]*dt);
}

void I_PC::Predictor(double t, double dt, int start, int stop)
{
	for(int j=start; j < stop; j++) y1[j]=(y0[j]+dt*source0[j])*expinv[j];
	Source(source,y1,t+dt);
}

int I_PC::Corrector(double dt, int dynamic, int start, int stop)
{
	const double halfdt=0.5*dt;
	if(dynamic) for(int j=start; j < stop; j++) {
		Var pred=y[j];
		y[j]=y0[j]*expinv[j]+halfdt*(source0[j]*expinv[j]+source[j]);
		if(!errmask || errmask[j]) 
			CalcError(y0[j]*expinv[j],y[j],pred,y[j]);
	} else for(int j=start; j < stop; j++) {
		y[j]=y0[j]*expinv[j]+halfdt*(source0[j]*expinv[j]+source[j]);
	}
	return 1;
}

void CE_PC::Allocate(int n)
{
	PC::Allocate(n);
	expinv=new Nu[nyprimary];
	onemexpinv=new Real[nyprimary];
	
	nuR_inv=new Real[nyprimary];
	nuI=new Real[nyprimary];

	for(int j=0; j < nyprimary; j++) {
		nuI[j]=imag(nu[j]);
		if(real(nu[j]) != 0.0) nuR_inv[j]=1.0/real(nu[j]);
		else nuR_inv[j]=1.0;
	}
}

void CE_PC::TimestepDependence(double dt)
{
	for(int j=0; j < nyprimary; j++) {
		onemexpinv[j]=-expm1(-real(nu[j])*dt);
		expinv[j]=exp(-nu[j]*dt);
		if(real(nu[j]) == 0.0) onemexpinv[j]=dt;
	}
	dtinv=1.0/dt;
}

void CE_PC::Predictor(double t, double, int start, int stop)
{
	for(int j=start; j < stop; j++)
		y[j]=expinv[j]*y0[j]+onemexpinv[j]*source0[j];
	Source(source,y,t+dt);
}

int CE_PC::Corrector(double, int dynamic, int start, int stop)
{
	int j;
	for(j=start; j < stop; j++)
		if(!Correct(expinv[j]*y0[j],y1[j],y[j],source0[j],source[j],
					onemexpinv[j])) return 0;
	if(dynamic) for(j=start; j < stop; j++) {
		Var corr=0.5*(source0[j]+source[j]);
		if(!errmask || errmask[j]) 
			CalcError(y0[j]*dtinv,corr,source0[j],corr);
	}
	return 1;
}

void I_RK2::Allocate(int n)
{
	RK2::Allocate(n);
	expinv=new Nu[nyprimary];
}

void I_RK2::TimestepDependence(double dt)
{
	RK2::TimestepDependence(dt);
	for(int j=0; j < nyprimary; j++) expinv[j]=exp(-nu[j]*halfdt);
}

void I_RK2::Predictor(double t, double, int start, int stop)
{
	for(int j=start; j < stop; j++) y1[j]=(y0[j]+halfdt*source0[j])*expinv[j];
	Source(source,y1,t+halfdt);
}

int I_RK2::Corrector(double dt, int dynamic, int start, int stop)
{
	int j;
	for(j=start; j < stop; j++)
		y[j]=y0[j]*expinv[j]+dt*source[j];
	if(dynamic) for(j=start; j < stop; j++)
		if(!errmask || errmask[j]) 
			CalcError(y0[j]*expinv[j],y[j],
					  (y0[j]+dt*source0[j])*expinv[j],y[j]);
	for(j=start; j < stop; j++) y[j] *= expinv[j];
	return 1;
}

inline int C_RK2::Correct(const Real Y0, const Real Y1, Real& Y,
						  const Real, const Real Source,
						  const double)
{
	Real temp=dt*Source;
	Real discr=Y0*Y0+2.0*Y1*temp;
	if(discr < 0.0) return 0;
	Y=sgn(Y0+temp)*sqrt(discr);
	return 1;
}

inline int C_RK2::Correct(const Complex Y0, const Complex Y1, Complex& Y,
						  const Complex Source0, const Complex Source,
						  const double dt)
{
	if(!Correct(Y0.re,Y1.re,Y.re,Source0.re,Source.re,dt)) return 0;
	if(!Correct(Y0.im,Y1.im,Y.im,Source0.im,Source.im,dt)) return 0;
	return 1;
}

int C_RK2::Corrector(double dt, int dynamic, int start, int stop)
{
	int j;
	for(j=start; j < stop; j++)
		if(!Correct(y0[j],y1[j],y[j],source0[j],source[j],dt)) return 0;
	if(dynamic)
		for(j=start; j < stop; j++)
			if(!errmask || errmask[j]) 
				CalcError(y0[j],y[j],y0[j]+dt*source0[j],y[j]);
	return 1;
}

void I_RK4::Allocate(int n)
{
	RK4::Allocate(n);
	expinv=new Nu[nyprimary];
}

void I_RK4::TimestepDependence(double dt)
{
	RK4::TimestepDependence(dt);
	thirddt=dt/3.0;
	for(int j=0; j < nyprimary; j++) expinv[j]=exp(-nu[j]*halfdt);
	dtinv=1.0/dt;
}

void I_RK4::Predictor(double t, double dt, int start, int stop)
{
	int j;
	for(j=start; j < stop; j++) y[j]=(y0[j]+halfdt*source0[j])*expinv[j];
	Source(source1,y,t+halfdt);
	for(j=start; j < stop; j++) y[j]=y0[j]*expinv[j]+halfdt*source1[j];
	Source(source2,y,t+halfdt);
	for(j=start; j < stop; j++) y[j]=(y0[j]*expinv[j]+dt*source2[j])*expinv[j];
	Source(source,y,t+dt);
}

int I_RK4::Corrector(double, int dynamic, int start, int stop)
{
	int j;
	for(j=start; j < stop; j++) {
		y[j]=((y0[j]+sixthdt*source0[j])*expinv[j]+
			  thirddt*(source1[j]+source2[j]))*expinv[j]+sixthdt*source[j];
	}
	if(dynamic)
		for(j=start; j < stop; j++)
			if(!errmask || errmask[j]) 
				CalcError(y0[j]*dtinv,source[j],source2[j]*expinv[j],
						  source[j]);
	return 1;
}

void C_RK4::TimestepDependence(double dt)
{
	halfdt=0.5*dt;
	sixthdt=dt/6.0;
	sixthdt2=sixthdt*dt;
}

inline int C_RK4::Correct(const Real Y0, Real& Y, const Real Source0,
						  const Real Source1, const Real Source2,
						  const Real Source, const double)
{			
	Real temp=sixthdt*(Source0+2.0*(Source1+Source2)+Source);
	Real discr=Y0*Y0+2.0*(Y0*temp+sixthdt2*(Source1*(Source0+Source2)+
											Source2*Source));
	if(discr < 0.0) return 0;
	Y=sgn(Y0+temp)*sqrt(discr);
	return 1;
}

inline int C_RK4::Correct(const Complex Y0, Complex& Y, const Complex Source0,
						  const Complex Source1, const Complex Source2,
						  const Complex Source, const double dt)
{			
	if(!Correct(Y0.re,Y.re,Source0.re,Source1.re,Source2.re,Source.re,dt))
		return 0;
	if(!Correct(Y0.im,Y.im,Source0.im,Source1.im,Source2.im,Source.im,dt))
		return 0;
	return 1;
}

int C_RK4::Corrector(double dt, int dynamic, int start, int stop)
{
	int j;
	if(dynamic) for(j=start; j < stop; j++) {
		Var pred=y[j];
		if(!Correct(y0[j],y[j],source0[j],source1[j],source2[j],source[j],dt))
			return 0;
		if(!errmask || errmask[j]) 
			CalcError(y0[j],y[j],pred,y[j]);
	}
	else for(j=start; j < stop; j++) {
		if(!Correct(y0[j],y[j],source0[j],source1[j],source2[j],source[j],dt))
			return 0;
	}
	return 1;
}

void I_RK5::Allocate(int n)
{
	RK5::Allocate(n);
	expinv=new Nu[nyprimary];
	expinv5=new Nu[nyprimary];
}

void I_RK5::Predictor(double t, double, int start, int stop)
{
	int j;
#pragma ivdep		
	for(j=start; j < stop; j++) {
		expinv[j]=exp(-nu[j]*a1);
		y[j]=(y0[j]+b10*source0[j])*expinv[j];
	}
	Source(source,y,t+a1);
#pragma ivdep		
	for(j=start; j < stop; j++) {
		source[j] /= expinv[j];
		expinv[j]=exp(-nu[j]*a2);
		y2[j]=(y0[j]+b20*source0[j]+b21*source[j])*expinv[j];
	}
	Source(source2,y2,t+a2);
#pragma ivdep		
	for(j=start; j < stop; j++) {
		source2[j] /= expinv[j];
		expinv[j]=exp(-nu[j]*a3);
		y3[j]=(y0[j]+b30*source0[j]+b31*source[j]+b32*source2[j])*expinv[j];
	}
	Source(source3,y3,t+a3);
#pragma ivdep		
	for(j=start; j < stop; j++) {
		source3[j] /= expinv[j];
		expinv[j]=exp(-nu[j]*a4);
		y4[j]=(y0[j]+b40*source0[j]+b41*source[j]+b42*source2[j]+
			   b43*source3[j])*expinv[j];
	}
	Source(source4,y4,t+a4);
#pragma ivdep		
	for(j=start; j < stop; j++) {
		source4[j] /= expinv[j];
		expinv5[j]=exp(-nu[j]*a5);
		y[j]=(y0[j]+b50*source0[j]+b51*source[j]+b52*source2[j]+b53*source3[j]+
			  b54*source4[j])*expinv5[j];
	}
	Source(source,y,t+a5);
#pragma ivdep	
	for(j=start; j < stop; j++) source[j] /= expinv5[j];
}

int I_RK5::Corrector(double dt, int dynamic, int start, int stop)
{
	RK5::Corrector(dt,dynamic,start,stop);
#pragma ivdep	
	for(int j=start; j < stop; j++) y[j] *= expinv[j];
	return 1;
}

inline void C_RK5::Correct(const Real Y0, Real& Y2, Real& Y3,
						   const Real Y4, Real& Y,
						   const Real Source0, const Real Source2, 
						   const Real Source3, const Real Source4,
						   const Real Source, const double, int& invertible)
{
	Real discr=Y0*Y0+2.0*(c0*Y0*Source0+c2*Y2*Source2+c3*Y3*Source3+
						  c5*Y*Source);
	if(discr < 0.0) invertible=0;
	else {
// Put discr in Y2, pred in Y3 for deferred error analysis
		Y3=Y0*Y0+2.0*(d0*Y0*Source0+d2*Y2*Source2+d3*Y3*Source3+
					  d4*Y4*Source4+d5*Y*Source);
		Y2=discr;
		Y=sgn(Y0+c0*Source0+c2*Source2+c3*Source3+c5*Source)*sqrt(discr);
	}
}

inline void C_RK5::Correct(const Complex Y0, Complex& Y2, Complex& Y3,
						   const Complex Y4, Complex& Y,
						   const Complex Source0, const Complex Source2, 
						   const Complex Source3, const Complex Source4,
						   const Complex Source, const double dt,
						   int& invertible)
{
	Correct(Y0.re,Y2.re,Y3.re,Y4.re,Y.re,Source0.re,Source2.re,
			Source3.re,Source4.re,Source.re,dt,invertible);
#if !_CRAYMVP		
	if(!invertible) return;
#endif			
	Correct(Y0.im,Y2.im,Y3.im,Y4.im,Y.im,Source0.im,Source2.im,
			Source3.im,Source4.im,Source.im,dt,invertible);
}

int C_RK5::Corrector(double, int dynamic, int start, int stop)
{
	int j,invertible=1;
#pragma ivdep
	for(j=start; j < stop; j++) {
		Correct(y0[j],y2[j],y3[j],y4[j],y[j],source0[j],source2[j],
				source3[j],source4[j],source[j],dt,invertible);
#if !_CRAYMVP		
		if(!invertible) return 0;
#endif			
	}
#if _CRAYMVP
    if(!invertible) return 0;
#endif
	
	if(dynamic) {
		for(j=start; j < stop; j++) 
			if(!errmask || errmask[j]) 
				CalcError(y0[j]*y0[j],y2[j],y3[j],y2[j]);
		ExtrapolateTimestep();
	}
	return 1;
}

void NWave::DisplayInvariants(Real E, Real Z, Real P)
{
	cout << "Energy = " << E << newl;
	cout << "Enstrophy = " << Z << newl;
	cout << "Palinstrophy = " << P << endl;
}

void NWave::FinalOutput()
{
	Real E,Z,P;
	int i;
	const int Nlimit=128;
	
	cout << newl << "FINAL VALUES:" << newl << endl;
	if(Nmode <= Nlimit || verbose > 1) {
		for(i=0; i < Nmode; i++) cout << "psi[" << i << "] = " << y[i] << newl;
		cout << endl;
	}
	
	ComputeInvariants(y,Nmode,E,Z,P);
	DisplayInvariants(E,Z,P);
	
	if(Nmoment > 2 && t) {
		cout << newl << "AVERAGED VALUES:" << newl << endl;
// We overwrite y+3*Nmode here, since it is no longer needed.
		Var *y2=y+3*Nmode;
		for(i=0; i < Nmode; i++) y2[i] = (real(y2[i])+imag(y2[i]))/t;
		
		if(Nmode <= Nlimit || verbose > 1) {
			for(i=0; i < Nmode; i++) {
				Real y2avg=real(y2[i]);
				cout << "|psi|^2 [" << i << "] = " << y2avg << newl;
			}
			cout << endl;
		}
		
		for(i=0; i < Nmode; i++) y2[i] = sqrt(real(y2[i]));
		ComputeInvariants(y2,Nmode,E,Z,P);
		DisplayInvariants(E,Z,P);
	}
}
