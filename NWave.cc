#include "NWave.h"
#include "Cartesian.h"

int Npsi;
int NpsiR;
int Ntotal;

// (0 => evolve all modes, 1 => evolve only half of the modes).
int reality=1; // Reality condition flag 

Var *psibuffer,*psibuffer0,*psibufferR,*psibufferStop,*pqbuffer,*psitemp;
Var *convolution,*convolution0;

Var **pqIndex;
int *qStart;

int Ntriad;
DynVector<Triad> triad;
TriadLimits *triadLimits;

Real *kinv2;

Nu *nu,*nu_inv;
Real *nuR_inv,*nuI;
Real *forcing;
extern Real tauforce;

void ComputeMoments(Var *source, Var *psi) {
	Var *k, *q=source+Npsi, *kstop=psi+Npsi;
	for(k=psi; k < kstop; k++, q++) *q=product(*(q-Npsi),*k); // S_k*psi_k
#pragma ivdep
	if(Nmoment > 1) {
		for(k=psi; k < kstop; k++, q++) *q=*k; // psi_k
		for(int n=2; n < Nmoment; n++) { // psi_k^n
#pragma ivdep
			for(k=psi; k < kstop; k++, q++) *q=product(*(q-Npsi),*k);
		}
	}
}

void PrimitiveNonlinearitySR(Var *source, Var *psi, double)
{
	set(psibuffer,psi,Npsi);
	
	// Compute reflected psi's
#pragma ivdep		
	for(Var *k=psibuffer; k < psibufferR; k++) conjugate(*(k+Npsi),*k);
	
#if (_AIX || __GNUC__) && COMPLEX
	Var *pq=pqbuffer,*p;
	for(p=psibuffer; p < psibufferR; p++) {
		Real psipre=p->re, psipim=p->im;
		for(Var *q=psibufferR; q < psibufferStop; q++, pq++) {
			pq->re=psipre*q->re-psipim*q->im;
			pq->im=psipre*q->im+psipim*q->re;
		}
	}
	for(p=psibufferR; p < psibufferStop; p++) {
		Real psipre=p->re, psipim=p->im;
		for(Var *q=p; q < psibufferStop; q++, pq++) {
			pq->re=psipre*q->re-psipim*q->im;
			pq->im=psipre*q->im+psipim*q->re;
		}
	}
#else
#pragma _CRI taskloop private(ip) value(psibuffer,Ntotal,pqIndex,qStart)
	for(int ip=0; ip < Ntotal; ip++) {
		Var psip=psibuffer[ip];
		Var *pq=pqIndex[ip]; 
#pragma ivdep		
		for(int iq=qStart[ip]; iq < Ntotal; iq++) pq[iq]=psip*psibuffer[iq];
	}
#endif		
	
#pragma _CRI taskloop private(i) value(source,triadLimits,Npsi)
	for(int i=0; i < Npsi; i++) {
		Triad *t=triadLimits[i].start, *tstop=triadLimits[i].stop;
#if (_AIX || __GNUC__) && COMPLEX && MCREAL
		Real sumre, sumim;
		for(sumre=sumim=0.0; t < tstop; t++) {
			sumre += t->Mkpq * (*(t->pq)).re;
			sumim -= t->Mkpq * (*(t->pq)).im;
		}
		source[i].re=sumre;
		source[i].im=sumim;
#else		
		Var sum;
		for(sum=0.0; t < tstop; t++) sum += t->Mkpq * conj(*(t->pq));
		source[i]=sum;
#endif		
	}
	if(Nmoment) ComputeMoments(source,psi);
}

Cartesian *CartesianMode;

void PrimitiveNonlinearity(Var *source, Var *psi, double)
{
	int i;
	set(psibuffer,psi,Npsi);
	
	// Compute reflected psi's
#pragma ivdep		
	for(Var *p=psibuffer; p < psibufferR; p++) conjugate(*(p+Npsi),*p);
	
	for(i=0; i < Npsi; i++) source[i]=0.0;
		
	int q=0, nn=2*Npsi;
	for(int k=0; k < Npsi; k++) {
		Real kx=CartesianMode[k].X();
		Real ky=CartesianMode[k].Y();
		for(int p=0; p < nn; p++)	{
			Real px=CartesianMode[p].X();
			Real py=CartesianMode[p].Y();
			Cartesian mq = CartesianMode[k]-CartesianMode[p];
			if(q < nn-1 && CartesianMode[q+1] == mq) q++;
			else if(q > 1 && CartesianMode[q-1] == mq) q--;
			else for(q=0; q < nn && CartesianMode[q] != mq; q++);
			if(q < nn) source[k] += (kx*py-ky*px)*(px*px+py*py)*
				psibuffer[p]*psibuffer[q];
		}
	}

	for(i=0; i < Npsi; i++) source[i] *= kinv2[i];
	if(Nmoment) ComputeMoments(source,psi);
}


void PrimitiveNonlinearityFFT(Complex *source, Complex *psi, double)
{
	int i;
	extern Cartesian *CartesianMode;
	
	*psibuffer0=0.0;
	CartesianPad(psibuffer,psi);
	
#pragma ivdep	
	for(i=0; i < Npsi; i++) {
		source[i]=psi[i]*I*CartesianMode[i].K2();
		psitemp[i]=source[i]*CartesianMode[i].Y();
	}
	*convolution0=0.0;
	CartesianPad(convolution,psitemp);
	
	convolve(convolution0,convolution0,psibuffer0,Npsibuffer+1,log2n);
	CartesianUnPad(psitemp,convolution);
	
#pragma ivdep	
	for(i=0; i < Npsi; i++) source[i] *= CartesianMode[i].X();
	*convolution0=0.0;
	CartesianPad(convolution,source);
	
	// Reuse FFT of psibuffer0.
	convolve0(convolution0,convolution0,psibuffer0,Npsibuffer+1,log2n);
	CartesianUnPad(psibuffer,convolution);
	
#pragma ivdep	
	for(i=0; i < Npsi; i++) {
		source[i].re = -kinv2[i]*(CartesianMode[i].Y()*psibuffer[i].im-
								  CartesianMode[i].X()*psitemp[i].im);
		source[i].im = kinv2[i]*(CartesianMode[i].Y()*psibuffer[i].re-
								 CartesianMode[i].X()*psitemp[i].re);
	}
	if(Nmoment) ComputeMoments(source,psi);
}

void PrimitiveNonlinearityFFT(Real *, Real *, double)
{
	msg(ERROR,"Pseudospectral approximation requires COMPLEX=1");
}
	
void StandardLinearity(Var *source, Var *psi, double)
{
#pragma ivdep
	for(int k=0; k < Npsi; k++) source[k] -= nu[k]*psi[k];
}

void ExponentialLinearity(Var *source, Var *, double)
{
#pragma ivdep
	for(int k=0; k < Npsi; k++) source[k] *= nu_inv[k];
}

void ConservativeExponentialLinearity(Real *source, Real *, double)
{
#pragma ivdep
	for(int k=0; k < Npsi; k++) source[k] *= nuR_inv[k];
}

void ConservativeExponentialLinearity(Complex *source, Complex *psi, double)
{
#pragma ivdep
	for(int k=0; k < Npsi; k++) {
		source[k].re += imag(nu[k])*imag(psi[k]);
		source[k].im -= imag(nu[k])*real(psi[k]);
		source[k] *= nuR_inv[k];
	}
}

static Real last_t=-REAL_MAX;
static Complex randomfactor=0.0;

void ConstantForcing(Var *source, Var *, double t)
{
	if(t-last_t > tauforce) {last_t=t; crand_gauss(&randomfactor);}
#pragma ivdep
	for(int k=0; k < Npsi; k++) source[k] += forcing[k]*randomfactor;
}

#if COMPLEX	
Solve_RC C_Euler::Solve(double, double)
{
	msg(ERROR,"Complex C_Euler has not been implemented");
	return UNSUCCESSFUL;
}
#else
Solve_RC C_Euler::Solve(double t, double dt)
{
	int j,iter,cont,jfix=-1;
	Real yj;
	double tau,mu,temp;
	
	tau=dt;
	Source(source,y0,t);
	
	mu=temp=0.0;
	for(j=0; j < nyprimary; j++) {
		if(y0[j] != 0.0) temp=tau*source[j]/y0[j];
		else jfix=j;
		if(fabs(temp) > fabs(mu)) {mu=temp;	if(fabs(mu) >= 0.5) jfix=j;}
	}
	mu += 2.0;
	
	if(jfix == -1)	// Evolve forwards.
		for(j=0; j < nyprimary; j++)
			y0[j]=sgn(y0[j])*sqrt(y0[j]*(y0[j]+mu*tau*source[j]));
	else {			// Iterate backwards.
		set(y,y0,nyprimary);
		iter=0;
		for(j=0; j < nyprimary; j++) lastdiff[j]=0.0;
		y0[jfix]=y[jfix]+tau*source[jfix];
		temp=(y[jfix]/y0[jfix]+1.0)*source[jfix];
		do {
			cont=0;
			Source(source,y0,t);
			mu=temp/source[jfix];
			for(j=0; j < nyprimary; j++) {
				yj=y0[j];
				if(j !=jfix) {
					Real discr=y[j]*y[j]+mu*tau*source[j]*y0[j];
					if(discr >= 0.0) y0[j]=sgn(y[j])*sqrt(discr);
					else msg(ERROR,"Negative discriminant encountered");
				}
				Real diff=abs(y0[j]-yj);
				if(diff != 0.0 && 
				   (diff < lastdiff[j] || lastdiff[j] == 0.0)) cont=1;
				lastdiff[j]=diff;
			}
			if(++iter == 100) msg(ERROR,"Iteration did not converge");
		} while (cont);
	}
	for(j=nyprimary; j < ny; j++) y0[j] += dt*source[j];
	return SUCCESSFUL;
}
#endif	

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

int C_PC::Corrector(double dt, double& errmax, int start, int stop)
{
	int j;
	for(j=start; j < stop; j++)
		if(!Correct(y0[j],y1[j],y[j],source0[j],source[j],dt)) return 0;
	if(dynamic) for(j=start; j < stop; j++)
		calc_error(y0[j],y[j],y1[j],y[j],errmax);
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

int E_PC::Corrector(double, double& errmax, int start, int stop)
{
	int j;
	for(j=start; j < stop; j++) {
		source[j]=0.5*(source0[j]+source[j]);
		y[j]=expinv[j]*y0[j]+onemexpinv[j]*source[j];
	}
	if(dynamic) for(j=start; j < stop; j++)
		calc_error(y0[j]*dtinv,source[j],source0[j],source[j],errmax);
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

int I_PC::Corrector(double dt, double& errmax, int start, int stop)
{
	const double halfdt=0.5*dt;
	if(dynamic) for(int j=start; j < stop; j++) {
		Var pred=y[j];
		y[j]=y0[j]*expinv[j]+halfdt*(source0[j]*expinv[j]+source[j]);
		calc_error(y0[j]*expinv[j],y[j],pred,y[j],errmax);
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

int CE_PC::Corrector(double, double& errmax, int start, int stop)
{
	int j;
	for(j=start; j < stop; j++)
		if(!Correct(expinv[j]*y0[j],y1[j],y[j],source0[j],source[j],
					onemexpinv[j])) return 0;
	if(dynamic) for(j=start; j < stop; j++) {
		Var corr=0.5*(source0[j]+source[j]);
		calc_error(y0[j]*dtinv,corr,source0[j],corr,errmax);
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

int I_RK2::Corrector(double dt, double& errmax, int start, int stop)
{
	int j;
	for(j=start; j < stop; j++)
		y[j]=y0[j]*expinv[j]+dt*source[j];
	if(dynamic) for(j=start; j < stop; j++)
		calc_error(y0[j]*expinv[j],y[j],
				   (y0[j]+dt*source0[j])*expinv[j],y[j],errmax);
	for(j=start; j < stop; j++) y[j] *= expinv[j];
	return 1;
}

inline int C_RK2::Correct(const Real y0, const Real y1, Real& y,
						  const Real, const Real source,
						  const double)
{
	Real temp=dt*source;
	Real discr=y0*y0+2.0*y1*temp;
	if(discr < 0.0) return 0;
	y=sgn(y0+temp)*sqrt(discr);
	return 1;
}

inline int C_RK2::Correct(const Complex y0, const Complex y1, Complex& y,
						  const Complex source0, const Complex source,
						  const double dt)
{
	if(!Correct(y0.re,y1.re,y.re,source0.re,source.re,dt)) return 0;
	if(!Correct(y0.im,y1.im,y.im,source0.im,source.im,dt)) return 0;
	return 1;
}

int C_RK2::Corrector(double dt, double& errmax, int start, int stop)
{
	int j;
	for(j=start; j < stop; j++)
		if(!Correct(y0[j],y1[j],y[j],source0[j],source[j],dt)) return 0;
	if(dynamic) for(j=start; j < stop; j++)
		calc_error(y0[j],y[j],y0[j]+dt*source0[j],y[j],errmax);
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

int I_RK4::Corrector(double, double& errmax, int start, int stop)
{
	int j;
	for(j=start; j < stop; j++) {
		y[j]=((y0[j]+sixthdt*source0[j])*expinv[j]+
			  thirddt*(source1[j]+source2[j]))*expinv[j]+sixthdt*source[j];
	}
	if(dynamic) for(j=start; j < stop; j++)
		calc_error(y0[j]*dtinv,source[j],source2[j]*expinv[j],source[j],
				   errmax);
	return 1;
}

void C_RK4::TimestepDependence(double dt)
{
	halfdt=0.5*dt;
	sixthdt=dt/6.0;
	sixthdt2=sixthdt*dt;
}

inline int C_RK4::Correct(const Real y0, Real& y, const Real source0,
						  const Real source1, const Real source2,
						  const Real source, const double)
{			
	Real temp=sixthdt*(source0+2.0*(source1+source2)+source);
	Real discr=y0*y0+2.0*(y0*temp+sixthdt2*(source1*(source0+source2)+
											source2*source));
	if(discr < 0.0) return 0;
	y=sgn(y0+temp)*sqrt(discr);
	return 1;
}

inline int C_RK4::Correct(const Complex y0, Complex& y, const Complex source0,
						  const Complex source1, const Complex source2,
						  const Complex source, const double dt)
{			
	if(!Correct(y0.re,y.re,source0.re,source1.re,source2.re,source.re,dt))
		return 0;
	if(!Correct(y0.im,y.im,source0.im,source1.im,source2.im,source.im,dt))
		return 0;
	return 1;
}

int C_RK4::Corrector(double dt, double& errmax, int start, int stop)
{
	int j;
	if(dynamic) for(j=start; j < stop; j++) {
		Var pred=y[j];
		if(!Correct(y0[j],y[j],source0[j],source1[j],source2[j],source[j],dt))
			return 0;
		calc_error(y0[j],y[j],pred,y[j],errmax);
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
	expinv1=new Nu[nyprimary];
	expinv2=new Nu[nyprimary];
	expinv3=new Nu[nyprimary];
	expinv4=new Nu[nyprimary];
	expinv5=new Nu[nyprimary];
}

void I_RK5::TimestepDependence(double dt)
{
	RK5::TimestepDependence(dt);
	for(int j=0; j < nyprimary; j++) {
		expinv1[j]=exp(-nu[j]*a1);
		expinv2[j]=exp(-nu[j]*a2);
		expinv3[j]=exp(-nu[j]*a3);
		expinv4[j]=exp(-nu[j]*a4);
		expinv5[j]=exp(-nu[j]*a5);
	}
}

void I_RK5::Predictor(double t, double, int start, int stop)
{
	int j;
#pragma ivdep		
	for(j=start; j < stop; j++) y[j]=(y0[j]+b10*source0[j])*expinv1[j];
	Source(source,y,t+a1);
#pragma ivdep		
	for(j=start; j < stop; j++) {
		source[j] /= expinv1[j];
		y2[j]=(y0[j]+b20*source0[j]+b21*source[j])*expinv2[j];
	}
	Source(source2,y2,t+a2);
#pragma ivdep		
	for(j=start; j < stop; j++) {
		source2[j] /= expinv2[j];
		y3[j]=(y0[j]+b30*source0[j]+b31*source[j]+b32*source2[j])*expinv3[j];
	}
	Source(source3,y3,t+a3);
#pragma ivdep		
	for(j=start; j < stop; j++) {
		source3[j] /= expinv3[j];
		y4[j]=(y0[j]+b40*source0[j]+b41*source[j]+b42*source2[j]+
			   b43*source3[j])*expinv4[j];
	}
	Source(source4,y4,t+a4);
#pragma ivdep		
	for(j=start; j < stop; j++) {
		source4[j] /= expinv4[j];
		y[j]=(y0[j]+b50*source0[j]+b51*source[j]+b52*source2[j]+b53*source3[j]+
			  b54*source4[j])*expinv5[j];
	}
	Source(source,y,t+a5);
	for(j=start; j < stop; j++) source[j] /= expinv5[j];
}

int I_RK5::Corrector(double dt, double& errmax, int start, int stop)
{
	RK5::Corrector(dt,errmax,start,stop);
	for(int j=start; j < stop; j++) y[j] *= expinv4[j];
	return 1;
}

inline void C_RK5::Correct(const Real y0, Real& y2, Real& y3,
						   const Real y4, Real& y,
						   const Real source0, const Real source2, 
						   const Real source3, const Real source4,
						   const Real source, const double, int& invertible)
{
	Real discr=y0*y0+2.0*(c0*y0*source0+c2*y2*source2+c3*y3*source3+
						  c5*y*source);
	if(discr < 0.0) invertible=0;
	else {
// Put discr in y2, pred in y3 for deferred error analysis
		y3=y0*y0+2.0*(d0*y0*source0+d2*y2*source2+d3*y3*source3+
					  d4*y4*source4+d5*y*source);
		y2=discr;
		y=sgn(y0+c0*source0+c2*source2+c3*source3+c5*source)*sqrt(discr);
	}
}

inline void C_RK5::Correct(const Complex y0, Complex& y2, Complex& y3,
						   const Complex y4, Complex& y,
						   const Complex source0, const Complex source2, 
						   const Complex source3, const Complex source4,
						   const Complex source, const double dt,
						   int& invertible)
{
	Correct(y0.re,y2.re,y3.re,y4.re,y.re,source0.re,source2.re,
			source3.re,source4.re,source.re,dt,invertible);
#ifndef _CRAY		
	if(!invertible) return;
#endif			
	Correct(y0.im,y2.im,y3.im,y4.im,y.im,source0.im,source2.im,
			source3.im,source4.im,source.im,dt,invertible);
}

int C_RK5::Corrector(double, double& errmax, int start, int stop)
{
	int j,invertible=1;
#pragma ivdep
	for(j=start; j < stop; j++) {
		Correct(y0[j],y2[j],y3[j],y4[j],y[j],source0[j],source2[j],
				source3[j],source4[j],source[j],dt,invertible);
#ifndef _CRAY		
		if(!invertible) return 0;
#endif			
	}
#if _CRAY
    if(!invertible) return 0;
#endif
	
	if(dynamic) {
		for(j=start; j < stop; j++) 
			calc_error(y0[j]*y0[j],y2[j],y3[j],y2[j],errmax);
		ExtrapolateTimestep(errmax);
	}
	return 1;
}

void display_invariants(Real E, Real Z, Real P)
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
	if(Npsi <= Nlimit || verbose > 1) {
		for(i=0; i < Npsi; i++) cout << "psi[" << i << "] = " << y[i] << newl;
		cout << endl;
	}
	
	compute_invariants(y,Npsi,E,Z,P);
	display_invariants(E,Z,P);
	
	if(Nmoment > 2 && t) {
		cout << newl << "AVERAGED VALUES:" << newl << endl;
// We overwrite y+Npsi here, since it is no longer needed.
		Var *y2=y+3*Npsi;
		for(i=0; i < Npsi; i++) y2[i] = (real(y2[i])+imag(y2[i]))/t;
		
		if(Npsi <= Nlimit || verbose > 1) {
			for(i=0; i < Npsi; i++) {
				Real y2avg=y2[i].re;
				cout << "|psi|^2 [" << i << "] = " << y2avg << newl;
			}
			cout << endl;
		}
		
		for(i=0; i < Npsi; i++) y2[i] = sqrt(y2[i].re);
		compute_invariants(y2,Npsi,E,Z,P);
		display_invariants(E,Z,P);
	}
}
