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
	
	// Compute moments
	if(average && Nmoment > 0) {
		Var *k, *q=source+Npsi, *kstop=psi+Npsi;
#pragma ivdep
		for(k=psi; k < kstop; k++, q++) *q=product(*k,*k);   // psi^2
		for(int n=1; n < Nmoment; n++) { // psi^n
#pragma ivdep
			for(k=psi; k < kstop; k++, q++) *q=product(*(q-Npsi),*k);
		}
	}
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
	
	// Compute moments
	if(average && Nmoment > 0) {
		Var *k, *q=source+Npsi, *kstop=psi+Npsi;
#pragma ivdep
		for(k=psi; k < kstop; k++, q++) *q=product(*k,*k);   // psi^2
		for(int n=1; n < Nmoment; n++) { // psi^n
#pragma ivdep
			for(k=psi; k < kstop; k++, q++) *q=product(*(q-Npsi),*k);
		}
	}
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
	
	// Compute moments
	if(average && Nmoment > 0) {
		Var *k, *q=source+Npsi, *kstop=psi+Npsi;
#pragma ivdep
		for(k=psi; k < kstop; k++, q++) *q=product(*k,*k);   // psi^2
		for(int n=1; n < Nmoment; n++) { // psi^n
#pragma ivdep
			for(k=psi; k < kstop; k++, q++) *q=product(*(q-Npsi),*k);
		}
	}
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

Solve_RC C_Euler::Solve(Real *y0, double t, double dt)
{
	int j,iter,cont,jfix=-1;
	Real *rsource=(Real *) source, *ry=(Real *) y;
	Var *vy0=(Var *) y0;
	Real yj;
	double tau,mu,temp;
	
	tau=dt;
	Source(source,vy0,t);
	
	mu=temp=0.0;
	for(j=0; j < nyconserve; j++) {
		if(y0[j] != 0.0) temp=tau*rsource[j]/y0[j];
		else jfix=j;
		if(fabs(temp) > fabs(mu)) {mu=temp;	if(fabs(mu) >= 0.5) jfix=j;}
	}
	mu += 2.0;
	
	if(jfix == -1)	// Evolve forwards.
		for(j=0; j < nyconserve; j++)
			y0[j]=sgn(y0[j])*sqrt(y0[j]*(y0[j]+mu*tau*rsource[j]));
	else {			// Iterate backwards.
		set(ry,y0,nyconserve);
		iter=0;
		for(j=0; j < nyconserve; j++) lastdiff[j]=0.0;
		y0[jfix]=ry[jfix]+tau*rsource[jfix];
		temp=(ry[jfix]/y0[jfix]+1.0)*rsource[jfix];
		do {
			cont=0;
			Source(source,vy0,t);
			mu=temp/rsource[jfix];
			for(j=0; j < nyconserve; j++) {
				yj=y0[j];
				if(j !=jfix) {
					Real discr=ry[j]*ry[j]+mu*tau*rsource[j]*y0[j];
					if(discr >= 0.0) y0[j]=sgn(ry[j])*sqrt(discr);
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
	for(j=nyconserve; j < ny; j++) y0[j] += dt*rsource[j];
	return SUCCESSFUL;
}

Solve_RC C_Euler::Solve(Complex *, double, double)
{
	msg(ERROR,"Complex C_Euler has not been implemented");
	return UNSUCCESSFUL;
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

int C_PC::Corrector(Var *y0, double dt, double& errmax, int start, int stop)
{
	int j;
	if(dynamic) for(j=start; j < stop; j++) {
		Var pred=y1[j];
		if(!Correct(y0[j],y1[j],y[j],source0[j],source[j],dt)) return 0;
		calc_error(y0[j],y[j],pred,y[j],errmax);
	}
	else for(j=start; j < stop; j++) {
		if(!Correct(y0[j],y1[j],y[j],source0[j],source[j],dt)) return 0;
	}
	return 1;
}

void E_PC::Allocate(int n)
{
	PC::Allocate(n);
	expinv=new Nu[Npsi];
	onemexpinv=new Nu[Npsi];
	
	nu_inv=new Nu[Npsi];

	for(int j=0; j < Npsi; j++) {
		if(nu[j] != 0.0) nu_inv[j]=1.0/nu[j];
		else nu_inv[j]=1.0;
	}
}

void E_PC::TimestepDependence(double dt)
{
	for(int j=0; j < Npsi; j++) {
		onemexpinv[j]=-expm1(-nu[j]*dt);
		expinv[j]=1.0-onemexpinv[j];
		if(nu[j] == 0.0) onemexpinv[j]=dt;
	}
}

void E_PC::Predictor(Var *y0, double t, double)
{
	for(int j=0; j < Npsi; j++)	y1[j]=expinv[j]*y0[j]+onemexpinv[j]*source0[j];
	for(int j=Npsi; j < ny; j++) y1[j]=y0[j]+dt*source0[j];
	Source(source,y1,t+dt);
}

int E_PC::Corrector(Var *y0, double, double& errmax, int start, int stop)
{
	int j;
	if(dynamic)	for(j=start; j < stop; j++) {
		Var corr=0.5*(source0[j]+source[j]);
		y[j]=expinv[j]*y0[j]+onemexpinv[j]*corr;
		calc_error(source0[j],corr,source0[j],corr,errmax);
	}
	else for(j=start; j < stop; j++) {
		y[j]=expinv[j]*y0[j]+onemexpinv[j]*0.5*(source0[j]+source[j]);
	}
	return 1;
}

void I_PC::Allocate(int n)
{
	PC::Allocate(n);
	expinv=new Nu[Npsi];
}

void I_PC::TimestepDependence(double dt)
{
	for(int j=0; j < Npsi; j++) expinv[j]=exp(-nu[j]*dt);
}

void I_PC::Predictor(Var *y0, double t, double dt)
{
	if(new_y0) for(int j=0; j < Npsi; j++) y0[j] *= expinv[j];
	for(int j=0; j < Npsi; j++) source0[j] *= expinv[j];
	PC::Predictor(y0,t,dt);
}

void I_RK2::Allocate(int n)
{
	RK2::Allocate(n);
	expinv=new Nu[Npsi];
}

void I_RK2::TimestepDependence(double dt)
{
	RK2::TimestepDependence(dt);
	for(int j=0; j < Npsi; j++) expinv[j]=exp(-nu[j]*halfdt);
}

void I_RK2::Predictor(Var *y0, double t, double dt)
{
	if(new_y0) for(int j=0; j < Npsi; j++) y0[j] *= expinv[j];
	else for(int j=0; j < Npsi; j++) y0[j] /= expinv[j];
	for(int j=0; j < Npsi; j++) source0[j] *= expinv[j];
	RK2::Predictor(y0,t,dt);
	if(new_y0) for(int j=0; j < Npsi; j++) y0[j] *= expinv[j];
}

void E_RK2::Allocate(int n)
{
	E_PC::Allocate(n);
	expinv1=new Nu[Npsi];
	onemexpinv1=new Nu[Npsi];
}

void E_RK2::TimestepDependence(double dt)
{
	for(int j=0; j < Npsi; j++) {
		onemexpinv[j]=-expm1(-nu[j]*dt);
		expinv[j]=1.0-onemexpinv[j];
		
		Nu factor=1.0/(expinv[j]+0.5*onemexpinv[j]);
		if(nu[j] == 0.0) onemexpinv[j]=dt;
		
		expinv1[j]=expinv[j]*factor;
		onemexpinv1[j]=0.5*onemexpinv[j]*factor;
	}
}

void E_RK2::Predictor(Var *y0, double t, double)
{
	for(int j=0; j < Npsi; j++)
		y1[j]=expinv1[j]*y0[j]+onemexpinv1[j]*source0[j];
	for(int j=Npsi; j < ny; j++) y1[j]=y0[j]+halfdt*source0[j];
	Source(source,y1,t+halfdt);
}

int E_RK2::Corrector(Var *y0, double, double& errmax, int start, int stop)
{
	int j;
	for(j=start; j < stop; j++) y[j]=expinv[j]*y0[j]+onemexpinv[j]*source[j];
	if(dynamic)	for(j=start; j < stop; j++) {
		calc_error(source0[j],source[j],source0[j],source[j],errmax);
	}
	return 1;
}

void E_RK4::Allocate(int n)
{
	E_RK2::Allocate(n);
	source1=new Var[n];
	source2=new Var[n];
}


void E_RK4::Predictor(Var *y0, double t, double dt)
{
	int j;
	
	for(j=0; j < Npsi; j++) y[j]=expinv1[j]*y0[j]+onemexpinv1[j]*source0[j];
	for(j=Npsi; j < ny; j++) y[j]=y0[j]+halfdt*source0[j];
	Source(source1,y,t+halfdt);
	
	for(j=0; j < Npsi; j++) y[j]=expinv1[j]*y0[j]+onemexpinv1[j]*source1[j];
	for(j=Npsi; j < ny; j++) y[j]=y0[j]+halfdt*source1[j];
	Source(source2,y,t+halfdt);
	
	for(j=0; j < Npsi; j++) y[j]=expinv[j]*y0[j]+onemexpinv[j]*source2[j];
	for(j=Npsi; j < ny; j++) y[j]=y0[j]+dt*source2[j];
	Source(source,y,t+dt);
}

int E_RK4::Corrector(Var *y0, double, double& errmax, int start, int stop)
{
	const double sixth=1.0/6.0;
	int j;
	
	if(dynamic)	for(j=start; j < stop; j++) {
		Var corr=sixth*(source0[j]+2.0*source1[j]+2.0*source2[j]+source[j]);
		y[j]=expinv[j]*y0[j]+onemexpinv[j]*corr;
		calc_error(source0[j],corr,source2[j],corr,errmax);
	}
	else for(j=start; j < stop; j++) {
		y[j]=expinv[j]*y0[j]+onemexpinv[j]*sixth*(source0[j]+2.0*source1[j]+
												  2.0*source2[j]+source[j]);
	}
	return 1;
}

void CE_PC::Allocate(int n)
{
	PC::Allocate(n);
	expinv=new Nu[n];
	onemexpinv=new Real[n];
	
	nuR_inv=new Real[Npsi];
	nuI=new Real[Npsi];

	for(int j=0; j < Npsi; j++) {
		nuI[j]=imag(nu[j]);
		if(real(nu[j]) != 0.0) nuR_inv[j]=1.0/real(nu[j]);
		else nuR_inv[j]=1.0;
	}
}

void CE_PC::TimestepDependence(double dt)
{
	int j;
	for(j=0; j < Npsi; j++) {
		onemexpinv[j]=-expm1(-real(nu[j])*dt);
		expinv[j]=exp(-nu[j]*dt);
		
		if(real(nu[j]) == 0.0) onemexpinv[j]=dt;
	}
	for(j=Npsi; j < ny; j++) {
		expinv[j]=1.0; onemexpinv[j]=dt;
	}
}

void CE_PC::Predictor(Var *y0, double t, double)
{
	for(int j=0; j < ny; j++)
		y[j]=expinv[j]*y0[j]+onemexpinv[j]*source0[j];
	Source(source,y,t+dt);
}

int CE_PC::Corrector(Var *y0, double, double& errmax, int start, int stop)
{
	int j;
	if(dynamic)	for(j=start; j < stop; j++) {
		Var y0j=expinv[j]*y0[j];
		Var corr=onemexpinv[j]*source[j];
		Var corr0=onemexpinv[j]*source0[j];
		if(!Correct(y0j,y1[j],y[j],corr0,corr,1.0)) return 0;
		calc_error(y0j,corr,corr0,corr,errmax);
	}
	else for(j=start; j < stop; j++) {
		if(!Correct(expinv[j]*y0[j],y1[j],y[j],source0[j],source[j],
					onemexpinv[j]))
			return 0;
	}
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

int C_RK2::Corrector(Var *y0, double dt, double& errmax, int start, int stop)
{
	int j;
	if(dynamic) for(j=start; j < stop; j++) {
		Var pred=y0[j]+dt*source0[j];
		if(!Correct(y0[j],y1[j],y[j],source0[j],source[j],dt)) return 0;
		calc_error(y0[j],y[j],pred,y[j],errmax);
	}
	else for(j=start; j < stop; j++) {
		if(!Correct(y0[j],y1[j],y[j],source0[j],source[j],dt)) return 0;
	}
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

int C_RK4::Corrector(Var *y0, double dt, double& errmax, int start, int stop)
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

int C_RK5::Corrector(Var *y0, double, double& errmax, int start, int stop)
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

void E_RK5::Allocate(int n)
{
	RK5::Allocate(n);
	nu_inv=new Nu[Npsi];

	for(int j=0; j < Npsi; j++) {
		if(nu[j] != 0.0) nu_inv[j]=1.0/nu[j];
		else nu_inv[j]=1.0;
	}
	
	expinv1=new Nu[Npsi];
	expinv2=new Nu[Npsi];
	expinv3=new Nu[Npsi];
	expinv4=new Nu[Npsi];
	expinv5=new Nu[Npsi];
	expinv=expinv4;
	onemexpinv1=new Nu[Npsi];
	onemexpinv2=new Nu[Npsi];
	onemexpinv3=new Nu[Npsi];
	onemexpinv4=new Nu[Npsi];
	onemexpinv5=new Nu[Npsi];
	onemexpinv=onemexpinv4;
}


inline void E_RK5::CalcExp(Real *expinvi, Real *onemexpinvi, double a) {
	for(int j=0; j < Npsi; j++) {
		Nu factor=1.0/(expinv[j]+a*onemexpinv[j]);
		expinvi[j]=expinv[j]*factor;
		onemexpinvi[j]=onemexpinv[j]*factor;
		if(nu[j] == 0.0) onemexpinvi[j]=dt;
	}
}

void E_RK5::TimestepDependence(double dt)
{
	RK5::TimestepDependence(dt);
	
	for(int j=0; j < Npsi; j++) {
		onemexpinv[j]=-expm1(-nu[j]*dt);
		expinv[j]=1.0-onemexpinv[j];
	}
		
	CalcExp(expinv1,onemexpinv1,0.2);
	CalcExp(expinv2,onemexpinv2,0.3);
	CalcExp(expinv3,onemexpinv3,0.6);
	CalcExp(expinv5,onemexpinv5,0.875);
	
	for(int j=0; j < Npsi; j++) {
		if(nu[j] == 0.0) onemexpinv[j]=dt;
	}
}

void E_RK5::Predictor(Var *y0, double t, double)
{
	int j;
	const double B10=0.2;
	const double B20=3.0/40.0, B21=9.0/40.0;
	const double B30=0.3, B31=-0.9, B32=1.2;
	const double B40=-11.0/54.0, B41=2.5, B42=-70.0/27.0, B43=35.0/27.0;
	const double B50=1631.0/55296.0, B51=175.0/512.0, B52=575.0/13824.0;
	const double B53=44275.0/110592.0, B54=253.0/4096.0;
	
#pragma ivdep		
	for(j=0; j < Npsi; j++)
		y[j]=(expinv[j]*y0[j]+(nu[j] ? onemexpinv[j]: dt)*
		B10*source0[j])/(expinv[j]+a*onemexpinv[j]);
	for(j=Npsi; j < ny; j++) y[j]=y0[j]+b10*source0[j];
	Source(source,y,t+a1);
#pragma ivdep		
	for(j=0; j < Npsi; j++)
		y2[j]=expinv2[j]*y0[j]+onemexpinv2[j]*(B20*source0[j]+B21*source[j]);
	for(j=Npsi; j < ny; j++) y2[j]=y0[j]+b20*source0[j]+b21*source[j];
	Source(source2,y2,t+a2);
#pragma ivdep		
	for(j=0; j < Npsi; j++)
		y3[j]=expinv3[j]*y0[j]+onemexpinv3[j]*(B30*source0[j]+B31*source[j]
											   +B32*source2[j]);
	for(j=Npsi; j < ny; j++) y3[j]=y0[j]+b30*source0[j]+b31*source[j]+
		b32*source2[j];
	Source(source3,y3,t+a3);
#pragma ivdep		
	for(j=0; j < Npsi; j++) 
		y4[j]=expinv4[j]*y0[j]+onemexpinv4[j]*(B40*source0[j]+B41*source[j]+
											   B42*source2[j]+B43*source3[j]);
	for(j=Npsi; j < ny; j++) y4[j]=y0[j]+b40*source0[j]+b41*source[j]+
		b42*source2[j]+b43*source3[j];
	Source(source4,y4,t+a4);
#pragma ivdep		
	for(j=0; j < Npsi; j++) 
		y[j]=expinv5[j]*y0[j]+onemexpinv5[j]*(B50*source0[j]+B51*source[j]+
											  B52*source2[j]+B53*source3[j]+
											  B54*source4[j]);
	for(j=Npsi; j < ny; j++) 
		y[j]=y0[j]+b50*source0[j]+b51*source[j]+b52*source2[j]+b53*source3[j]+
			b54*source4[j];
	Source(source,y,t+a5);
}


inline void E_RK5::Correct(const Real y0, Real& y,
						   const Real expinv, const Real onemexpinv,
						   const Real source0, const Real source2, 
						   const Real source3, const Real,
						   const Real source, const double)
{
	const double c0=37.0/378.0, c2=250.0/621.0, c3=125.0/594.0, c5=512.0/1771.0;
	y=expinv*y0+onemexpinv*(c0*source0+c2*source2+c3*source3+c5*source);
}

inline void E_RK5::CalcError(const Real y0, Real& y,
							 const Real source0, const Real source2, 
							 const Real source3, const Real source4,
							 const Real source, const double)
{
	y=y0+d0*source0+d2*source2+d3*source3+d4*source4+d5*source;
}

inline void E_RK5::Correct(const Complex y0, Complex& y,
						   const Real expinv, const Real onemexpinv,
						   const Complex source0, const Complex source2, 
						   const Complex source3, const Complex source4,
						   const Complex source, const double)
{
	Correct(y0.re,y.re,expinv,onemexpinv,
			source0.re,source2.re,source3.re,source4.re,source.re,dt);
	Correct(y0.im,y.im,expinv,onemexpinv,
			source0.im,source2.im,source3.im,source4.im,source.im,dt);
}

inline void E_RK5::CalcError(const Complex y0, Complex& y,
						   const Complex source0, const Complex source2, 
						   const Complex source3, const Complex source4,
						   const Complex source, const double)
{
	CalcError(y0.re,y.re,source0.re,source2.re,source3.re,source4.re,
			  source.re,dt);
	CalcError(y0.im,y.im,source0.im,source2.im,source3.im,source4.im,
			  source.im,dt);
}

int E_RK5::Corrector(Var *y0, double, double& errmax, int start, int stop)
{
	int j;
	Var pred;

#pragma ivdep		
	for(j=start; j < stop; j++)
		Correct(y0[j],y[j],expinv[j],onemexpinv[j],source0[j],source2[j],
				source3[j],source4[j],source[j],dt);
	if(dynamic) {
#pragma ivdep		
		for(j=start; j < stop; j++)
			CalcError(y0[j],y3[j],source0[j],source2[j],source3[j],source4[j],
					  source[j],dt);
		for(j=start; j < stop; j++)
			calc_error(y0[j],y[j],y3[j],y[j],errmax);
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
	
	if(average && t) {
		cout << newl << "AVERAGED VALUES:" << newl << endl;
// We overwrite y+Npsi here, since it is no longer needed.
		Var *y2=y+Npsi;
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

