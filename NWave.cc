#include "NWave.h"
#include "Cartesian.h"

int Npsi;
int NpsiR;
int Ntotal;

// (0 => evolve all modes, 1 => evolve only half of the modes).
int reality=1; // Reality condition flag 

Var *psibuffer,*psibuffer0,*psibufferR,*psibufferStop,*pqbuffer;
Var *convolution,*convolution0;
int pseudospectral=0;
unsigned int log2n; // Number of FFT levels

Var **pqIndex;
int *qStart;

int Ntriad;
DynVector<Triad> triad;
TriadLimits *triadLimits;

int Nchainp,Nchainn;
DynVector<Chain> chainp,chainn;
Chain *chainpBase,*chainnBase;
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
		Var *k, *q=source+Npsi, *kstop=psibuffer+Npsi;
#pragma ivdep
		for(k=psibuffer; k < kstop; k++, q++) *q=product(*k,*k);   // psi^2
		for(int n=1; n < Nmoment; n++) { // psi^n
#pragma ivdep
			for(k=psibuffer; k < kstop; k++, q++) *q=product(*(q-Npsi),*k);
		}
	}
}

void PrimitiveNonlinearity(Var *source, Var *psi, double)
{
	int i;
	Cartesian mk,mp;
	extern Cartesian *CartesianMode;

	set(psibuffer,psi,Npsi);
	
	// Compute reflected psi's
#pragma ivdep		
	for(Var *k=psibuffer; k < psibufferR; k++) conjugate(*(k+Npsi),*k);
	
	for(i=0; i < Npsi; i++) source[i]=0.0;
		
	for(i=0; i < Nchainp; i++) {
		Chain *c=chainpBase+i;
		mk=CartesianMode[c->k]; mp=CartesianMode[c->p];
		Real kx=mk.Kx(), ky=mk.Ky(), py=mp.Ky();
		Real kxpy=kx*py, py2=py*py;
		Var sum=0.0, *psip=psibuffer+c->p, *psiq=c->psiq;
		int stop=c->stop;
		for(int px=mp.Column(); px < stop; px++,psiq++)
			sum += (ky*px-kxpy)*(px*px+py2)*conj(psip[px]*(*psiq));
		source[c->k] += sum;
	}
	
	for(i=0; i < Nchainn; i++) {
		Chain *c=chainnBase+i;
		mk=CartesianMode[c->k]; mp=CartesianMode[c->p];
		Real kx=mk.Kx(), ky=mk.Ky(), py=mp.Ky();
		Real kxpy=kx*py, py2=py*py;
		Var sum=0.0, *psip=psibuffer+c->p, *psiq=c->psiq;
		int stop=c->stop;
		for(int px=mp.Column(); px < stop; px++,psiq--)
			sum += (ky*px-kxpy)*(px*px+py2)*conj(psip[px]*(*psiq));
		source[c->k] += sum;
	}
	
	for(i=0; i < Npsi; i++) source[i] *= kinv2[i];

	// Compute moments
	if(average && Nmoment > 0) {
		Var *k, *q=source+Npsi, *kstop=psibuffer+Npsi;
#pragma ivdep
		for(k=psibuffer; k < kstop; k++, q++) *q=product(*k,*k);   // psi^2
		for(int n=1; n < Nmoment; n++) { // psi^n
#pragma ivdep
			for(k=psibuffer; k < kstop; k++, q++) *q=product(*(q-Npsi),*k);
		}
	}
}

	
void convolve(Complex *H, Complex *F, Complex *G, unsigned int m, unsigned
			  int log2n);
void convolve0(Complex *H, Complex *F, Complex *g, unsigned int m, unsigned
			   int log2n);
	
void PrimitiveNonlinearityFFT(Var *source, Var *psi, double)
{
	extern Cartesian *CartesianMode;
	int i;

	*psibuffer0=0.0;
	*convolution0=0.0;
	set(psibuffer,psi,Npsi);
	set(convolution,psi,Npsi);
	for(i=0; i < Npsi; i++)
		convolution[i] *= CartesianMode[i].K2()*CartesianMode[i].Ky();
	convolve(convolution0,convolution0,psibuffer0,Npsi+1,log2n);
	set(source,convolution,Npsi);
							
	*convolution0=0.0;
	set(convolution,psi,Npsi);
	for(i=0; i < Npsi; i++)
		convolution[i] *= CartesianMode[i].K2()*CartesianMode[i].Kx();
	// Reuse FFT of psibuffer0.
	convolve0(convolution0,convolution0,psibuffer0,Npsi+1,log2n);
		
	for(i=0; i < Npsi; i++)
		source[i] = kinv2[i]*(CartesianMode[i].Kx()*source[i]-
							  CartesianMode[i].Ky()*convolution[i]);
		
	// Compute moments
	if(average && Nmoment > 0) {
		Var *k, *q=source+Npsi, *kstop=psibuffer+Npsi;
#pragma ivdep
		for(k=psibuffer; k < kstop; k++, q++) *q=product(*k,*k);   // psi^2
		for(int n=1; n < Nmoment; n++) { // psi^n
#pragma ivdep
			for(k=psibuffer; k < kstop; k++, q++) *q=product(*(q-Npsi),*k);
		}
	}
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
	ny=n;
	PC::Allocate(n);
	expinv=new Nu[n];
	onemexpinv=new Nu[n];
	
	nu_inv=new Nu[Npsi];

	for(int j=0; j < Npsi; j++) {
		if(nu[j] != 0.0) nu_inv[j]=1.0/nu[j];
		else nu_inv[j]=1.0;
	}
}

void E_PC::TimestepDependence(double dt)
{
    int j;
	for(j=0; j < Npsi; j++) {
		onemexpinv[j]=-expm1(-nu[j]*dt);
		expinv[j]=1.0-onemexpinv[j];
		
		if(nu[j] == 0.0) onemexpinv[j]=dt;
	}
	for(j=Npsi; j < ny; j++) {
		expinv[j]=1.0; onemexpinv[j]=dt;
	}
}

void E_PC::Predictor(Var *y0, double t, double)
{
	for(int j=0; j < ny; j++)
		y[j]=expinv[j]*y0[j]+onemexpinv[j]*source0[j];
	Source(source,y,t+dt);
}

int E_PC::Corrector(Var *y0, double, double& errmax, int start, int stop)
{
	Var y0j,corr;
	int j;
	if(dynamic)	for(j=start; j < stop; j++) {
		y0j=expinv[j]*y0[j];
		corr=onemexpinv[j]*0.5*(source0[j]+source[j]);
		y[j]=y0j+corr;
		calc_error(y0j,corr,onemexpinv[j]*source0[j],corr,errmax);
	}
	else for(j=start; j < stop; j++) {
		y[j]=expinv[j]*y0[j]+onemexpinv[j]*0.5*(source0[j]+source[j]);
	}
	return 1;
}

void E_RK2::Allocate(int n)
{
	E_PC::Allocate(n);
	expinv1=new Nu[n];
	onemexpinv1=new Nu[n];
}

void E_RK2::TimestepDependence(double dt)
{
	Nu factor;
	int j;
	
	halfdt=0.5*dt;
	for(j=0; j < Npsi; j++) {
		onemexpinv[j]=-expm1(-nu[j]*dt);
		expinv[j]=1.0-onemexpinv[j];
		
		factor=1.0/(2.0-onemexpinv[j]);
		if(nu[j] == 0.0) onemexpinv[j]=dt;
		
		expinv1[j]=2.0*expinv[j]*factor;
		onemexpinv1[j]=onemexpinv[j]*factor;
	}
	
	for(j=Npsi; j < ny; j++) {
		expinv1[j]=expinv[j]=1.0;
		onemexpinv1[j]=onemexpinv[j]=dt;
	}
}

void E_RK2::Predictor(Var *y0, double t, double)
{
	for(int j=0; j < ny; j++)
		y[j]=expinv1[j]*y0[j]+onemexpinv1[j]*source0[j];
	Source(source,y,t+halfdt);
}

int E_RK2::Corrector(Var *y0, double, double& errmax, int start, int stop)
{
	Var y0j,corr;
	int j;
	if(dynamic)	for(j=start; j < stop; j++) {
		y0j=expinv[j]*y0[j];
		corr=onemexpinv[j]*source[j];
		y[j]=y0j+corr;
		calc_error(y0j,corr,onemexpinv[j]*source0[j],corr,errmax);
	}
	else for(j=start; j < stop; j++) {
		y[j]=expinv[j]*y0[j]+onemexpinv[j]*source[j];
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
	
	for(j=0; j < ny; j++)
		y[j]=expinv1[j]*y0[j]+onemexpinv1[j]*source0[j];
	Source(source1,y,t+halfdt);
	
	for(j=0; j < ny; j++) 
		y[j]=expinv1[j]*y0[j]+onemexpinv1[j]*source1[j];
	Source(source2,y,t+halfdt);
	
	for(j=0; j < ny; j++)
		y[j]=expinv[j]*y0[j]+onemexpinv[j]*source2[j];
	Source(source,y,t+dt);
}

int E_RK4::Corrector(Var *y0, double, double& errmax, int start, int stop)
{
	Var y0j,corr;
	const double sixth=1.0/6.0;
	int j;
	
	if(dynamic)	for(j=start; j < stop; j++) {
		y0j=expinv[j]*y0[j];
		corr=onemexpinv[j]*sixth*
			(source0[j]+2.0*source1[j]+2.0*source2[j]+source[j]);
		y[j]=y0j+corr;
		calc_error(y0j,corr,onemexpinv[j]*source2[j],corr,errmax);
	}
	else for(j=start; j < stop; j++) {
		y[j]=expinv[j]*y0[j]+onemexpinv[j]*sixth*(source0[j]+2.0*source1[j]+
												  2.0*source2[j]+source[j]);
	}
	return 1;
}

void CE_PC::Allocate(int n)
{
	ny=n;
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


inline int C_RK5::Correct(const Real y0, const Real y2, const Real y3,
						  const Real y4, const Real y5, Real& y,
						  const Real source0, const Real source2, 
						  const Real source3, const Real source4,
						  const Real source, const double,
						  double& errmax)
{
	Real discr=y0*y0+2.0*(c0*y0*source0+c2*y2*source2+c3*y3*source3+
						  c5*y5*source);
	if(discr < 0.0) return 0;
	Var pred=y0*y0+2.0*(d0*y0*source0+d2*y2*source2+d3*y3*source3+
						d4*y4*source4+d5*y5*source);
	y=sgn(y0+c0*source0+c2*source2+c3*source3+c5*source)*sqrt(discr);
	calc_error(y0*y0,discr,pred,discr,errmax);
	return 1;
}

inline int C_RK5::Correct(const Complex y0, const Complex y2, const Complex y3,
						  const Complex y4, const Complex y5, Complex& y,
						  const Complex source0, const Complex source2, 
						  const Complex source3, const Complex source4,
						  const Complex source, const double dt,
						  double& errmax)
{
	if(!Correct(y0.re,y2.re,y3.re,y4.re,y5.re,y.re,source0.re,source2.re,
				source3.re,source4.re,source.re,dt,errmax)) return 0;
	if(!Correct(y0.im,y2.im,y3.im,y4.im,y5.im,y.im,source0.im,source2.im,
				source3.im,source4.im,source.im,dt,errmax)) return 0;
	return 1;
}

int C_RK5::Corrector(Var *y0, double, double& errmax, int start, int stop)
{
	for(int j=start; j < stop; j++) {
		if(!Correct(y0[j],y2[j],y3[j],y4[j],y5[j],y[j],source0[j],
					source2[j],source3[j],source4[j],source[j],dt,errmax))
			return 0;
	}
	if(dynamic) ExtrapolateTimestep(errmax);
	return 1;
}

void compute_invariants(Var *y, int Npsi, Real& E, Real& Z, Real& P)
{
	Real Ek,Zk,Pk,k2;
	E=Z=P=0.0;
	for(int i=0; i < Npsi; i++) {
		k2=Geometry->K2(i);
		Ek=Geometry->Normalization(i)*abs2(y[i])*Geometry->Area(i);
		Zk=k2*Ek;
		Pk=k2*Zk;
		E += Ek;
		Z += Zk;
		P += Pk;
	}
	
	Real factor=(reality ? 1.0: 0.5)*continuum_factor;
	E *= factor;
	Z *= factor;
	P *= factor;
}	

void display_invariants(Real E, Real Z, Real P)
{
	cout << "Energy = " << E << endl;
	cout << "Enstrophy = " << Z << endl;
	cout << "Palinstrophy = " << P << endl;
}

void NWave::FinalOutput()
{
	Real E,Z,P;
	int i;
	
	cout << endl << "FINAL VALUES:" << endl << endl;
	for(i=0; i < Npsi; i++) cout << "psi[" << i << "] = " << y[i] << endl;
	cout << endl;
	
	compute_invariants(y,Npsi,E,Z,P);
	display_invariants(E,Z,P);
	
	if(average && t) {
		cout << endl << "AVERAGED VALUES:" << endl << endl;
// We overwrite part of y here, since it is no longer needed.
		Var *y2=y+Npsi;
		for(i=0; i < Npsi; i++) {
			Real y2avg = (real(y2[i])+imag(y2[i]))/t;
			y2[i] = sqrt(y2avg);
			cout << "|psi|^2 [" << i << "] = " << y2avg << endl;
		}
		cout << endl;
		compute_invariants(y2,Npsi,E,Z,P);
		display_invariants(E,Z,P);
	}
}
