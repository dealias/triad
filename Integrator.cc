#include "options.h"
#include "kernel.h"

inline void IntegratorBase::ChangeTimestep(double& dt, const double dtnew,
										   const double t)
{
	if(verbose > 1) cout << newl << "Time step changed from " << dt <<
						" to " << dtnew << " at t=" << t << "." << endl;
	if(dtnew == 0.0) msg(ERROR,"Zero time step encountered");
	TimestepDependence(dt=dtnew);
}

static clock_t realtime,lasttime=0;
static const int nperline=10;

void IntegratorBase::Integrate(Var *const y, double& t, double tmax,
							   double& dt, const double sample)
	// Don't dump or microprocess if sample is negative.
{
	double dtold=0.0, dtorig=0.0;
	int it,itx,cont;
	int nout=0;
	const int forwards=(tmax >= t);
	const double sign=(forwards ? 1.0 : -1.0);
	
	y0=y;
	
	double tstop=((sample > 0.0) ? 0.0 : tmax);
	dt=min(dt/Microfactor(),dtmax);
	dt *= sign;
	if(dt == 0.0) msg(ERROR,"Zero time step encountered");
	
	TimestepDependence(dt);
	microprocess=(sample >= 0.0) ? Problem->Microprocess() : 0;
	
	// Main integration loop
	for(it=0; it < itmax && (forwards ? t < tmax : t > tmax); it++) {
		
		if(verbose) {
			if((it % nperline) == 0) cout << newl;
			cout << "[" << it << flush;
		}
		
		if(sample > 0.0 && (forwards ? t >= tstop : t <= tstop)) {
			dump(it,0,tmax);
			tstop=sign*min(++nout*sample,sign*tmax);
			if(dtorig) {ChangeTimestep(dt,dtorig,t); dtorig=0.0;}
		}
		else if(sample == 0.0) dump(it,0,tmax);
		
		for(itx=0; itx < microsteps; itx++) {
			if(microprocess) Problem->Microprocess();
			if(polltime) {
				realtime=time(NULL);
				if(realtime-lasttime > polltime) {
					if (poll()) {tmax=tstop=t; exit_signal=CONTINUE;}
					lasttime=realtime;
				}
			}
			
			if(forwards ? t+dt > tstop : t+dt < tstop) {
				if(t==tstop) break;
				dtorig=dt;
				ChangeTimestep(dt,tstop-t,t);
				itx=microsteps-1; // This is the final iteration.
			}
			cont=1;
			do {
				switch(Solve(t,dt))	{
				case ADJUST:
					t += dt;
					ChangeTimestep(dt,sign*min(sign*dt*stepfactor,dtmax),t);
					cont=0;	break;
				case SUCCESSFUL:
					t += dt;
					if(dtold) {ChangeTimestep(dt,dtold,t); dtold=0.0;}   
					cont=0;	break;
				case NONINVERTIBLE:
					if(!dtold) dtold=dtorig ? dtorig : dt;
					invert_cnt++;
					ChangeTimestep(dt,dt*stepnoninverse,t); break;
				case UNSUCCESSFUL:
					ChangeTimestep(dt,dt*stepinverse,t);
				}
			} while(cont);
			iteration++;
		}
		if(verbose) cout << "] ";
	}
	
	if(verbose) cout << endl;
	if(dtorig) ChangeTimestep(dt,dtorig,t);
	if(sample >= 0.0) dump(it,1,tmax);
	dt *= sign;
}	

void IntegratorBase::Allocate(int n)
{
	ny=n; yi=NULL; source=new(0) Var[n];
	nyprimary=ny/(Nmoment+1);
	if(nyprimary*(Nmoment+1) != ny) 
		msg(ERROR, "ny=%d is incompatible with Nmoment=%d",ny,Nmoment);
}

Solve_RC Euler::Solve(double t, double dt)
{
	Source(source,y0,t);
	Problem->Transform(y0,t,dt,yi);
	for(int j=0; j < ny; j++) y0[j] += dt*source[j];
	Problem->BackTransform(y0,t+dt,dt,yi);
	return SUCCESSFUL;
}

Solve_RC PC::Solve(double t, double dt)
{
	double errmax=0.0;
	
	if(new_y0) Source(source0,y0,t);
	Problem->Transform(y0,t,dt,yi);
	
	Predictor(t,dt,0,nyprimary);
	if(!Corrector(dt,errmax,0,nyprimary)) {
		if(hybrid) StandardCorrector(dt,errmax,0,nyprimary);
		else return NONINVERTIBLE;
	}
	
	// Disregard averaging error
	if(Nmoment) {
		StandardPredictor(t,dt,nyprimary,ny);
		int dynamic_value=dynamic;
		dynamic=0;
		StandardCorrector(dt,errmax,nyprimary,ny);
		dynamic=dynamic_value;
	}

	Solve_RC flag=(dynamic ? CheckError(errmax) : SUCCESSFUL);
	new_y0=(flag != UNSUCCESSFUL);
	if(new_y0) {
		set(y0,y,ny);
		Problem->BackTransform(y0,t+dt,dt,yi);
	} else if(yi) set(y0,yi,ny);
	return flag;
}

void PC::Predictor(double t, double dt, int start, int stop)
{
	for(int j=start; j < stop; j++) y1[j]=y0[j]+dt*source0[j];
	Problem->BackTransform(y1,t+dt,dt,y0);
	Source(source,y1,t+dt);
}

int PC::Corrector(double dt, double& errmax, int start, int stop)
{
	const double halfdt=0.5*dt;
	for(int j=start; j < stop; j++) 
		y[j]=y0[j]+halfdt*(source0[j]+source[j]);
	if(dynamic) for(int j=start; j < stop; j++) 
		calc_error(y0[j],y[j],y0[j]+dt*source0[j],y[j],errmax);
	return 1;
}

void LeapFrog::Predictor(double t, double, int start, int stop)
{
	if(new_y0) {oldhalfdt=lasthalfdt;}
	else set(yp,yp0,ny);
	double dtprime=halfdt+oldhalfdt;
	Problem->Transform(yp,t-oldhalfdt,dtprime,yp0);
	for(int j=start; j < stop; j++) yp[j] += dtprime*source0[j];
	Problem->BackTransform(yp,t+halfdt,dtprime,yp0);
	Source(source,yp,t+halfdt);
	lasthalfdt=halfdt;
}
	
int LeapFrog::Corrector(double dt, double& errmax, int start, int stop)
{
	for(int j=start; j < stop; j++) y[j]=y0[j]+dt*source[j];
	if(dynamic) for(int j=start; j < stop; j++)
		calc_error(y0[j],y[j],y0[j]+dt*source0[j],y[j],errmax);
	return 1;
}

void RK2::TimestepDependence(double dt)
{
	halfdt=0.5*dt;
}

void RK2::Predictor(double t, double, int start, int stop)
{
	for(int j=start; j < stop; j++) y1[j]=y0[j]+halfdt*source0[j];
	Problem->BackTransform(y1,t+halfdt,halfdt,y0);
	Source(source,y1,t+halfdt);
}

int RK2::Corrector(double dt, double& errmax, int start, int stop)
{
	int j;
	for(j=start; j < stop; j++) y[j]=y0[j]+dt*source[j];
	if(dynamic) for(j=start; j < stop; j++)
		calc_error(y0[j],y[j],y0[j]+dt*source0[j],y[j],errmax);
	return 1;
}

void RK4::TimestepDependence(double dt)
{
	halfdt=0.5*dt;
	sixthdt=dt/6.0;
}

void RK4::Predictor(double t, double dt, int start, int stop)
{
	int j;
	for(j=start; j < stop; j++) y[j]=y0[j]+halfdt*source0[j];
	Problem->BackTransform(y,t+halfdt,halfdt,y0);
	Source(source1,y,t+halfdt);
	for(j=start; j < stop; j++) y[j]=y0[j]+halfdt*source1[j];
	Problem->BackTransform(y,t+halfdt,halfdt,y0);
	Source(source2,y,t+halfdt);
	for(j=start; j < stop; j++) y[j]=y0[j]+dt*source2[j];
	Problem->BackTransform(y,t+dt,dt,y0);
	Source(source,y,t+dt);
}

int RK4::Corrector(double, double& errmax, int start, int stop)
{
	int j;
	for(j=start; j < stop; j++) 
		y[j]=y0[j]+sixthdt*(source0[j]+2.0*(source1[j]+source2[j])+source[j]);
	if(dynamic) for(j=start; j < stop; j++)
		calc_error(y0[j],y[j],y0[j]+dt*source2[j],y[j],errmax);
	return 1;
}

void RK5::TimestepDependence(double dt)
{
	a1=0.2*dt; a2=0.3*dt; a3=0.6*dt; a4=dt; a5=0.875*dt;
	b10=0.2*dt;
	b20=3.0/40.0*dt; b21=9.0/40.0*dt;
	b30=0.3*dt; b31=-0.9*dt; b32=1.2*dt;
	b40=-11.0/54.0*dt; b41=2.5*dt; b42=-70.0/27.0*dt; b43=35.0/27.0*dt;
	b50=1631.0/55296.0*dt; b51=175.0/512.0*dt; b52=575.0/13824.0*dt;
	b53=44275.0/110592.0*dt; b54=253.0/4096.0*dt;
	c0=37.0/378.0*dt; c2=250.0/621.0*dt; c3=125.0/594.0*dt;	c5=512.0/1771.0*dt;
	d0=2825.0/27648.0*dt; d2=18575.0/48384.0*dt; d3=13525.0/55296.0*dt;
	d4=277.0/14336.0*dt; d5=0.25*dt;
}

#if COMPLEX

void RK5::Predictor(double t, double, int start, int stop)
{
	int j;
#pragma ivdep		
	for(j=start; j < stop; j++) {
		y[j].re=y0[j].re+b10*source0[j].re;
		y[j].im=y0[j].im+b10*source0[j].im;
	}
	Problem->BackTransform(y,t+a1,a1,y0);
	Source(source,y,t+a1);
#pragma ivdep		
	for(j=start; j < stop; j++) {
		y2[j].re=y0[j].re+b20*source0[j].re+b21*source[j].re;
		y2[j].im=y0[j].im+b20*source0[j].im+b21*source[j].im;
	}
	Problem->BackTransform(y2,t+a2,a2,y0);
	Source(source2,y2,t+a2);
#pragma ivdep		
	for(j=start; j < stop; j++) {
		y3[j].re=y0[j].re+b30*source0[j].re+b31*source[j].re+
			b32*source2[j].re;
		y3[j].im=y0[j].im+b30*source0[j].im+b31*source[j].im+
			b32*source2[j].im;
	}
	Problem->BackTransform(y3,t+a3,a3,y0);
	Source(source3,y3,t+a3);
#pragma ivdep		
	for(j=start; j < stop; j++) {
		y4[j].re=y0[j].re+b40*source0[j].re+b41*source[j].re+
			b42*source2[j].re+b43*source3[j].re;
		y4[j].im=y0[j].im+b40*source0[j].im+b41*source[j].im+
			b42*source2[j].im+b43*source3[j].im;
	}
	Problem->BackTransform(y4,t+a4,a4,y0);
	Source(source4,y4,t+a4);
#pragma ivdep		
	for(j=start; j < stop; j++) {
		y[j].re=y0[j].re+b50*source0[j].re+b51*source[j].re+
			b52*source2[j].re+b53*source3[j].re+b54*source4[j].re;
		y[j].im=y0[j].im+b50*source0[j].im+b51*source[j].im+
			b52*source2[j].im+b53*source3[j].im+b54*source4[j].im;
	}
	Problem->BackTransform(y,t+a5,a5,y0);
	Source(source,y,t+a5);
}
	
#else // COMPLEX
	
void RK5::Predictor(double t, double, int start, int stop)
{
	int j;
#pragma ivdep		
	for(j=start; j < stop; j++) y[j]=y0[j]+b10*source0[j];
	Problem->BackTransform(y,t+a1,a1,y0);
	Source(source,y,t+a1);
#pragma ivdep		
	for(j=start; j < stop; j++) y2[j]=y0[j]+b20*source0[j]+b21*source[j];
	Problem->BackTransform(y2,t+a2,a2,y0);
	Source(source2,y2,t+a2);
#pragma ivdep		
	for(j=start; j < stop; j++) y3[j]=y0[j]+b30*source0[j]+b31*source[j]+
									b32*source2[j];
	Problem->BackTransform(y3,t+a3,a3,y0);
	Source(source3,y3,t+a3);
#pragma ivdep		
	for(j=start; j < stop; j++) y4[j]=y0[j]+b40*source0[j]+b41*source[j]+
									b42*source2[j]+b43*source3[j];
	Problem->BackTransform(y4,t+a4,a4,y0);
	Source(source4,y4,t+a4);
#pragma ivdep		
	for(j=start; j < stop; j++) y[j]=y0[j]+b50*source0[j]+b51*source[j]+
									b52*source2[j]+b53*source3[j]+
									b54*source4[j];
	Problem->BackTransform(y,t+a5,a5,y0);
	Source(source,y,t+a5);
}

#endif // COMPLEX


int RK5::Corrector(double, double& errmax, int start, int stop)
{
	int j;
#pragma ivdep		
	for(j=start; j < stop; j++) {
#if COMPLEX	
		y[j].re=y0[j].re+c0*source0[j].re+c2*source2[j].re+c3*source3[j].re+
			c5*source[j].re;
		y[j].im=y0[j].im+c0*source0[j].im+c2*source2[j].im+c3*source3[j].im+
			c5*source[j].im;
#else	
		y[j]=y0[j]+c0*source0[j]+c2*source2[j]+c3*source3[j]+c5*source[j];
#endif // COMPLEX
	}
	
	if(dynamic) {
#pragma ivdep		
		for(j=start; j < stop; j++) {
#if COMPLEX	
			y3[j].re=y0[j].re+d0*source0[j].re+d2*source2[j].re+
				d3*source3[j].re+d4*source4[j].re+d5*source[j].re;
			y3[j].im=y0[j].im+d0*source0[j].im+d2*source2[j].im+
				d3*source3[j].im+d4*source4[j].im+d5*source[j].im;
#else
			y3[j]=y0[j]+d0*source0[j]+d2*source2[j]+
				d3*source3[j]+d4*source4[j]+d5*source[j];
#endif
		}
		for(j=start; j < stop; j++)
			calc_error(y0[j],y[j],y3[j],y[j],errmax);
		ExtrapolateTimestep(errmax);
	}
	return 1;
}
