#include "options.h"

const int nest=100;

typedef struct {
  int lorr;
  Real psum, f1t, f2t, f3t, dat, estr;
} TABLE;

int 				// Returns 1 if successful, 0 otherwise.
simpfast(Real (*f)(Real),	// Pointer to function to be integrated.
	 Real a, Real b,// Lower and upper limits of integration (a <= b).
	 Real acc, 	  	// Desired relative accuracy of sum. Routine tries to
	 // make fabs(error) <= acc*abs(sum).
	 Real& sum,		// Approximate value of the integral.
	 Real dxmax,	// Maximum limit on the width of a subinterval. For
	 // periodic functions, dxmax should be set to the period or smaller
	 // to prevent premature convergence of Simpson's rule.
		 
	 int& iflag 	// Error code:
	 // 0=successful,
	 // 1=nesting capacity exceeded,
	 // Note: This routine ignores underflow.
	 )
{
  Real diff, area, estl, estr, alpha, da, dx, wt, est, arg, fv[5];
  TABLE table[nest],*p,*pstop;
	
  iflag=0;
  p=table;
  pstop=table+nest-1;
  p->lorr=1;
  p->psum=0.0;
  alpha=a;
  da=b-a;
  fv[0]=(*f)(alpha);
  fv[2]=(*f)(alpha + 0.5*da);
  fv[4]=(*f)(alpha + da);
  wt=da/6.0;
  est=wt*(fv[0] + 4.0*fv[2] + fv[4]);
  area=est;

  //  Have estimate est of integral on (alpha, alpha+da).
  //  Bisect and compute estimates on left and right half intervals.
  //  Sum is better value for integral.

  for(;;) {
    dx=0.5*da;
    arg=alpha+0.5*dx;
    fv[1]=(*f)(arg);
    fv[3]=(*f)(arg+dx);
    wt=dx/6.0;
    estl=wt*(fv[0] + 4.0*fv[1] + fv[2]);
    estr=wt*(fv[2] + 4.0*fv[3] + fv[4]);
    sum=estl + estr;
    diff=est - sum;
    area -= diff;

    if(p >= pstop) iflag=1;
    if(iflag || fabs(diff) <= acc*fabs(area) && da <= dxmax) {
      //  Accept approximate integral sum. If it was a right interval,
      //  add results to finish at this level.  If it was a left
      //  interval, process right interval. Array lorr indicates left
      //  or right interval at each level.

      for(;;) {
	if (p->lorr == 0) { //  process right-half interval
	  alpha += da;
	  p->lorr=1;
	  p->psum=sum;
	  fv[0]=p->f1t;
	  fv[2]=p->f2t;
	  fv[4]=p->f3t;
	  da=p->dat;
	  est=p->estr;
	  break;
	}
	sum=p->psum + sum;
	if(--p <= table) return iflag;
      }

    } else {
      //  Raise level and store information to process right-half
      //  interval later. Initialise for basic step.
      ++p;
      da=dx;
      est=estl;
      p->lorr=0;
      p->f1t=fv[2];
      p->f2t=fv[3];
      p->f3t=fv[4];
      p->dat=dx;
      p->estr=estr;
      fv[4]=fv[2];
      fv[2]=fv[1];
    }
  }
}
