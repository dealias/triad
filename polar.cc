// POLARAVG

// INTRODUCTION. This documentation describes the polar coordinate
// calculation of the geometric coefficients that arise in the Fourier space
//integration of a 2-D statistical closure. The input arguments to the
function |polar| are the six coordinates specifying the magnitude and angular
limits of the three bins which are integrated over.

@ FWEB organization.
@a
@<Common@>
@<Utilities@>
@<Calc@>

#define ALLOCATE 0	/* Set to 1 for testing. */
#if(ALLOCATE)
@<Test@>
#endif

@ Link this code with
\begintt
link/deb polar,simp
\endtt

@i formats.web
@i mixlang.web 

@<Common@>=
#define ALLOCATE 0
#define NULL 0

#include "common.c"	/* Standard include */

#define DOUBLE double

@* UTILITIES. First we define some utility routines and global constants.


@ We will need to interchange the input arguments to |polar| to achieve the
desired ordering $\Dg\ge\Db\ge\Da$. Subject to this ordering, we
would like to further order $\qg-\ql\ge\kg-\kl$ if possible for reasons
of computational efficiency.
@<Utilities@>=

typedef struct 
{
	DOUBLE l,g;		/* magnitude limits of bin */
	DOUBLE tl,tg;		/* angular limits of bin */
} XK;

swap_ptr(p,q)
	XK **p,**q;
{
	XK *temp;
	temp=*p;
	*p=*q;
	*q=temp;
}

volatile double diff1,diff2; /* Don't optimize these away. */

#define SORT(xk,xq,level) \
if( diff1=(xq->tg - xq->tl), diff2=(xk->tg - xk->tl), \
   diff1 < diff2 || (diff1==diff2 && (xq->g - xq->l) < (xk->g - xk->l)) \
   ) {swap_ptr(&xk,&xq); cnt+=level;}

@* CALCULATIONS. This is the user entry point for the calculation. We sort
the arguments and then pass them on to |polar0|, which actually does the work.
@<Calc@>=

static double acc;		/* desired relative accuracy */
static double (*f)();		/* pointer to function $\fkpq$ */
static double eps,epsplus1;

static int perm;		/* index to our permutation of the original
						   input parameters. */

static XK xk0,xp0,xq0,*xk,*xp,*xq;

DOUBLE polar(kl,kg,pl,pg,ql,qg,al,ag,bl,bg,gl,gg,f0,acc0)
	DOUBLE kl,kg,pl,pg,ql,qg,al,ag,bl,bg,gl,gg,acc0;
	double (*f0)();
{
	double polar0();
	int cnt;/* Index to keep track of which variables have been swapped. */
	int permcnt[]={0,3,1,2,5,4,2,1};

	f=f0; acc=acc0;
	twopi=2.0*pi;

	eps=4.0*DBL_EPSILON;
	epsplus1=1.0+eps;

	xk=&xk0; xp=&xp0; xq=&xq0;

	xk->l=kl; xk->g=kg; xk->tl=al; xk->tg=ag;
	xp->l=pl; xp->g=pg; xp->tl=bl; xp->tg=bg;
	xq->l=ql; xq->g=qg; xq->tl=gl; xq->tg=gg;

	cnt=0;
	SORT(xk,xq,4);
	SORT(xp,xq,2);

	if(xk->tg - xk->tl > xp->tg - xp->tl) {swap_ptr(&xk,&xp); cnt+=1;}
	/* w.l.o.g. $\Db\ge\Da$ */

	perm=permcnt[cnt];/* index to permutation resulting from all of these sorts. */

	return polar0(xk->l,xk->g,xp->l,xp->g,xq->l,xq->g,xk->tl,xk->tg,
				  xp->tl,xp->tg,xq->tl,xq->tg);
}

@ This routine checks the validity of the input arguments.
@<Calc@>=

checklim(al,ag,offset)
	double *al,*ag,offset;
{
	*al-=offset; *ag-=offset;
	if(*al > *ag)errorfrom("checklim","Invalid limits (%.15g > %.15g)",*al,*ag);
	if(*ag-*al > twopi)*ag=*al+twopi;
}

@ Calculate the geometric coefficients for the special case where 
$\Dg\ge\Db\ge\Da$.
@<Calc@>=

static double kl,kg,pl,pg,ql,qg,al,ag,bl,bg,gl,gg;	/* input parameters:
													   $\kl,\kg,\pl,\pg,\ql,\qg,\al,\ag,\bl,\bg,\gl,\gg$ */	
static double dalpha,dbeta;  /* $\Da,\Db$ */
static double alpb,roff,ql2,qg2,dxmax3,a0;
static int count,iflag;

double polar0(kl0,kg0,pl0,pg0,ql0,qg0,al0,ag0,bl0,bg0,gl0,gg0)
	DOUBLE kl0,kg0,pl0,pg0,ql0,qg0,al0,ag0,bl0,bg0,gl0,gg0;
{
	static double a,b,ans,error,dxmax;
	double pint(),offset,kc,pc,qc,ac,bc,gc;
	double low,high;
	int iflag;
	FORTRAN(void, simp)();

	if(kg0 == kl0 || pg0 == pl0)return(0.0);
	if(kg0 <  kl0 || pg0 <  pl0)
	errorfrom("polar0","Invalid wavenumber magnitude limits");

	kl=kl0; kg=kg0; al=al0; ag=ag0;
	pl=pl0; pg=pg0; bl=bl0; bg=bg0;
	ql=ql0; qg=qg0; gl=gl0; gg=gg0;

	offset=gl;
	checklim(&al,&ag,offset); checklim(&bl,&bg,offset); checklim(&gl,&gg,offset);
	dalpha=ag-al;
	dbeta=bg-bl;
	alpb=pbranch(al);
	gl-=alpb;	/* form $\gl-\al$, where $\al$ is principal angle */
	gg-=alpb;	/* form $\gg-\al$ for use in \Eq(angular)*/
	roff=pbranch(bl-al);
	a0=al;		/* Used in |rint0| */
	
	dxmax3=pi;		/* |dxmax| for |rint0| */

	ql2=ql*ql;		
	qg2=qg*qg;
	count=0;

	a=pl;
	b=pg;
	dxmax=0.4*(qg-ql);

	simp(pint,&a,&b,&acc,&ans,&error,&dxmax,&iflag);
	/* This is the adaptive Simpson integration routine (Fortran) */
	if(iflag != 1) errorfrom("polar0","simp returned code %d",iflag);
	/* see documentation in |simp| for a description of possible errors */
	return ans;
}

@ Evaluate the $p$-integrand.
@<Calc@>=

static double p,p2,ql2mp2,qg2mp2;
double pint(p0)
	double *p0;
{
	static double a,b,ans,error,dxmax;
	double kint();
	FORTRAN(void, simp2)();

	p=*p0;
	p2=p*p;
	ql2mp2=ql2-p2;		/* Used later in solving \Eq(rrestrict) */
	qg2mp2=qg2-p2;

	a=kl;
	b=kg;
	dxmax=0.4*(qg-ql);

	simp2(kint,&a,&b,&acc,&ans,&error,&dxmax,&iflag);
	if(iflag != 1) errorfrom("pint","simp2 returned code %d",iflag);
	return p*ans;
}
@ Evaluate the $k$-integrand. Note that the integrand is evaluated only for
values of $\rb$ satisfying the magnitude restriction.
@<Calc@>=

static double k,k2,twokp,term;
static int coeffa,coeffb;

double kint(k0)
	double *k0;
{
	double a,b,ans,error;
	double rint0(),sum,cosr1,cosr2,r1,r2,r3,r4,r1true;
	int i,j,loopcount;
	FORTRAN(void, simp3)();

	k=*k0;
	k2=k*k;
	twokp=2.0*k*p;
	if(twokp != 0.0)
	{
		cosr1=(qg2mp2-k2)/twokp;	/* solve \Eq(rrestrict) in the
									   interval $[0,\pi]$*/
		cosr2=(ql2mp2-k2)/twokp;
		if(fabs(cosr1) <= 1.0)r1=acos(cosr1);
		else r1= (cosr1 > 0.0) ? 0.0 : pi;   /* If |cosr1| is out of range
												then $r_1$ is set to $0$ or $\pi$ */
		if(fabs(cosr2) <= 1.0)r2=acos(cosr2);
		else r2= (cosr2 > 0.0) ? 0.0 : pi;
	}				/* Note $0\le r_1 \le r_2 \le \pi$ */
	if(p==0.0)				
	{r1=0.0;			
	 if(ql <= k && k < qg)r2=pi;	/* Special case of the
									   magnitude restriction */
	 else r2=r1;		/* magnitude restriction is violated */
 }
	if(k==0.0)
	{r1=0.0;
	 if(ql <= p && p < qg)r2=pi;
	 else r2=r1;
 }

	r3=twopi-r2; r4=twopi-r1;	/* $r_{3,4}$ are from interval $[\pi,2\pi]$ */
	r1-=roff; r2-=roff;		/* Subtract $r_{off}$ to obtain $\rb$ values */
	r3-=roff; r4-=roff;
	sum=0.0;

	for(i=1;i<=2;++i)	/* Loop twice: First for $[\rb_1,\rb_2]$
						   and then for $[\rb_3,\rb_4]$ */
	{
		if(r2 < 0.0){r1+=twopi; r2+=twopi;}
		r1true=r1;
		loopcount=1;
		if(r1 < 0.0){r1=0.0; loopcount=2;}	/* Split the interval into two
											   subintervals and do the following procedure for
											   each subinterval */

		for(j=1;j<=loopcount;++j)	/* Now evaluate \Eq(rbarint) */
		{
			b=MIN(dbeta,r2);
			a=MAX(dbeta-dalpha,r1);
			coeffa=0; term=dbeta; coeffb=-1; /* $\abgr=\Db-\rb$ */
			/* $\ablr=0$ */
#define EPS ((a > 0) ? eps : -eps)	

			if(b > a*(1.0+EPS))
			{
				simp3(rint0,&a,&b,&acc,&ans,&error,&dxmax3,&iflag);
				if(iflag != 1) errorfrom("kint",
										 "First call to simp3 returned code %d",iflag);
				sum+=ans;
			}

			b=MIN(dbeta-dalpha,r2);
			a=MAX(0.0,r1);
			term=dalpha; coeffb=0;		 /* $\abgr=\Da$ */
			/* $\ablr=0$ */
			if(b > a*(1.0+EPS))
			{
				simp3(rint0,&a,&b,&acc,&ans,&error,&dxmax3,&iflag);
				if(iflag != 1) errorfrom("kint",
										 "Second call to simp3 returned code %d",iflag);
				sum+=ans;
			}

			b=MIN(0.0,r2-twopi);
			a=MAX(-dalpha,r1-twopi);
			coeffa=-1;		 /* $\abgr=\Da$ */
			/* $\ablr=-r$ */
			if(b > a*(1.0+EPS))
			{
				simp3(rint0,&a,&b,&acc,&ans,&error,&dxmax3,&iflag);
				if(iflag != 1) errorfrom("kint",
										 "Third call to simp3 returned code %d",iflag);
				sum+=ans;
			}
			r1=twopi+r1true;	/* Setup for second sub-interval (in case */
			r2=twopi;		/* of a split $[\rb_1,\rb_2]$ interval) */
		}
		r1=r3;				/* Treat $[\rb_3,\rb_4]$ interval similarly */
		r2=r4;
	}
	return k*sum;                                  
}
@ Evaluate the $\rb$-integrand.
@<Calc@>=

#ifndef sincos
#define sincos(x,sinx,cosx) (*sinx=sin(x), *cosx=cos(x))
#endif

double rint0(rbar0)
	double *rbar0;
{
	double q,theta,x,y,sum;
	double alphabl,alphabg,r,rbar,cosr,sinr,b0,g0,lo,hi;
	int n;
	rbar=*rbar0;
	count+=1;	  /* counter for total number of function evaluations */
	r=rbar+roff;	  /* $r=\rb+r_{off}$ */
	alphabg=term+coeffb*rbar; /* $\abgr=|term|+|coeffb|*\rb$ */
	alphabl=coeffa*rbar;	  /* $\ablr=|term|+|coeffb|*\rb$ */

	sincos(r,&sinr,&cosr);

	q=sqrt(k2+p2+twokp*cosr); /* $q=\abs{\vk+\vp}$ */
	if(q < ql*0.99999999 || q >= qg*1.00000001)
	errorfrom("rint0","Magnitude of q is out of bounds");

	y=-p*sinr; x=-k-p*cosr;
	if(x==0.0) theta=0.0; 		/* Note |atan2(y,0.0)| is undefined. */
	else theta=atan2(y,x);		/* Returns angle in $[-\pi,\pi]$ */
	if(theta < 0.0)theta+=twopi;	/* We want $\th\in[0,2\pi]$ */

	b0=rbar+bl; g0=theta+al;
	/* $\a=\ab+\a_0, \b=\ab+b_0, \g=\ab+g_0$ */ 

	sum=0.0;
	for(n=0;n<=2;n++)	/* n indexes offsets of $0$, $2\pi$, $4\pi$ */
	{
		/* In the following, |gg| is actually $\gg-\al$ */
		hi = MIN(alphabg, gg+n*twopi-theta);
		lo = MAX(alphabl, gl+n*twopi-theta);

		if(hi > lo*epsplus1) switch(perm)	/* \Eq(angular) */
		{
		case 0: sum +=  (*f)(&k,&p,&q,&a0,&b0,&g0,&lo,&hi); break;
			case 1: sum +=  (*f)(&k,&q,&p,&a0,&g0,&b0,&lo,&hi); break;
			case 2: sum +=  (*f)(&p,&q,&k,&b0,&g0,&a0,&lo,&hi); break;
			case 3: sum +=  (*f)(&p,&k,&q,&b0,&a0,&g0,&lo,&hi); break;
			case 4: sum +=  (*f)(&q,&k,&p,&g0,&a0,&b0,&lo,&hi); break;
			case 5: sum +=  (*f)(&q,&p,&k,&g0,&b0,&a0,&lo,&hi); break;
			}
	}
	return sum;
}

@ Mode-coupling functions used for testing.
@<Test@>=

double unity(k0,p0,q0,a0,b0,g0,lo,hi)
	/* Unity mode-coupling function $f$ (for testing). */
	double *k0,*p0,*q0,*a0,*b0,*g0,*lo,*hi;
{
	return (*hi)-(*lo); /* Unity weight factor */
}

double Leithv(k0,p0,q0,a0,b0,g0,lo,hi)
	/* Leith's mode-coupling factor $\Abs{sin(\b-\g)}/2k$. */
	double *k0,*p0,*q0,*a0,*b0,*g0,*lo,*hi;
{
	double k,p,q,a,b,g;
	k=*k0; p=*p0; q=*q0;
	a=*a0; b=*b0; g=*g0;
	return ((*hi)-(*lo))*fabs(sin(b-g))/(2.0*k); /* Triangle volume fraction. */
}

double Leiths(k0,p0,q0,a0,b0,g0,lo,hi)
	/* Leith's mode-coupling function $\Abs{sin(\b-\g)}^2\max(k,p,q)/(2k^2)$. */
	double *k0,*p0,*q0,*a0,*b0,*g0,*lo,*hi;
{
	double k,p,q,a,b,g;
	double mc,sina;
	NO_REF(k0;p0;q0;a0;b0;g0);
	k=*k0; p=*p0; q=*q0;
	a=*a0; b=*b0; g=*g0;
	sina=fabs(sin(b-g));
	mc=sina/(2.0*k);
	mc*=sina*MAX(MAX(k,p),q)/k; /* Used to evaluate Leith's subsidiary
								   correction factor (in his notation, $sin\ab$). */
	return ((*hi)-(*lo))*mc;
}

double Leithexact(k0,p0,q0,a0,b0,g0,lo,hi)
	/* Leith's mode-coupling function $\Abs{sin(\b-\g)}^2\max(k,p,q)/(2k^2)$. */
	double *k0,*p0,*q0,*a0,*b0,*g0,*lo,*hi;
{
	double k,p,q,a,b,g;
	double sina;
	NO_REF(k0;p0;q0;a0;b0;g0);
	k=*k0; p=*p0; q=*q0;
	a=*a0; b=*b0; g=*g0;
	sina=fabs(sin(b-g));
	return ((*hi)-(*lo))*p*p*q*q*sina*sina;
}

@* TESTING. Here is a main program which accepts a command line of the form
\begintt
q kl kg pl pg ql qg al ag bl bg gl gg acc
\endtt
where you must replace~\.kl, \.kg, etc. with numerical values and where
the code has been defined to be a foreign procedure \via\
\begintt
q :== $sys$login:polar.exe
\endtt
If~\.{acc}, the relative accuracy to which the integral is to be computed,
is not specified then a default value will be used. If the value $-1$ is
given for one of the angular variables, the value $\2\pi$ will be used.
Alternatively, if no arguments are specified
the table from Leith's paper ('72) will be generated with the default accuracy.
@<Test@>=

static double fpar;
double unity();
FORTRAN(void, precis)(); 	/* Fortran machine precision. */

main(num_args,arg)
	UNSIGNED num_args;
	char *arg[];
{
	DOUBLE  kl,kg,pl,pg,ql,qg;	/* Magnitude input parameters. */
	DOUBLE  al,ag,bl,bg,gl,gg;	/* Angular input parameters. */
	DOUBLE  kl1,kg1,pl1,pg1,ql1,qg1; /* Used by checksum test for splitting each */
	DOUBLE  al1,ag1,bl1,bg1,gl1,gg1; /*   main interval into |nrad| radial and */
	/*   |nang| angular sub-intervals. */
	double atof(),area;
	double ratio,totarea,sumarea;	/* Used for checksum test */
	int sumcount;
	double v();
	int i,j,k,l,m,n;		/* Counters for partitioning intervals */
	int nrad,nang;

	pi=acos(-1.0);
	precis();			/* Set up internal Fortran precision. */

	twopi=2.0*pi;
	acc=1.0E-2;			/* Default relative accuracy of the answer */
	nrad=1;
	nang=1;
	if(num_args != 1)		
	{
		if((num_args != 13) & (num_args != 14))
		errorfrom("main",
				  "Invalid number (%d) of command-line arguments",num_args);

		/* Process command line. */
		kl = atof(arg[1]);
		kg = atof(arg[2]);
		pl = atof(arg[3]);
		pg = atof(arg[4]);
		ql = atof(arg[5]);
		qg = atof(arg[6]);
		al = atof(arg[7]);
		ag = atof(arg[8]);
		bl = atof(arg[9]);
		bg = atof(arg[10]);
		gl = atof(arg[11]);
		gg = atof(arg[12]);
		if(num_args == 14)acc= atof(arg[13]);	/* Override default accuracy */

		if(ag==-1.0)ag=twopi;
		if(bg==-1.0)bg=twopi;
		if(gg==-1.0)gg=twopi;
		
		/* Calculate geometric factor for given bin boundary parameters */

		printf(" acc = %.10e\n\n",acc);	
		printf(" kl = %.5e, kg = %.5e, al = %.5f, ag = %.5f\n",kl,kg,al,ag);
		printf(" pl = %.5e, pg = %.5e, bl = %.5f, bg = %.5f\n",pl,pg,bl,bg);
		printf(" ql = %.5e, qg = %.5e, gl = %.5f, gg = %.5f\n",ql,qg,gl,gg);
		area=polar(kl,kg,pl,pg,ql,qg,al,ag,bl,bg,gl,gg,unity,acc);
		printf("area = %.10e,   count = %8d\n\n\n",area,count);

		if(nrad==1 && nang==1)exit();
		totarea=area;			/* Setup for checksum calculation */
		sumarea=0.0;
		sumcount=0.0;

#define LOOPP(i,nrad) \
		for(i=0;i<=nrad-1;++i)
#define PART(kl1,kg1,kl,kg,i,nrad) \
		kl1=kl+i*(kg-kl)/nrad; kg1=kl1+(kg-kl)/nrad; 

		/* PARTition the intervals */
		LOOPP(i,nrad) { PART(kl1,kg1,kl,kg,i,nrad);
						LOOPP(j,nrad) { PART(pl1,pg1,pl,pg,j,nrad);
										LOOPP(k,nrad) { PART(ql1,qg1,ql,qg,k,nrad);
														LOOPP(l,nang) { PART(al1,ag1,al,ag,l,nang);
																		LOOPP(m,nang) { PART(bl1,bg1,bl,bg,m,nang);
																						LOOPP(n,nang) { PART(gl1,gg1,gl,gg,n,nang);

																										printf(" i =%3d  j =%3d  k =%3d  l =%3d  m =%3d  n =%3d\n",i,j,k,l,m,n);
																										printf(" kl = %.5e, kg = %.5e, al = %.5f, ag = %.5f\n",kl1,kg1,al1,ag1);
																										printf(" pl = %.5e, pg = %.5e, bl = %.5f, bg = %.5f\n",pl1,pg1,bl1,bg1);
																										printf(" ql = %.5e, qg = %.5e, gl = %.5f, gg = %.5f\n",ql1,qg1,gl1,gg1);

																										area=polar(kl1,kg1,pl1,pg1,ql1,qg1,al1,ag1,bl1,bg1,gl1,gg1,unity,acc);
																										sumarea+=area;
																										sumcount+=count;
																										printf("area = %.10e,   count = %8d\n\n",area,count);
																									}}}}}}   

		printf("sumarea = %.10e,   sumcount = %8d\n",sumarea,sumcount);
		ratio=(totarea==0.0) ? 1.0 : sumarea/totarea; 
		printf("ratio = %.5e\n",ratio);
		if(fabs(1.0-ratio) > acc)printf(" WARNING: CHECKSUM ERROR!");
		exit();
	}
	/* Generate Leith's Table 1. */
	fpar=4.0;			/* $F=4$ */
	{
		printf(" q ");
		for(j=0;j<=4;j++)	/* Note that |j, l| represent Leith's $p, q$ indices. */
		{
			printf("  p = %2d ",j);
		}
		printf("\n");

		for(l=0;l<=15;l++)
		{
			printf("%2d",l);
			for(j=0;j<=4;j++)
			{
				printf("  %1.5f",v(j,l));
			}
			printf("\n");
		}
	}
}
@ Calculate Leith's function $\bar v(l-m,l)$, (taking $n=0$).
  Note that |j, l| represent Leith's $p, q$ indices.
  @<Test@>=

  double v(j,l)		/* Leith's function $\bar v$ */
	int j,l;
{
	DOUBLE polar();
	double area;
	DOUBLE areas,areav,k0,p0,q0;
	DOUBLE kl,kg,dk,pl,pg,dp,ql,qg,dq,ag;
	int m;
	m=l-j;

	kl=pow(2.0,(l-0.5)/fpar);
	kg=pow(2.0,(l+0.5)/fpar);
	pl=pow(2.0,(m-0.5)/fpar);
	pg=pow(2.0,(m+0.5)/fpar);
	ql=pow(2.0,-0.5/fpar);
	qg=pow(2.0,0.5/fpar);

	dk=kg-kl;
	dp=pg-pl;
	dq=qg-ql;

	ag=twopi;
	k0=sqrt(kl*kg);
	p0=sqrt(pl*pg);
	q0=sqrt(ql*qg);

	areas=polar(kl,kg,pl,pg,ql,qg,0.0,ag,0.0,ag,0.0,ag,Leithexact,acc);
	areav=polar(kl,kg,pl,pg,ql,qg,0.0,ag,0.0,ag,0.0,ag,Leithv,acc);
	if(areav==0.0)
	if(areas==0.0) area=0.0;
	else errorfrom("v","Division by zero attempted");
	else area=areas/areav*MAX(MAX(k0,p0),q0)/(2.0*p0*p0*q0*q0*k0*k0);
	return area;
	/* return area/(twopi*dk*dp*dq); */
}
