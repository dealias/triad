/* Test of PolarAverage routine */

#include "options.h"
#include "Polar.h"
#include <stdio.h>

static double unity(double k, double p, double q,
					double a, double b, double g, double lo, double hi)
	// Unity mode-coupling function $f$ (for testing).
{
	return hi-lo; // Unity weight factor
}

static double Leithv(double k, double p, double q,
					 double a, double b, double g, double lo, double hi)
	// Leith's mode-coupling factor $\Abs{sin(\b-\g)}/2k$.
{
	return (hi-lo)*fabs(sin(b-g))/(2.0*k); // Triangle volume fraction.
}

static double Leiths(double k, double p, double q,
					 double a, double b, double g, double lo, double hi)
	// Leith's mode-coupling function $\Abs{sin(\b-\g)}^2\max(k,p,q)/(2k^2)$.
{
	double mc,sina;
	sina=fabs(sin(b-g));
	mc=sina/(2.0*k);
	// Evaluate Leith's subsidiary correction factor
	// (in his notation, $sin\ab$).
	mc*=sina*max(max(k,p),q)/k; 
	return (hi-lo)*mc;
}

static double Leithexact(double k, double p, double q,
						 double a, double b, double g, double lo, double hi)
	// Leith's mode-coupling function $\Abs{sin(\b-\g)}^2\max(k,p,q)/(2k^2)$.
{
	double sina;
	sina=fabs(sin(b-g));
	return (hi-lo)*p*p*q*q*sina*sina;
}

inline void make_bin(Bin<Polar> *k, char *a, char *b, char *c, char *d) {
	build(k,atof(a),atof(b),atof(c),atof(d));
	if(k->max.th == -1.0) k->max.th=twopi;
}

inline void output_bin(Bin<Polar> *k, char *s, char *t) {
	printf(" %sl = %.5e, %sg = %.5e, %sl = %.5f, %sg = %.5f\n",
		   s,k->min.r,s,k->max.r,t,k->min.th,t,k->max.th);
}


static Real acc;
static Real fpar;
static Real v(int j, int l);
extern int PolarAverageCount;

int main(int argc, char *argv[])
{
	Bin<Polar> k,p,q,k1,p1,q1;
	Real area;
	Real ratio,totarea,sumarea; // Used for checksum test
	int nrad,nang,sumcount;
	int i,j,l,m,n,r,s; // Counters for partitioning intervals

	acc=1.0E-4; // Default relative accuracy of the answer
	nrad=1; // number of radial divisions
	nang=1; // number of angular divisions

	if(argc != 1) {
		if(argc < 13 && argc > 16) {
			msg(ERROR,"Invalid number (%d) of command-line arguments",
				argc);
		}
		// Process command line.
		
		make_bin(&k,argv[1],argv[2],argv[7],argv[8]);
		make_bin(&p,argv[3],argv[4],argv[9],argv[10]);
		make_bin(&q,argv[5],argv[6],argv[11],argv[12]);

		// Override defaults:
		if(argc >= 14) acc=atof(argv[13]); 
		if(argc >= 15) nrad=atoi(argv[14]);
		if(argc == 16) nang=atoi(argv[15]);


		printf(" acc = %.10e\n\n",acc);
		
		output_bin(&k,"k","a");
		output_bin(&p,"p","b");
		output_bin(&q,"q","g");
		
		// Calculate geometric factor for given bin boundary parameters
		area=BinAverage(&k,&p,&q,unity,acc);
		printf("area = %.10e,   count = %8d\n\n\n",area,PolarAverageCount);

		if(nrad == 1 && nang == 1) exit(1);
		totarea=area; // Setup for checksum calculation
		sumarea=0.0;
		sumcount=0;

		// partition the intervals
		Polar ndiv((double) nrad, (double) nang);
		
		for(i=0; i < nrad; ++i) for(j=0; j < nang; ++j) {
			Polar delta=k.Delta(), incr((double) i,(double) j);
			k1.min=k.min+incr*delta/ndiv;
			k1.max=k1.min+delta/ndiv;
			
			for(m=0; m < nrad; ++m) for(n=0; n < nang; ++n) {
				Polar delta=p.Delta(), incr((double) m,(double) n);
				p1.min=p.min+incr*delta/ndiv;
				p1.max=p1.min+delta/ndiv;
				
				for(r=0; r < nrad; ++r) for(s=0; s < nang; ++s) {
					Polar delta=q.Delta(), incr((double) r,(double) s);
					q1.min=q.min+incr*delta/ndiv;
					q1.max=q1.min+delta/ndiv;
					
					printf(" i =%3d  j =%3d  m =%3d  n =%3d  r =%3d  s =%3d\n",
						   i,j,m,n,r,s);
					
					output_bin(&k1,"k","a");
					output_bin(&p1,"p","b");
					output_bin(&q1,"q","g");
					
					area=BinAverage(&k1,&p1,&q1,unity,acc);
					printf("area = %.10e, count = %8d\n\n", area,
						   PolarAverageCount);

					sumarea+=area;
					sumcount+=PolarAverageCount;
				}
			}
		}

		printf("sumarea = %.10e,   sumcount = %8d\n",sumarea,sumcount);
		ratio=(totarea == 0.0) ? 1.0: sumarea/totarea;
		printf("ratio = %.5e\n",ratio);
		if(fabs(1.0-ratio) > acc) printf(" WARNING: CHECKSUM ERROR!");
		exit(1);
	}
	
	// Generate Leith's Table 1.
	fpar=4.0; // $F=4$
	printf(" q ");
	// j, l represent Leith's $p, q$ indices.
	for(j=0; j <= 4; j++) printf("  p = %2d ",j);
	
	printf("\n");

	for(l=0; l <= 15; l++) {
		printf("%2d",l);
		for(j=0; j <= 4; j++) printf("  %1.5f",v(j,l));
		printf("\n");
	}
}


static Real v(int j, int l) // Leith's function $\bar v$
{
	Bin<Polar> k,p,q;
	Real area,areas,areav,k0,p0,q0,dk,dp,dq;
	int m=l-j;

	build(&k,pow(2.0,(l-0.5)/fpar),pow(2.0,(l+0.5)/fpar),0.0,twopi);
	build(&p,pow(2.0,(m-0.5)/fpar),pow(2.0,(m+0.5)/fpar),0.0,twopi);
	build(&q,pow(2.0,( -0.5)/fpar),pow(2.0,( +0.5)/fpar),0.0,twopi);

	dk=k.Delta().r;
	dp=k.Delta().r;
	dq=k.Delta().r;

	k0=sqrt(k.min.r*k.max.r);
	p0=sqrt(p.min.r*p.max.r);
	q0=sqrt(q.min.r*q.max.r);

	areas=BinAverage(&k,&p,&q,Leithexact,acc);
	areav=BinAverage(&k,&p,&q,Leithv,acc);
	
	if(areav == 0.0)
		if(areas == 0.0)area=0.0;
		else msg(ERROR,"Division by zero attempted");
	else area=areas/areav*max(max(k0,p0),q0)/(2.0*p0*p0*q0*q0*k0*k0);
	return area;
	// return area/(twopi*dk*dp*dq);
}
