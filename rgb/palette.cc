#include "utils.h"
#include "rgb.h"

int NColors=65536;

Array1<u_char> Red;
Array1<u_char> Blue;
Array1<u_char> Green;

extern int verbose;
extern int two;
extern int gradient;

static int k,incr;

int ColorByte(double r) {
	return (int)(255.0*r+0.5);
}

void AddColor(double r, double g, double b) {
	Real factor=gradient ? ((double) k)/NColors : 1.0;
	Red[k]=ColorByte(r*factor);
	Green[k]=ColorByte(g*factor);
	Blue[k]=ColorByte(b*factor);
	k += incr;
}

void MakePalette(int palette)
{
	int i;
	int num,nintervals=1,offset=0;
	int divisor=1;
		
	switch(palette) {
	case RAINBOW:
		offset=1;
		nintervals=5;
		break;
	case BRAINBOW:
		offset=1;
		nintervals=5;
		divisor=3;
		break;
	case WRAINBOW:
	case WHEEL:
		offset=0;
		nintervals=6;
		break;
	case BWRAINBOW:
		offset=1;
		nintervals=6;
		divisor=3;
		break;
	case RGreyB:
		two=0; // Can't double this palette
		offset=1;
		nintervals=2;
		break;
	default:
		msg(ERROR,"Invalid palette: %d",palette);
	}
	
	if(two) nintervals += 6;
	
	num=NColors-offset;
	int n=(num/(nintervals*divisor))*divisor;
	NColors=n*nintervals+offset;
		
	if(verbose) cout << "Constructing color palette with " << NColors 
					 << " colors." << endl;
	
	int allcolors=NColors+FirstColor;

	Red.Allocate(allcolors);
	Blue.Allocate(allcolors);
	Green.Allocate(allcolors);

	// Define extra colors;
	Red[0]=Blue[0]=Green[0]=0; // BLACK
	Red[1]=Blue[1]=Green[1]=255; // WHITE
	
	double ninv=1.0/n;

	incr=-1;
	k=allcolors-1;
	if(two) gradient=1;

	if(palette == WHEEL) {
		for(i=0; i < n; i++) AddColor(1.0,0.0,(n-i)*ninv);
	}
	
	if(palette == RGreyB) {
		ninv *= 0.5;
		for(i=0; i < n; i++) AddColor((2*n-i)*ninv,i*ninv,i*ninv);
		for(i=0; i <= n; i++) AddColor((n-i)*ninv,(n-i)*ninv,(n+i)*ninv);
		return;
	}
	
	if(palette == WRAINBOW) {
		for(i=0; i < n; i++) AddColor(1.0,(n-i)*ninv,(n-i)*ninv);
	}
	
	if(palette == BWRAINBOW) {
		int n3=n/3;
		int n23=2*n3;
		for(i=0; i < n3; i++) AddColor(1.0,(n23-i)*1.5*ninv,(n-i)*ninv);
		for(i=0; i < n3; i++) AddColor(1.0,(n3-i)*1.5*ninv,(n23-i)*ninv);
		for(i=0; i < n3; i++) AddColor(1.0,0.0,(n3-i)*ninv);
	}
	
	for(i=0; i < n; i++) AddColor(1.0,i*ninv,0.0);
	for(i=0; i < n; i++) AddColor((n-i)*ninv,1.0,0.0);

	for(i=0; i < n; i++) AddColor(0.0,1.0,i*ninv);
	for(i=0; i < n; i++) AddColor(0.0,(n-i)*ninv,1.0);

	if((palette == BRAINBOW || palette == BWRAINBOW) && !gradient) {
		int n3=n/3;
		int n23=2*n3;
		for(i=0; i < n3; i++) AddColor(i*ninv,0.0,(n-i)*ninv);
		for(i=0; i < n3; i++) AddColor(n3*ninv,0.0,(n23-i)*ninv);
		for(i=0; i <= n3; i++) AddColor((n3-i)*ninv,0.0,(n3-i)*ninv);
	}

	if(palette == RAINBOW || palette == WRAINBOW || palette == WHEEL 
	   || gradient) {
		for(i=0; i <= n; i++) AddColor(i*ninv,0.0,1.0);
	}
	
	if(two) {
		for(i=0; i < n; i++) AddColor(1.0,0.0,(n-i)*ninv);
		for(i=0; i < n; i++) AddColor(1.0,i*ninv,0.0);
		for(i=0; i < n; i++) AddColor((n-i)*ninv,1.0,0.0);

		for(i=0; i < n; i++) AddColor(0.0,1.0,i*ninv);
		for(i=0; i < n; i++) AddColor(0.0,(n-i)*ninv,1.0);
		for(i=0; i <= n; i++) AddColor(i*ninv,0.0,1.0);
	}
}


