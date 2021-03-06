#include "utils.h"
#include "rgb.h"

using namespace Array;

int NColors=65536;

Array1<u_char> Red, Blue, Green;
Array1<u_char> Y, U, V;

extern int verbose;
extern int two;
extern int gradient;
extern int damp;
extern double r1,g1,b1;
extern double r2,g2,b2;
extern int reverse;

static int k,incr;

void AddColor(double r, double g, double b) 
{
  if(damp) {
    Real factor=sqrt(r*r+g*g+b*b);
    if(factor != 0.0) factor=1.0/factor;
    r *= r*factor;
    g *= g*factor;
    b *= b*factor;
  }
	
  if(gradient) {
    Real factor=((double) k)/NColors;
    r *= factor;
    g *= factor;
    b *= factor;
  }

  Red[k]=RGBByte(r);
  Green[k]=RGBByte(g);
  Blue[k]=RGBByte(b);
	
  double y=cr*r+cg*g+cb*b;

  Y[k]=YByte(y);
  U[k]=UVByte(cu*(b-y));
  V[k]=UVByte(cv*(r-y));
  
  k += incr;
}

void MakePalette(int palette)
{
  int i;
  int num,nintervals=1,offset=0;
  int divisor=1;
		
  switch(palette) {
  case RED:
  case GREEN:
  case BLUE:
  case YELLOW:
  case CYAN:
  case MAGENTA:
  case REDBLUE:
  case REDGREEN:
  case GREENBLUE:
  case GENERAL:
    offset=0;
    nintervals=1;
    break;
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
  case RGREYB:
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
		
  int allcolors=NColors+FirstColor;

  Red.Allocate(allcolors);
  Blue.Allocate(allcolors);
  Green.Allocate(allcolors);
  
  Y.Allocate(allcolors);
  U.Allocate(allcolors);
  V.Allocate(allcolors);

  // Define extra colors;
  Red[0]=Blue[0]=Green[0]=RGBByte(0.0); // BLACK
  Red[1]=Blue[1]=Green[1]=RGBByte(1.0); // WHITE
	
  Y[0]=YByte(0.0); U[0]=V[0]=UVByte(0.0); // BLACK
  Y[1]=YByte(1.0); U[0]=V[0]=UVByte(0.0);  // WHITE
  
  double ninv=1.0/n;

  incr=-1;
  k=allcolors-1;
  if(two) gradient=1;

  switch(palette) {
  case RGREYB:
    ninv *= 0.5;
    for(i=0; i < n; i++) AddColor((2*n-i)*ninv,i*ninv,i*ninv);
    for(i=0; i <= n; i++) AddColor((n-i)*ninv,(n-i)*ninv,(n+i)*ninv);
    return;
		
  case RED: 
    for(i=0; i <= n; i++) AddColor((n-i)*ninv,0.0,0.0);
    return;

  case GREEN: 
    for(i=0; i <= n; i++) AddColor(0.0,(n-i)*ninv,0.0);
    return;
		
  case BLUE: 
    for(i=0; i <= n; i++) AddColor(0.0,0.0,(n-i)*ninv);
    return;
		
  case YELLOW: 
    for(i=0; i <= n; i++) AddColor((n-i)*ninv,(n-i)*ninv,0.0);
    return;

  case CYAN:
    for(i=0; i <= n; i++) AddColor(0.0,(n-i)*ninv,(n-i)*ninv);
    return;
		
  case MAGENTA: 
    for(i=0; i <= n; i++) AddColor((n-i)*ninv,0.0,(n-i)*ninv);
    return;
		
  case REDBLUE: 
    for(i=0; i <= n; i++) AddColor((n-i)*ninv,0.0,i*ninv);
    return;
		
  case REDGREEN: 
    for(i=0; i <= n; i++) AddColor((n-i)*ninv,i*ninv,0.0);
    return;
		
  case GREENBLUE: 
    for(i=0; i <= n; i++) AddColor(0.0,(n-i)*ninv,i*ninv);
    return;
		
  case GENERAL:
    reverse=!reverse;
    for(i=0; i <= n; i++) AddColor(r1+r2*i*ninv,g1+g2*i*ninv,b1+b2*i*ninv);
    return;
  }

	
  if(palette == WHEEL) {
    for(i=0; i < n; i++) AddColor(1.0,0.0,(n-i)*ninv);
  }
	
  if(palette == WRAINBOW || palette == BWRAINBOW) {
    for(i=0; i < n; i++) AddColor(1.0,(n-i)*ninv,(n-i)*ninv);
  }
	
#if 0	
  if(palette == BWRAINBOW) {
    int n3=n/3;
    int n23=2*n3;
    for(i=0; i < n3; i++) AddColor(1.0,(n23-i)*1.5*ninv,(n-i)*ninv);
    for(i=0; i < n3; i++) AddColor(1.0,(n3-i)*1.5*ninv,(n23-i)*ninv);
    for(i=0; i < n3; i++) AddColor(1.0,0.0,(n3-i)*ninv);
  }
#endif	
	
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
  
  if(verbose) cout << "Constructed color palette with " << NColors 
		   << " colors." << endl;
}


