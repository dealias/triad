#include "Array.h"

using namespace Array;

enum Palette {RAINBOW,BRAINBOW,WRAINBOW,BWRAINBOW,WHEEL,RGREYB,
			  RED,GREEN,BLUE,YELLOW,CYAN,MAGENTA,REDBLUE,REDGREEN,GREENBLUE,
			  GENERAL,NPalette};

enum AllColors {BLACK,WHITE,FirstColor};
extern int NColors;

const int ExtraColors=2; // Black, White

typedef unsigned char u_char;

extern Array1<u_char> Red, Blue, Green;
extern Array1<u_char> Y, U, V;

// RGB to YUV conversion factors

static const double cr=0.299;
static const double cg=0.587;
static const double cb=0.114;
static const double cu=0.5/(1.0-cb);
static const double cv=0.5/(1.0-cr);
  
inline unsigned char RGBByte(double r) // Normalize to interval [0,255]
{
  int a=(int)(256.0*r);
  if(a == 256) a=255;
  if(a < 0 || a > 255) msg(ERROR,"Invalid color: %d",a);
  return (unsigned char) a;
}

inline unsigned char YByte(double y) // Normalize to interval [16,235]
{ 
  int a=(int)(220.0*y)+16;
  if(a == 236) a=235;
  if(a < 16 || a > 235) msg(ERROR,"Invalid Y value: %d",a);
  return (unsigned char) a;
}

inline unsigned char UVByte(double u) // Normalize to interval [16,240]
{
  int a=(int)(225.0*(u+0.5))+16;
  if(a == 241) a=240;
  if(a < 16 || a > 240) msg(ERROR,"Invalid U or V value: %d",a);
  return (unsigned char) a;
}

