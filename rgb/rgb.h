#undef __ExternalArrayExit
#include "Array.h"

enum Palette {RAINBOW,BRAINBOW,WRAINBOW,BWRAINBOW,WHEEL,RGREYB,
			  RED,GREEN,BLUE,YELLOW,CYAN,MAGENTA,REDBLUE,REDGREEN,GREENBLUE,
			  GENERAL,NPalette};

enum AllColors {BLACK,WHITE,FirstColor};
extern int NColors;

const int ExtraColors=2; // Black, White

typedef unsigned char u_char;

extern Array1<u_char> Red;
extern Array1<u_char> Blue;
extern Array1<u_char> Green;

