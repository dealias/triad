const int Ncolors=252;
const int PaletteMax=209;
typedef unsigned char u_char;

u_char red[Ncolors]={255,255,255,255,255,255,255,255,255,255,255,255,255,255,
					 255,255,255,255,255,255,255,255,255,255,255,255,255,255,
					 255,255,255,255,255,255,255,255,255,255,255,255,255,255,
					 255,248,242,236,230,224,218,212,206,200,194,188,182,176,
					 170,163,157,151,145,139,133,127,121,115,109,103,97,91,85,
					 78,72,66,60,54,48,42,36,30,24,18,12,6,0,0,0,0,0,0,0,0,0,0,
					 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
					 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
					 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,12,18,24,30,36,42,48,
					 54,60,66,72,78,85,91,97,103,109,115,121,127,133,139,145,
					 151,157,163,170,176,182,188,194,200,206,212,218,224,230,
					 236,242,248,255,255,255,255,255,255,255,255,255,255,255,
					 255,255,255,255,255,255,255,255,255,255,255,255,255,255,
					 255,255,255,255,255,255,255,255,255,255,255,255,255,255,
					 255,255,255};

u_char green[Ncolors]={0,6,12,18,24,30,36,42,48,54,60,66,72,78,85,91,97,103,
					   109,115,121,127,133,139,145,151,157,163,170,176,182,188,
					   194,200,206,212,218,224,230,236,242,248,255,255,255,255,
					   255,255,255,255,255,255,255,255,255,255,255,255,255,255,
					   255,255,255,255,255,255,255,255,255,255,255,255,255,255,
					   255,255,255,255,255,255,255,255,255,255,255,255,255,255,
					   255,255,255,255,255,255,255,255,255,255,255,255,255,255,
					   255,255,255,255,255,255,255,255,255,255,255,255,255,255,
					   255,255,255,255,255,255,255,255,255,255,255,248,242,236,
					   230,224,218,212,206,200,194,188,182,176,170,163,157,151,
					   145,139,133,127,121,115,109,103,97,91,85,78,72,66,60,54,
					   48,42,36,30,24,18,12,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
					   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
					   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
					   0,0,0,0,0,0,0,0,0,0,0,0};

u_char blue[Ncolors]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
					  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
					  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
					  0,6,12,18,24,30,36,42,48,54,60,66,72,78,85,91,97,103,
					  109,115,121,127,133,139,145,151,157,163,170,176,182,188,
					  194,200,206,212,218,224,230,236,242,248,255,255,255,255,
					  255,255,255,255,255,255,255,255,255,255,255,255,255,255,
					  255,255,255,255,255,255,255,255,255,255,255,255,255,255,
					  255,255,255,255,255,255,255,255,255,255,255,255,255,255,
					  255,255,255,255,255,255,255,255,255,255,255,255,255,255,
					  255,255,255,255,255,255,255,255,255,255,255,255,255,255,
					  255,255,255,255,255,255,255,255,255,255,255,248,242,236,
					  230,224,218,212,206,200,194,188,182,176,170,163,157,151,
					  145,139,133,127,121,115,109,103,97,91,85,78,72,66,60,54,
					  48,42,36,30,24,18,12,6};
