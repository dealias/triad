/* RGB:  A movie production utility
Copyright (C) 1998 John C. Bowman (bowman@math.ualberta.ca)

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. */

const char PROGRAM[]="RGB";
const char VERSION[]="1.06J";

#define NEW_IMAGEMAGICK 1

#if NEW_IMAGEMAGICK
char yuvformat[]="yuv";
char yuvinterlace[]="-interlace partition ";
#else
char yuvformat[]="yuv3";
char yuvinterlace[]="";
#endif 

#include "xstream.h"
#include <iostream.h>
#include <limits.h>
#include <errno.h>
#include <strstream.h>
#include <iomanip.h>
#include <unistd.h>
#include <string.h>
#include <sys/wait.h>

#include "getopt.h"
#include "DynVector.h"
#include "Array.h"
#include "rgb.h"

#ifndef PI
#define PI 3.14159265358979323846264338327950288
#endif

const double pi=PI;
const double twopi=2.0*PI;

// Number of digits used to index intermediate files.
const int NDIGITS=4;

static double *vminf, *vmaxf;
static char *rgbdir;
static strstream rgbdirbuf;
static int xsize,ysize;

static int verbose=0;
static int floating_scale=0;
static int floating_section=0;
static int preserve=0;
static int remote=0;
static int label=0;
static int pointsize=0;
static int grey=0;
static char *convertprog;
static strstream option;
static int Nx,Ny;
static int kmin,kmax;
static Real Pz;
static int a0,a;
static int R0,R,Rp;

int nx1=1,ny1=1,nz1=1;
int nx,ny,nz;
int sx=1,sy=1;
int byte=0;
int implicit=1;
int symmetric=0;
int invert=0;

enum Parameters {RFACTOR=256,THETA,PHI,YXASPECT,ZXASPECT,POINTSIZE,AVGX,AVGY,\
				 NXFINE,NYFINE,NZFINE};

Real Rfactor=2.0;
Real Theta=0.9;
Real Phi=0.9;
Real yxaspect=1.0;
Real zxaspect=1.0;

int Nxfine=4;
int Nyfine=5;
int Nzfine=10;

const int Undefined=-2;

void cleanup();
int system (char *command);

class Ivec {
public:
	float i;
	float j;
	float k;
	Ivec() {}
	Ivec(Real i0, Real j0, Real k0) : i(i0), j(j0), k(k0) {}
};

enum Transforms {IDENTITY,CIRCLE,TORUS};
typedef void Transform(Array2<Ivec>);
Transform Circle,Torus;
Transform *transform[]={NULL,Circle,Torus};
unsigned int Ntransform=sizeof(transform)/sizeof(Transform *);

template<class T>
void openfield(T& fin, char *fieldname, int& nx, int& ny, int& nz)
{
	fin.open(fieldname);
	if(!fin) {cleanup(); msg(ERROR,"Cannot open input file %s",fieldname);}
	if(implicit) {
		fin >> nx1 >> ny1 >> nz1;
		if(fin.eof()) {
			cleanup(); msg(ERROR,"End of file during processing");
		}
	}
	
	div_t d;
	d=div(nx1,sx);
	nx=(d.rem == 0) ? d.quot : d.quot+1;
	d=div(ny1,sy);
	ny=(d.rem == 0) ? d.quot : d.quot+1;
	nz=nz1;
				
}

int readframe(ixstream& xin, int nx, int ny, int nz, Array3<float> value,
			  double& gmin, double& gmax, double *vmink, double *vmaxk)
{
	gmin=DBL_MAX; gmax=-DBL_MAX;
	double vmin=DBL_MAX, vmax=-DBL_MAX;
	int j0;
	
	errno=0;
	for(int k=0; k < nz; k++) {
		Array2<float> valuek=value[k];
		int start,stop,incr;
		if(invert) {
			start=ny1-1;
			stop=-1;
			incr=-1;
			j0=ny;
		} else {
			start=0;
			stop=ny1;
			incr=1;
			j0=-1;
		}
		
		int rx=div(nx1,sx).rem;
		int ry=div(ny1,sy).rem;
		Real weightx=1.0/sx;
		Real weighty=1.0/sy;
		Real rweightx,rweighty;
		
		rweightx = (rx == 0.0) ? 1.0 : 1.0/rx;
		rweighty = (ry == 0.0) ? 1.0 : 1.0/ry;
		
		for(int j=start; j != stop; j += incr) {
			
			int init=0;
			if((j-start) % sy == 0) {j0 += incr; init=1;}
			Array1(float) valuekj=valuek[j0];
			if(init) for(int i0=0; i0 < nx; i0++) valuekj[i0]=0.0;
			
			Real sumv=0.0;
			int i0=0;
			for(int i=0; i < nx1; i++) {
				float v;
				if(byte) {
					xbyte x;
					xin >> x;
					v=x;
				}
				else xin >> v;
				
				if(xin.eof()) {
					if(implicit || i > 0)
						msg(WARNING,"End of file during processing");
					return EOF;
				}
				
				sumv += v*(i < nx1-rx ? weightx : rweightx)*
					((j-start)*incr < ny1-ry ? weighty : rweighty);
				
				if((i+1) % sx == 0 || (i+1) == nx1) {
					sumv += valuekj[i0];
					valuekj[i0++]=sumv;
					if((j-start+incr) % sy == 0) {
						if(sumv < vmin) vmin=sumv;
						if(sumv > vmax) vmax=sumv;
					}
					sumv=0.0;
				}
			}
		}
		
		if(symmetric && vmin < 0 && vmax > 0) {
			vmax=max(-vmin,vmax);
			vmin=-vmax;
		}
		vmink[k]=vmin;
		vmaxk[k]=vmax;
		if(vmin < gmin) gmin=vmin;
		if(vmax > gmax) gmax=vmax;
	}
	
	if(symmetric && gmin < 0 && gmax > 0) {
		gmax=max(-gmin,gmax);
		gmin=-gmax;
	}
	
	if(implicit) {
		int nx0,ny0,nz0;
		xin >> nx0 >> ny0 >> nz0;
		if(xin.eof()) return 1;
		if(nx0 != nx1 || ny0 != ny1 || nz0 != nz1) {
			cleanup();
			msg(ERROR,"Inconsistent image size");
		}
	}
	return 0;
}

void montage(int nfiles, char *const argf[], int n, char *const format,
			 char *const type);
void identify(int argc, char *const argf[], int n, char *const type,
			  int& xsize, int& ysize);
void mpeg(int argc, char *const argf[], int n, char *const type,
		  int xsize, int ysize);
void animate(int argc, char *const argf[], int n, char *const type,
			 const char *pattern, int xsize, int ysize);
	
#ifdef __i386__
extern "C" void putenv(const char *);
#endif

void usage(char *program)
{
	cerr << PROGRAM << " version " << VERSION
		 << " [(C) John C. Bowman <bowman@math.ualberta.ca> 1998]" << endl
		 << endl << "Usage: " << program << " [options] file1 [file2 ...]"
		 << endl;
}

void options()
{
	cerr << endl;
	cerr << "Options: " << endl;
	cerr << "-b\t\t single-byte (unsigned char instead of float) input"
		 << endl;
	cerr << "-f\t\t use a floating scale for each frame" << endl;
	cerr << "-g\t\t produce grey-scale output" << endl;
	cerr << "-h\t\t help" << endl;
	cerr << "-i\t\t invert vertical axis (y-origin at bottom)" << endl;
	cerr << "-m\t\t generate mpeg (.mpg) file" << endl;
	cerr << "-o option\t option to pass to convert" << endl;
	cerr << "-p\t\t preserve temporary output files" << endl;
	cerr << "-r\t\t remote X-server (substitute Postscript fonts)" << endl;
	cerr << "-v\t\t verbose output" << endl;
	cerr << "-l\t\t label frames with file names and values" << endl;
	cerr << "-F\t\t use a floating scale for each section" << endl;
	cerr << "-x mag\t\t overall magnification factor" << endl;
	cerr << "-H hmag\t\t horizontal magnification factor" << endl;
	cerr << "-V vmag\t\t vertical magnification factor" << endl;
	cerr << "-B begin\t first frame to process" << endl;
	cerr << "-E end\t\t last frame to process" << endl;
	cerr << "-L lower\t last section to process" << endl;
	cerr << "-U upper\t first section to process" << endl;
	cerr << "-S skip\t\t interval between processed frames" << endl;
	cerr << "-X xsize\t explicit horizontal size" << endl;
	cerr << "-Y ysize\t explicit vertical size" << endl;
	cerr << "-Z zsize\t explicit number of sections/frame" << endl;
	cerr << "-pointsize size\t point size to use with labels" << endl;
	cerr << "-avgx points\t number of data points per pixel in x direction " 
		 << "[default " << sx << "]" << endl;
	cerr << "-avgy points\t number of data points per pixel in y direction " 
		 << "[default " << sy << "]" << endl;
	cerr << "-symmetric\t make color palette symmetric about zero"
		 << " (if possible)" << endl;
	cerr << endl;
	cerr << "Available color palettes:" << endl;
	cerr << "-rainbow\t rainbow [default]" << endl;
	cerr << "-wheel\t\t full color wheel" << endl;
	cerr << "-idl39\t\t IDL palette 39" << endl;
	cerr << "-rgreyb\t\t red-grey-blue" << endl;
	cerr << endl;
	cerr << "Available transforms:" << endl;
	cerr << "-identity\t no transformation [default]" << endl;
	cerr << "-circle\t\t circle" << endl;
	cerr << "-torus\t\t projected torus (3D)" << endl;
	cerr << endl;
	cerr << "Available 3D options:" << endl;
	cerr << "-Rfactor factor\t viewpoint distance factor [default " << Rfactor
		 << "]" << endl;
	cerr << "-Theta radians\t angle to rotate about X axis [default " << Theta
		 << "]" << endl;
	cerr << "-Phi radians\t angle to rotate about Z axis [default " << Phi
		 << "]" << endl;
	cerr << "-yxaspect ratio\t Y/X aspect ratio [default " << yxaspect 
		 << "]" << endl;
	cerr << "-zxaspect ratio\t Z/X aspect ratio [default " << zxaspect 
		 << "]" << endl;
}

int main(int argc, char *const argv[])
{
	int nset=0, mx=1, my=1;
	int n,begin=0, skip=1, end=INT_MAX;
	int lower=0, upper=INT_MAX;
	int make_mpeg=0;
	u_char *red,*green,*blue;
	int trans=0;
	int palette=0;
	
	int syntax=0;
	extern int optind;
	extern char *optarg;
	int option_index = 0;
	
	static struct option long_options[] =
             {
               {"verbose", 0, 0, 'v'},
               {"rainbow", 0, &palette, RAINBOW},
               {"wheel", 0, &palette, WHEEL},
               {"idl39", 0, &palette, IDL39},
               {"rgreyb", 0, &palette, RGreyB},
               {"identity", 0, &trans, IDENTITY},
               {"circle", 0, &trans, CIRCLE},
               {"torus", 0, &trans, TORUS},
               {"symmetric", 0, &symmetric, 1},
               {"avgx", 1, 0, AVGX},
               {"avgy", 1, 0, AVGY},
               {"Nxfine", 1, 0, NXFINE},
               {"Nyfine", 1, 0, NYFINE},
               {"Nzfine", 1, 0, NZFINE},
               {"pointsize", 1, 0, POINTSIZE},
               {"Rfactor", 1, 0, RFACTOR},
               {"Theta", 1, 0, THETA},
               {"Phi", 1, 0, PHI},
               {"yxaspect", 1, 0, YXASPECT},
               {"zxaspect", 1, 0, ZXASPECT},
               {0, 0, 0, 0}
             };
	
#ifdef __GNUC__	
	optind=0;
#endif	
	errno=0;
	while (1) {
		int c = getopt_long_only(argc,argv,
								 "bfghbilmprvFo:x:H:V:B:E:L:U:S:X:Y:Z:",
								 long_options,&option_index);
		if (c == -1) break;
		switch (c) {
		case 0:
			break;
		case 'b':
			byte=1;
			break;
		case 'f':
			floating_scale=1;
			break;
		case 'F':
			floating_section=1;
			floating_scale=1;
			break;
		case 'g':
			grey=1;
			break;
		case 'h':
			usage(argv[0]);
			options();
			cerr << endl;
			exit(0);
		case 'i':
			invert=1;
			break;
		case 'm':
			make_mpeg=1;
			break;
		case 'o':
			option << " " << optarg;
			break;
		case 'p':
			preserve=1;
			break;
		case 'r':
			remote=1;
			break;
		case 'v':
			verbose=1;
			break;
		case 'l':
			label=1;
			break;
		case 'x':
			mx=my=atoi(optarg);
			break;
		case 'H':
			mx=atoi(optarg);
			break;
		case 'V':
			my=atoi(optarg);
			break;
		case 'B':
			begin=atoi(optarg);
			break;
		case 'E':
			end=atoi(optarg);
			break;
		case 'L':
			lower=atoi(optarg);
			break;
		case 'U':
			upper=atoi(optarg);
			break;
		case 'S':
			skip=atoi(optarg);
			if(skip <= 0) msg(ERROR,"skip value must be positive");
			break;
		case 'X':
			nx1=atoi(optarg);
			implicit=0;
			break;
		case 'Y':
			ny1=atoi(optarg);
			implicit=0;
			break;
		case 'Z':
			nz1=atoi(optarg);
			implicit=0;
			break;
		case AVGX:
			sx=atoi(optarg);
			break;
		case AVGY:
			sy=atoi(optarg);
			break;
		case NXFINE:
			Nxfine=atoi(optarg);
			break;
		case NYFINE:
			Nyfine=atoi(optarg);
			break;
		case NZFINE:
			Nzfine=atoi(optarg);
			break;
		case POINTSIZE:
			pointsize=atoi(optarg);
		    break;
		case RFACTOR:
			Rfactor=atof(optarg);
		    break;
		case THETA:
			Theta=atof(optarg);
		    break;
		case PHI:
			Phi=atof(optarg);
		    break;
		case YXASPECT:
			yxaspect=atof(optarg);
		    break;
		case ZXASPECT:
			zxaspect=atof(optarg);
		    break;
		default:
			syntax=1;
		}
	}
	
	errno=0;
	int nfiles=argc-optind;
	if(syntax || nfiles < 1) {
		usage(argv[0]);
		cerr << endl << "Type '" << argv[0]
			 << " -h' for a descriptions of options." << endl;
		exit(1);
	}
	
	option << ends;
	convertprog=(nfiles > 1 || label) ? "montage" : "convert";
	
	char *const *argf=argv+optind;
		
	if(!floating_scale) {
		vminf=new double[nfiles];
		vmaxf=new double[nfiles];
	}
	
	strstream dirbuf;
	rgbdir=getenv("RGB_DIR");
	if(rgbdir) dirbuf << rgbdir;
	else dirbuf << "/tmp/" << getenv("USER");
	unsigned int process=getpid();
	dirbuf << "/rgb." << process << ends;
	rgbdirbuf << dirbuf.str() << "/" << ends;
	strstream buf;
	buf << "mkdirhier " << dirbuf.str() << ends;
	char *cmd=buf.str();
	if(verbose) cout << cmd << endl;
	system(cmd);
	rgbdir=rgbdirbuf.str();
	
	char *const format=grey ? "gray" : "rgb";
	int PaletteMin=grey ? 0 : ColorPaletteMin[palette];
	int PaletteRange=grey ? 255 : (ColorPaletteMax[palette]-PaletteMin);
	
	red=Red[palette];
	green=Green[palette];
	blue=Blue[palette];
	
	for(int f=0; f < nfiles; f++) {
		char *fieldname=argf[f];
		ixstream xin;
		
		openfield(xin,fieldname,nx,ny,nz);
	
		double *vmink=new double [nz], *vmaxk=new double [nz];
		
		kmin=0;
		kmax=nz-1;
		if(upper < lower || lower < 0) 
			msg(ERROR, "Invalid section range (%d,%d)",lower,upper);

		if(kmin < lower) kmin=lower;
		if(kmax > upper) kmax=upper;
		
		int mpal=max(5,my);
		int msep=max(2,my);
		
		switch(trans) {
		case 0: 
			Nx=nx; Ny=ny;
			break;
		case 1: 
			a0=(int) (ceil(ny/twopi)+0.5);
			a=a0+nx;
			Nx=Ny=2*a+1;
			break;
		case 2:
			kmin=kmax=0;
			a0=(int) (ceil(ny/twopi)+0.5);
			a=a0+nx;
			R0=(int) (zxaspect*a+0.5);
			R=R0+a;
			Pz=R*Rfactor;
			Rp=R;
			Nx=Ny=(2*Rp+1)*5/4;
			break;
		}
		
		int nz0=kmax-kmin+1;
		Array2<Ivec> Index;
		
		if(trans) {
			Index.Allocate(Nx,Ny);
			transform[trans](Index);
			}
		
		xsize=mx*Nx;
		ysize=my*Ny*nz0+msep*nz0+mpal;
		
		Array3<float> value(nz,ny,nx);
		double gmin=DBL_MAX, gmax=-DBL_MAX; // Global min and max
		
		if(!floating_scale)	{
			n=0;
			int rc;
			int s=1;
			int set=0;
			do {
				double vmin,vmax;
				rc=readframe(xin,nx,ny,nz,value,vmin,vmax,vmink,vmaxk);
				if(rc == EOF) break;
				if(n < begin) continue;
				if(--s) continue;
				s=skip;
				if(vmin < gmin) gmin=vmin;
				if(vmax > gmax) gmax=vmax;
				set++;
			} while (n++ < end && rc == 0);
			nset=nset ? min(nset,set) : set;
			
			if(symmetric && gmin < 0 && gmax > 0) {
				gmax=max(-gmin,gmax);
				gmin=-gmax;
			}
			vminf[f]=gmin;
			vmaxf[f]=gmax;
			
			if(verbose && f==0) cout << nset << " frames found." << endl;
		}
		
		xin.close();
		openfield(xin,fieldname,nx,ny,nz);
		
		n=0;
		int rc;
		int s=1;
		int set=0;
		do {
			double vmin,vmax;
			rc=readframe(xin,nx,ny,nz,value,vmin,vmax,vmink,vmaxk);
			if(rc == EOF) break;
			if(n < begin) continue;
			if(--s) continue;
			s=skip;
			
			if(!floating_scale) {vmin=gmin; vmax=gmax;}
			
			strstream buf;
			buf << rgbdir << fieldname << setfill('0') << setw(NDIGITS)
				<< set << "." << format << ends;
			char *oname=buf.str();
			ofstream fout(oname);
			if(!fout) {
				cleanup();
				msg(ERROR,"Cannot open output file %s",oname);
			}
			
			for(int k=kmin; k <= kmax; k++) {
				if(floating_section) {vmin=vmink[k]; vmax=vmaxk[k];}
				double step=(vmax == vmin) ? 0.0 : PaletteRange/(vmax-vmin);
				for(int j=0; j < Ny; j++)  {
					for(int j2=0; j2 < my; j2++) {
						for(int i=0; i < Nx; i++)  {
							Ivec x;
							if(trans) x=Index(i,j);
							else x=Ivec(i,j,k);
							int index;
							
							if (step == 0.0 ||
								x.i == Undefined || 
								x.j == Undefined ||
								x.k == Undefined) index=-1;
							else {
								Real val;
								int ci=(int) floor(x.i);
								int cj=(int) floor(x.j);
								int ck=(int) floor(x.k);
								
								if(ci < -1 || ci >= nx ||
								   cj < -1 || cj >= ny ||
								   ck < -1 || ck >= nz)
									msg(ERROR,"Index out of range");
								
								int x1,x2,y1,y2,z1,z2;
								
								if(ci == -1 || ci == nx-1) {x1=nx-1; x2=0;}
								else {x1=ci; x2=ci+1;}
								
								if(cj == -1 || cj == ny-1) {y1=ny-1; y2=0;}
								else {y1=cj; y2=cj+1;}
								
								if(ck == -1 || ck == nz-1) {z1=nz-1; z2=0;}
								else {z1=ck; z2=ck+1;}
								
								val=(ck+1-x.k)*((cj+1-x.j)*((ci+1-x.i)*
															value(z1,y1,x1)+
															(x.i-ci)*
															value(z1,y1,x2))+
												(x.j-cj)*((ci+1-x.i)*
														  value(z1,y2,x1)+
														  (x.i-ci)*
														  value(z1,y2,x2)))+
									(x.k-ck)*((cj+1-x.j)*((ci+1-x.i)*
														  value(z2,y1,x1)+
														  (x.i-ci)*
														  value(z2,y1,x2))+
											  (x.j-cj)*((ci+1-x.i)*
														value(z2,y2,x1)+
														(x.i-ci)*
														value(z2,y2,x2)));


								index=(int) ((val-vmin)*step+0.5)+PaletteMin;
							}
							
							if(grey) {
								if(index == -1) 
									index=PaletteMin+PaletteRange; 
								for(int i2=0; i2 < mx; i2++)
									fout << (unsigned char) index;
							} else {
								unsigned char r=red[index+1],
									g=green[index+1], b=blue[index+1];
								for(int i2=0; i2 < mx; i2++)
									fout << r << g << b;
							}
						}
					}
				}
				unsigned char black=0;
				int xmsep=xsize*msep;
				for(int i=0; i < xmsep; i++) // Output separator
					if(grey) fout << black;
					else fout << black << black << black;
			}
			
			for(int j2=0; j2 < mpal; j2++) { // Output palette
				int Nxmx=Nx*mx;
				double step=1.0/Nxmx;
				for(int i=0; i < Nxmx; i++)  {
					int index;
					index=PaletteMin+(int) (PaletteRange*i*step+0.5);
					if(grey) fout << (unsigned char) index;
					else fout << red[index] << green[index] << blue[index];
				}
			}
			
			fout.close();
			if(!fout) {
				cleanup();
				msg(ERROR,"Cannot write to output file %s",oname);
			}
			set++;
		} while (n++ < end && rc == 0);
		nset=nset ? min(nset,set) : set;
		if(verbose && f==0 && floating_scale)
			cout << nset << " frames found." << endl;
	}
	
	if(nset == 1) msg(ERROR, "More than one frame required");
	
	if((remote || !label) && make_mpeg) putenv("DISPLAY=");

	if(nset) {
		if(label || make_mpeg) { 
			if(make_mpeg) montage(nfiles,argf,0,format,"miff");
			for(n=0; n < nset; n++) 
				montage(nfiles,argf,n,format,make_mpeg ? yuvformat : "miff");
			identify(nfiles,argf,0,"miff",xsize,ysize);
			
			if(make_mpeg) mpeg(nfiles,argf,nset-1,"mpg",xsize,ysize);
			else animate(nfiles,argf,nset-1,"miff","%d",xsize,ysize);
		} else {
			strstream buf;
			buf << "%0" << NDIGITS << "d" << ends;
			animate(nfiles,argf,nset-1,format,buf.str(),xsize,ysize);
		}
	}
	
	cleanup();
}

void cleanup()
{
	if(!preserve) {
		strstream buf;
		buf << "rm -r " << rgbdir << " > /dev/null 2>&1" << ends;
		char *cmd=buf.str();
		if(verbose) cout << cmd << endl;
		system(cmd);
	}
}

#if sun
char *separator="_______________________________________________";
#else
char *separator="                                               ";
#endif

void montage(int nfiles, char *const argf[], int n, char *const format,
			 char *const type)
{
	strstream buf;
	
	buf << convertprog << " -size " << xsize << "x" << ysize
		<< " -geometry " << xsize << "x" << ysize
		<< option.str() << " -interlace none ";
	if(pointsize) buf << "-pointsize " << pointsize << " ";
	for(int f=0; f < nfiles; f++) {
		char *fieldname=argf[f];
		if(label) {
			buf << " -label \"";
			if(!(floating_scale || byte)) 
				buf << setprecision(2) << vminf[f]
					<< separator << setprecision(2) << vmaxf[f] << "\\n";
			buf << fieldname << "\" ";
		}
		buf << format << ":" << rgbdir << fieldname
			<< setfill('0') << setw(NDIGITS) << n << "." << format << " ";
	}
	buf << yuvinterlace << type << ":" << rgbdir << argf[0] << n;
#if !NEW_IMAGEMAGIK	
	if(strcmp(type,"yuv3") != 0)
#endif
		buf << "." << type;
	if(!verbose) buf << "> /dev/null 2>&1";
	buf << ends;
	char *cmd=buf.str();
	if(verbose) cout << cmd << endl;
	system(cmd);
	
	if(n > 0 && !preserve) {
		strstream buf;
		buf << "rm ";
		for(int f=0; f < nfiles; f++) {
			char *fieldname=argf[f];
			buf << rgbdir << fieldname << setfill('0') << setw(NDIGITS)
				<< n << "." << format << " ";
		}
		if(!verbose) buf << " > /dev/null 2>&1";
		buf << ends;
		cmd=buf.str();
		if(verbose) cout << cmd << endl;
		system(cmd);
	}
}

void identify(int, char *const argf[], int n, char *const type,
			  int& xsize, int& ysize)
{
	strstream buf;
	char *iname=".identify";
	buf << "identify " << rgbdir << argf[0] << n << "." << type
		<< " > " << iname << ends;
	char *cmd=buf.str();
	if(verbose) cout << cmd << endl;
	system(cmd);
	ifstream fin(iname);
	if(!fin) {cleanup(); msg(ERROR,"Cannot open identify file %s",iname);}
	char c='x';
	while(fin && c != ' ') fin.get(c);
	fin >> xsize;
	fin.get(c);
	if(c != 'x') {
		cleanup();
		msg(ERROR,"Parse error in reading combined image size");
	}
	fin >> ysize;
}

void mpeg(int, char *const argf[], int n, char *const type,
		  int xsize, int ysize)
{
	strstream buf;
	buf << "mpeg -a 0 -b " << n << " -h " << xsize << " -v " << ysize
		<< " -PF " << rgbdir << argf[0] << " -s " << argf[0]
		<< "." << type;
	if(!verbose) buf << " > /dev/null";
	buf << ends;
	char *cmd=buf.str();
	if(verbose) cout << cmd << endl;
	system(cmd);
}

void animate(int, char *const argf[], int n, char *const type, 
			 const char *pattern, int xsize, int ysize)
{
	strstream buf;
	buf << "animate -scene 0-" << n << " -size " << xsize << "x" << ysize
		<< " " << type << ":" << rgbdir << argf[0] << pattern << "." << type
		<< ends;
	char *cmd=buf.str();
	if(verbose) cout << cmd << endl;
	system(cmd);
}

extern char **environ;

int system (char *command) {
	int pid, status;
	static int cleaning=0;

	if (command == 0) return 1;
	pid = fork();
	if (pid == -1) return -1;
	if (pid == 0) {
		char *argv[4];
		argv[0] = "sh";
		argv[1] = "-c";
		argv[2] = command;
		argv[3] = 0;
		execve("/bin/sh",argv,environ);
		exit(127);
	}
	
	do {
		if (waitpid(pid, &status, 0) == -1) {
			if (errno != EINTR) return -1;
		} else {
			if(status != 0) {
				if(cleaning) return status;
				cleaning=1; cleanup();
				msg(ERROR,"%s\nReceived signal %d",command,status);
			}
			return status;
		}
	} while(1);
}

void Circle(Array2<Ivec> Index)
{
	for(int u=0; u < Nx; u++)  {
		for(int v=0; v < Ny; v++)  {
			Index(u,v)=Ivec(Undefined,Undefined,Undefined);
		}
	}
	
	for(int u=0; u < Nx; u++)  {
		for(int v=0; v < Ny; v++)  {
			int i=u-a;
			int j=v-a;
			Real xi=sqrt(i*i+j*j)-a0;
			Real xj=(atan2(j,i)+pi)*ny/twopi;
			if(xi < 0 || xi > nx-1) xi=Undefined;
			if(xj > ny-1) xj -= ny;
			for(int k=kmin; k <= kmax; k++)
				Index(u,v)=Ivec(xi,xj,0.0);
		}
	}
}

void Torus(Array2<Ivec> Index)
{
	for(int u=0; u < Nx; u++)  {
		for(int v=0; v < Ny; v++)  {
			Index(u,v)=Ivec(Undefined,Undefined,Undefined);
		}
	}
	
	const Real xfinestep=1.0/(2.0*Nxfine+1.0);
	const Real yfinestep=1.0/(2.0*Nyfine+1.0);
	const Real zfinestep=1.0/(2.0*Nzfine+1.0);
	
	Real twopibyny=twopi/ny;
	Real twopibynz=twopi/nz;
	
	Array2<Real> zmax(Nx,Ny);
	for(int u=0; u < Nx; u++)
		for(int v=0; v < Ny; v++)
			zmax(u,v)=-REAL_MAX;
									
	Real cosPhi,sinPhi;
	Real cosTheta,sinTheta;
	sincos(-Phi,&sinPhi,&cosPhi);
	sincos(-Theta,&sinTheta,&cosTheta);
	
	// Rotate about z axis by -Phi, then about new x axis by -Theta.
	
	Real Axx=cosPhi;          Real Axy=-sinPhi;         Real Axz=0.0;
	Real Ayx=cosTheta*sinPhi; Real Ayy=cosTheta*cosPhi; Real Ayz=-sinTheta;
	Real Azx=sinTheta*sinPhi; Real Azy=sinTheta*cosPhi; Real Azz=cosPhi;
	
	Real cutoff=1.75*pi;
	
	Real xoffset=0.5+Rp*1.2;
	Real yoffset=0.5+Rp*1.5;
	
	for(int j=0; j < ny; j++)  {
		for(int j2=-Nyfine; j2 <= Nyfine; j2++)  {
			Real xj=j+j2*yfinestep;
			Real theta=xj*twopibyny;
			Real sintheta,costheta;
			sincos(theta,&sintheta,&costheta);
			for(int k=0; k < nz; k++)  {
				for(int k2=-Nzfine;	k2 <= Nzfine; k2++)  {
					Real xk=k+k2*zfinestep;
					Real phi=xk*twopibynz;
					if(phi >= 0 && phi <= cutoff) {
						Real cosphi,sinphi;
						sincos(phi,&sinphi,&cosphi);
						for(int i=0; i < nx; i++)  {
							for(int i2=(i == 0 ? 0 : -Nxfine);
								i2 <= (i == nx-1 ? 0 : Nxfine); i2++)  {
								Real xi=i+i2*xfinestep;
								Real r=a0+xi;
								Real rperp=R0+r*costheta;
								Real x=rperp*cosphi;
								Real y=rperp*sinphi;
								Real z=r*sintheta;
							
								Real xp=Axx*x+Axy*y+Axz*z;
								Real yp=Ayx*x+Ayy*y+Ayz*z;
								Real zp=Azx*x+Azy*y+Azz*z;
							
								Real projection=Pz/(Pz-zp);
								int u=(int)(xp*projection+xoffset);
								int v=Ny-((int)(yp*projection+yoffset));
								if(u >= 0 && u < Nx && v >=0 && v < Ny
								   && zp > zmax(u,v)) {
									zmax(u,v)=zp;
									Index(u,v)=Ivec(xi,xj,xk);
								}	
							}
						}
					}
				}
			}
		}
	}
	return;
}
