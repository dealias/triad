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
const char VERSION[]="1.04";

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

#include "DynVector.h"
#include "Array.h"
#include "rgb.h"

#ifndef PI
#define PI 3.14159265358979323846264338327950288
#endif

const double pi=PI;
const double twopi=2.0*PI;

static double *vminf, *vmaxf;
static char *rgbdir;
static strstream rgbdirbuf;
static int xsize,ysize;

static int verbose=0;
static int floating_scale=0;
static int floating_section=0;
static int preserve=0;
static int remote=0;
static int pointsize=0;
static int gray=0;
static char *convertprog;
static strstream option;
static int Nx,Ny;
static int kmin,kmax;
static Real Pz;
static int a0,a;
static int R0,R,Rp;

int nx=1,ny=1,nz=1;
int byte=0;
int implicit=1;
int zero=0;
int invert=0;

Real alpha=1.0;
Real Rfactor=1.0;
Real Theta=0.0;
Real Phi=0.0;

void cleanup();
int system (char *command);

class Ivec {
public:
	int i;
	int j;
	int k;
	Ivec() {}
	Ivec(int i0, int j0, int k0) : i(i0), j(j0), k(k0) {}
};

typedef void Transform(Array3<Ivec>);
Transform Circle,ProjectedTorus;

Transform *transform[]={NULL,Circle,ProjectedTorus};
unsigned int Ntransform=sizeof(transform)/sizeof(Transform *);

template<class T>
void openfield(T& fin, char *fieldname, int& nx, int& ny, int& nz)
{
	fin.open(fieldname);
	if(!fin) {cleanup(); msg(ERROR,"Cannot open input file %s",fieldname);}
	if(implicit) {
		fin >> nx >> ny >> nz;
		if(fin.eof()) {
			cleanup(); msg(ERROR,"End of file during processing");
		}
	}
}

int readframe(ixstream& xin, int nx, int ny, int nz, float **value,
			  double& gmin, double& gmax, double *vmink, double *vmaxk)
{
	gmin=DBL_MAX; gmax=-DBL_MAX;
	double vmin=DBL_MAX, vmax=-DBL_MAX;
	
	errno=0;
	for(int k=0; k < nz; k++) {
		float *valuek=value[k];
		int start,stop,incr;
		if(invert) {
			start=ny-1;
			stop=-1;
			incr=-1;
		} else {
			start=0;
			stop=ny;
			incr=1;
		}
		
		for(int j=start; j != stop; j += incr) {
			int nxj=nx*j;
			for(int i=0; i < nx; i++) {
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
				if(v < vmin) vmin=v;
				if(v > vmax) vmax=v;
				valuek[nxj+i]=v;
			}
		}
		if(zero && vmin < 0 && vmax > 0) {
			vmax=max(-vmin,vmax);
			vmin=-vmax;
		}
		vmink[k]=vmin;
		vmaxk[k]=vmax;
		if(vmin < gmin) gmin=vmin;
		if(vmax > gmax) gmax=vmax;
	}
	
	if(zero && gmin < 0 && gmax > 0) {
		gmax=max(-gmin,gmax);
		gmin=-gmax;
	}
	
	if(implicit) {
		int nx0,ny0,nz0;
		xin >> nx0 >> ny0 >> nz0;
		if(xin.eof()) return 1;
		if(nx0 != nx || ny0 != ny || nz0 != nz) {
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
			 int xsize, int ysize);
void manimate(int argc, char *const argf[], int n, char *const type,
			  int xsize, int ysize);
	
extern "C" int getopt(int argc, char *const argv[], const char *optstring);
#ifdef __i386__
extern "C" void putenv(const char *);
#endif

void usage(char *program)
{
	cerr << PROGRAM << " version " << VERSION
		 << " [(C) John C. Bowman <bowman@math.ualberta.ca> 1998]" << endl
		 << endl << "Usage: " << program
		 << " [-bfghimprvzF] [-a ratio] [-c theta] [-d phi]" << endl
	     << "           [-l pointsize] [-o option]" << endl
		 << "           [-x mag] [-H hmag] [-V vmag]" << endl
		 << "           [-B begin] [-E end] [-L lower] [-U upper]" << endl
         << "           [-P palette] [-R factor]" << endl
		 << "           [-S skip] [-T transform]" << endl
		 << "           [-X xsize -Y ysize [-Z zsize]] file1 [file2 ...]"
		 << endl;
}

void options()
{
	cerr << "Options: " << endl;
	cerr << "-a ratio\t Z/X and Z/Y aspect ratio (for 3D plots)" << endl;
	cerr << "-b\t\t single-byte (unsigned char instead of float) input"
		 << endl;
	cerr << "-c theta\t radians to rotate about X axis (for 3D plots)" << endl;
	cerr << "-d phi\t radians to rotate about Z axis (for 3D plots)" << endl;
	cerr << "-f\t\t use a floating scale for each frame" << endl;
	cerr << "-g\t\t produce gray-scale output" << endl;
	cerr << "-h\t\t help" << endl;
	cerr << "-i\t\t invert vertical axis (y-origin at bottom)" << endl;
	cerr << "-m\t\t generate mpeg (.mpg) file" << endl;
	cerr << "-o option\t option to pass to convert" << endl;
	cerr << "-p\t\t preserve temporary output files" << endl;
	cerr << "-r\t\t remote X-server (substitute Postscript fonts)" << endl;
	cerr << "-v\t\t verbose output" << endl;
	cerr << "-z\t\t make color palette symmetric about zero"
		 << " (if possible)" << endl;
	cerr << "-l pointsize\t label frames with file names and values" << endl;
	cerr << "-x mag\t\t overall magnification factor" << endl;
	cerr << "-F\t\t use a floating scale for each section" << endl;
	cerr << "-H hmag\t\t horizontal magnification factor" << endl;
	cerr << "-V vmag\t\t vertical magnification factor" << endl;
	cerr << "-B begin\t first frame to process" << endl;
	cerr << "-E end\t\t last frame to process" << endl;
	cerr << "-L lower\t last section to process" << endl;
	cerr << "-U upper\t first section to process" << endl;
	cerr << "-P palette\t palette (integer between 0 and " << NPalette-1
		 << ")" << endl;
	cerr << "-R factor\t viewpoint distance factor (for 3D plots)" << endl;
	cerr << "-S skip\t\t interval between processed frames" << endl;
	cerr << "-T transform\t 2D transformation (integer between 0 and " 
		 << Ntransform-1 << ")" << endl;
	cerr << "-X xsize\t explicit horizontal size" << endl;
	cerr << "-Y ysize\t explicit vertical size" << endl;
	cerr << "-Z zsize\t explicit number of sections/frame" << endl;
}

int main(int argc, char *const argv[])
{
	int nset=0, mx=1, my=1;
	int n,begin=0, skip=1, end=INT_MAX;
	int lower=0, upper=INT_MAX;
	int make_mpeg=0;
	int syntax=0;
	extern int optind;
	extern char *optarg;
	int palette=0;
	u_char *red,*green,*blue;
	unsigned int trans=0;
	
#ifdef __GNUC__	
	optind=0;
#endif	
	errno=0;
	while (1) {
		int c = getopt(argc,argv,
					   "bfghimprvzFa:c:d:l:o:x:H:V:B:E:L:U:P:R:S:T:X:Y:Z:");
		if (c == -1) break;
		switch (c) {
		case 'a':
			alpha=atof(optarg);
		    break;
		case 'b':
			byte=1;
			break;
		case 'c':
			Theta=atof(optarg);
		    break;
		case 'd':
			Phi=atof(optarg);
		    break;
		case 'f':
			floating_scale=1;
			break;
		case 'F':
			floating_section=1;
			floating_scale=1;
			break;
		case 'g':
			gray=1;
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
		case 'z':
			zero=1;
			break;
		case 'l':
			pointsize=atoi(optarg);
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
		case 'P':
			palette=atoi(optarg);
			if(palette < 0 || palette >= NPalette)
				msg(ERROR,"Invalid palette (%d)",palette);
			break;
		case 'R':
			Rfactor=atof(optarg);
		    break;
		case 'S':
			skip=atoi(optarg);
			if(skip <= 0) msg(ERROR,"skip value must be positive");
			break;
		case 'T':
			trans=atoi(optarg);
			if(trans >= Ntransform)
				msg(ERROR, "Invalid transform (%d)", trans);
			break;
		case 'X':
			nx=atoi(optarg);
			implicit=0;
			break;
		case 'Y':
			ny=atoi(optarg);
			implicit=0;
			break;
		case 'Z':
			nz=atoi(optarg);
			implicit=0;
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
	convertprog=(nfiles > 1 || pointsize) ? "montage" : "convert";
	
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
	
	char *const format=gray ? "gray" : "rgb";
	int PaletteMin=gray ? 0 : ColorPaletteMin[palette];
	int PaletteRange=gray ? 255 : (ColorPaletteMax[palette]-PaletteMin);
	
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
			R0=(int) (alpha*a+0.5);
			R=R0+a;
			Pz=R*Rfactor;
			Rp=R;
			Nx=Ny=(2*Rp+1)*5/4;
			break;
		}
		
		int nz0=kmax-kmin+1;
		Array3<Ivec> Index;
		
		if(trans) {
			Index.Allocate(Nx,Ny,kmax-kmin+1);
			transform[trans](Index);
			}
		
		xsize=mx*Nx;
		ysize=my*Ny*nz0+msep*nz0+mpal;
		
		float **value=new float* [nz];
		double gmin=DBL_MAX, gmax=-DBL_MAX; // Global min and max
		for(int k=0; k < nz; k++) value[k]=new float[nx*ny];
		
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
			
			if(zero && gmin < 0 && gmax > 0) {
				gmax=max(-gmin,gmax);
				gmin=-gmax;
			}
			vminf[f]=gmin;
			vmaxf[f]=gmax;
			
			if(verbose && f==0) cout << nset << " frames found." << endl;
		}
		
		xin.close(); openfield(xin,fieldname,nx,ny,nz);
		
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
			buf << rgbdir << fieldname << setfill('0') << setw(4)
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
							if(trans) x=Index(i,j,k);
							else x=Ivec(i,j,k);
							int index=(step == 0.0 || 
									   x.i < 0 || x.i >= nx ||
									   x.j < 0 || x.j >= ny ||
								       x.k < 0 || x.k >= nz) ? 
								-1 : (int) ((value[x.k][x.i+nx*x.j]-vmin)*step
											+0.5)+PaletteMin;
							if(gray) {
								if(index == -1) 
									index=PaletteMin+PaletteRange/2; 
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
					if(gray) fout << black;
					else fout << black << black << black;
			}
			
			for(int j2=0; j2 < mpal; j2++) { // Output palette
				int Nxmx=Nx*mx;
				double step=1.0/Nxmx;
				for(int i=0; i < Nxmx; i++)  {
					int index;
					index=PaletteMin+(int) (PaletteRange*i*step+0.5);
					if(gray) fout << (unsigned char) index;
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
	
	if((remote || !pointsize) && make_mpeg) putenv("DISPLAY=");

	if(nset) {
		if(pointsize || make_mpeg) { 
			if(make_mpeg) montage(nfiles,argf,0,format,"miff");
			for(n=0; n < nset; n++) 
				montage(nfiles,argf,n,format,make_mpeg ? yuvformat : "miff");
			identify(nfiles,argf,0,"miff",xsize,ysize);
			
			if(make_mpeg) mpeg(nfiles,argf,nset-1,"mpg",xsize,ysize);
			else manimate(nfiles,argf,nset-1,"miff",xsize,ysize);
		} else animate(nfiles,argf,nset-1,format,xsize,ysize);
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
		if(pointsize) {
			buf << " -label \"";
			if(!(floating_scale || byte)) 
				buf << setprecision(2) << vminf[f]
					<< separator << setprecision(2) << vmaxf[f] << "\\n";
			buf << fieldname << "\" ";
		}
		buf << format << ":" << rgbdir << fieldname
			<< setfill('0') << setw(4) << n << "." << format << " ";
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
			buf << rgbdir << fieldname << setfill('0') << setw(4) << n << "."
				<< format << " ";
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

void animate(int, char *const argf[], int, char *const type,
			 int xsize, int ysize)
{
	strstream buf;
	buf << "animate -size " << xsize << "x" << ysize
		<< " -interlace none " << rgbdir << argf[0] << "*."	<< type << ends;
	char *cmd=buf.str();
	if(verbose) cout << cmd << endl;
	system(cmd);
}

void manimate(int, char *const argf[], int n, char *const type,
			  int xsize, int ysize)
{
	strstream buf;
	buf << "animate -scene 0-" << n << " -size " << xsize << "x" << ysize
		<< " " << type << ":" << rgbdir << argf[0] << "%d" << "." << type
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

void Circle(Array3<Ivec> Index)
{
	for(int i=0; i < Nx; i++)  {
		for(int j=0; j < Ny; j++)  {
			int i0=i-a;
			int j0=j-a;
			int ip=int(sqrt(i0*i0+j0*j0)-a0+0.5);
			int jp=int((atan2(j0,i0)+pi)*ny/twopi+0.5);
			if(jp == ny) jp=0;
			for(int k=kmin; k <= kmax; k++) Index(i,j,k)=Ivec(ip,jp,k);
		}
	}
}

void ProjectedTorus(Array3<Ivec> Index)
{
	for(int u=0; u < Nx; u++)  {
		for(int v=0; v < Ny; v++)  {
			Index(u,v,0)=Ivec(-1,-1,-1);
		}
	}
	
	const int Nxfine=4;
	const int Nyfine=5;
	const int Nzfine=100;

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
	
	// Rotation matrix; Rotate about x axis by Theta, then about z axis by Phi.
	Real Axx=cosPhi;          Real Axy=-sinPhi;         Real Axz=0.0;
	Real Ayx=cosTheta*sinPhi; Real Ayy=cosTheta*cosPhi; Real Ayz=-sinTheta;
	Real Azx=sinTheta*sinPhi; Real Azy=sinTheta*cosPhi; Real Azz=cosPhi;
	
	Real cutoff=1.75*pi;
	
	Real xoffset=0.5+Rp*1.2;
	Real yoffset=0.5+Rp*1.5;
	
	for(int j=0; j < ny; j++)  {
		for(int j2=-Nyfine; j2 <= Nyfine; j2++)  {
			Real theta=(j+j2*yfinestep)*twopibyny;
			Real sintheta,costheta;
			sincos(theta,&sintheta,&costheta);
			for(int k=0; k < nz; k++)  {
				for(int k2=-Nzfine; k2 <= Nzfine; k2++)  {
					Real phi=(k+k2*zfinestep)*twopibynz;
					if(phi >= 0 && phi <= cutoff) {
						Real cosphi,sinphi;
						sincos(phi,&sinphi,&cosphi);
						for(int i=0; i < nx; i++)  {
							Real a0i=a0+i;
							for(int i2=-Nxfine; i2 <= Nxfine; i2++)  {
								Real r=a0i+i2*xfinestep;
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
									Index(u,v,0)=Ivec(i,j,k);
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
