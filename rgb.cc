/* RGB:  A movie production utility
Copyright (C) 1997 John C. Bowman (bowman@ipp-garching.mpg.de)

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
const char VERSION[]="1.0";

#ifdef NEW_IMAGEMAGICK
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

#include "DynVector.h"
#include "rgb.h"

static double *vminf, *vmaxf;
static char *rgbdir;
static strstream rgbdirbuf;
static int xsize,ysize;

static int verbose=0;
static int floating_scale=0;
static int floating_section=0;
static int byte=0;
static int preserve=0;

int implicit=1;
int zero=0;
int invert=0;
int gray=0;

void cleanup();

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

inline float get_value(ifstream& fin)
{
	unsigned char c;
	fin.get(c);
	return (float) c;
}

inline float get_value(ixstream& fin)
{
	float v;
	fin >> v;
	return v;
}

template<class T>
int readframe(T& fin, int nx, int ny, int nz, float **value,
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
				float v=get_value(fin);
				if(fin.eof()) {
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
		if(byte) {
			char c;
			fin >> c;
		}
		int nx0,ny0,nz0;
		fin >> nx0 >> ny0 >> nz0;
		if(fin.eof()) return 1;
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

void usage(char *program)
{
	cerr << PROGRAM << " version " << VERSION
		 << " [(C) John C. Bowman <bowman@ipp-garching.mpg.de> 1997]" << endl
		 << endl << "Usage: " << program
		 << " [-bfFghilmpvz] [-x mag] [-H hmag] [-V vmag] " << endl
		 << "           [-B begin] [-E end] [-P palette] [-S skip]" << endl 
		 << "           [-X xsize -Y ysize [-Z zsize]] file1 [file2 ...]"
		 << endl;
}

void options()
{
	cerr << "Options: " << endl;
	cerr << "-b\t\t single-byte (unsigned char instead of float) input"
		 << endl;
	cerr << "-f\t\t use a floating scale for each frame" << endl;
	cerr << "-F\t\t use a floating scale for each section" << endl;
	cerr << "-g\t\t produce gray-scale output" << endl;
	cerr << "-h\t\t help" << endl;
	cerr << "-i\t\t invert vertical axis (y-origin at bottom)" << endl;
	cerr << "-l\t\t label frames with file names and values" << endl;
	cerr << "-m\t\t generate mpeg (.mpg) file" << endl;
	cerr << "-p\t\t preserve temporary output files" << endl;
	cerr << "-v\t\t verbose output" << endl;
	cerr << "-z\t\t make color palette symmetric about zero" <<
		" (if possible)" << endl;
	cerr << "-x mag\t\t overall magnification factor" << endl;
	cerr << "-H hmag\t\t horizontal magnification factor" << endl;
	cerr << "-V vmag\t\t vertical magnification factor" << endl;
	cerr << "-B begin\t first frame to process" << endl;
	cerr << "-E end\t\t last frame to process" << endl;
	cerr << "-P palette\t palette (integer between 0 and " << NPalette-1 <<
		")" << endl;
	cerr << "-S skip\t\t interval between processed frames" << endl;
	cerr << "-X xsize\t explicit horizontal size" << endl;
	cerr << "-Y ysize\t explicit vertical size" << endl;
	cerr << "-Z zsize\t explicit number of sections/frame" << endl;
}

int main(int argc, char *const argv[])
{
	int nx=1,ny=1,nz=1;
	int nset=0, mx=1, my=1;
	int n,begin=0, skip=1, end=INT_MAX;
	int label=0;
	int make_mpeg=0;
	int syntax=0;
	extern int optind;
	extern char *optarg;
	int palette=0;
	u_char *red,*green,*blue;
	
#ifdef __GNUC__	
	optind=0;
#endif	
	errno=0;
	while (1) {
		int c = getopt(argc,argv,"bfghilmpvzFx:H:V:B:E:P:S:X:Y:Z:");
		if (c == -1) break;
		switch (c) {
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
		case 'l':
			label=1;
			break;
		case 'm':
			make_mpeg=1;
			break;
		case 'p':
			preserve=1;
			break;
		case 'v':
			verbose=1;
			break;
		case 'z':
			zero=1;
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
		case 'P':
			palette=atoi(optarg);
			if(palette < 0 || palette >= NPalette)
				msg(ERROR,"Invalid palette (%d)",palette);
			break;
		case 'S':
			skip=atoi(optarg);
			if(skip <= 0) msg(ERROR,"skip value must be positive");
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
	
	if(nfiles > 1) label=1;
	
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
		ifstream fin;
		ixstream xin;
		
		if(byte) openfield(fin,fieldname,nx,ny,nz);
		else openfield(xin,fieldname,nx,ny,nz);
	
		double *vmink=new double [nz], *vmaxk=new double [nz];
		
		int mpal=max(5,my);
		int msep=max(2,my);
		xsize=mx*nx;
		ysize=my*ny*nz+msep*nz+mpal;
		
		float **value=new float* [nz];
		double gmin=DBL_MAX, gmax=-DBL_MAX; // Global min and max
		for(int k=0; k < nz; k++) value[k]=new float[nx*ny];
		
		if(!floating_scale)	{
			n=0;
			int rc;
			int s=skip;
			int set=0;
			do {
				double vmin,vmax;
				rc=byte ? readframe(fin,nx,ny,nz,value,vmin,vmax,vmink,vmaxk) :
					readframe(xin,nx,ny,nz,value,vmin,vmax,vmink,vmaxk);
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
		}
		
		if(byte) {fin.close(); openfield(fin,fieldname,nx,ny,nz);}
		else {xin.close(); openfield(xin,fieldname,nx,ny,nz);}
		
		n=0;
		int rc;
		int s=skip;
		int set=0;
		do {
			double vmin,vmax;
			rc=byte ? readframe(fin,nx,ny,nz,value,vmin,vmax,vmink,vmaxk) :
				readframe(xin,nx,ny,nz,value,vmin,vmax,vmink,vmaxk);
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
			
			for(int k=0; k < nz; k++) {
				if(floating_section) {vmin=vmink[k]; vmax=vmaxk[k];}
				double step=(vmax == vmin) ? 0.0 : PaletteRange/(vmax-vmin);
				for(int j=0; j < ny; j++)  {
					for(int j2=0; j2 < my; j2++) {
						for(int i=0; i < nx; i++)  {
							int index=(step == 0.0) ? PaletteRange/2 : 
								(int) ((value[k][i+nx*j]-vmin)*step+0.5);
							index += PaletteMin;
							if(gray) {
								for(int i2=0; i2 < mx; i2++)
									fout << (unsigned char) index;
							} else {
								unsigned char r=red[index],
									g=green[index], b=blue[index];
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
				int nxmx=nx*mx;
				double step=1.0/nxmx;
				for(int i=0; i < nxmx; i++)  {
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
	}
	
	if(nset) {
		if(label || make_mpeg) { 
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
	buf << "montage -size " << xsize << "x" << ysize
	    << " -geometry " << xsize << "x" << ysize << " -interlace none";
	for(int f=0; f < nfiles; f++) {
		char *fieldname=argf[f];
		buf << " -label \"";
		if(!(floating_scale || byte)) 
			buf << setprecision(2) << vminf[f]
				<< separator << setprecision(2) << vmaxf[f] << "\\n";
		buf << fieldname << "\" " << rgbdir
			<< fieldname << setfill('0') << setw(4) << n << "." << format;
	}
	buf << " " << yuvinterlace << type << ":" << rgbdir << argf[0] << n;
#ifndef NEW_IMAGEMAGIK	
	if(strcmp(type,"yuv3") != 0)
#endif
		buf << "." << type;
	if(!verbose) buf << "> /dev/null 2>&1";
	buf << ends;
	char *cmd=buf.str();
	if(verbose) cout << cmd << endl;
	system(cmd);
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
