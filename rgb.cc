#include "xstream.h"
#include <iostream.h>
#include <limits.h>
#include <errno.h>
#include <string.h>
#include <strstream.h>
#include <iomanip.h>
#include <unistd.h>

#include "DynVector.h"
#include "rgb.h"

static double *vmin;
static double *vmax;
static char *rgbdir;
static int xsize,ysize;

void montage(int argc, char *const argf[], int n, char *const type);
void identify(int argc, char *const argf[], int n, char *const type,
			  int& xsize, int& ysize);
void mpeg(int argc, char *const argf[], int n, char *const type,
		  int xsize, int ysize);
void animate(int argc, char *const argf[], int n, char *const type,
			 int xsize, int ysize);
void manimate(int argc, char *const argf[], int n, char *const type,
			  int xsize, int ysize);

int verbose=0;

extern "C" int getopt(int argc, char *const argv[], const char *optstring);

int main(int argc, char *const argv[])
{
	int nx,ny,nz;
	int nset=0, mx=1, my=1;
	int c;
	int label=0;
	int make_mpeg=0;
	extern int optind;
	extern char *optarg;
	
	optind=0;
	while (1) {
		c = getopt(argc,argv,"lmvx:H:V:");
		if (c == -1) break;
		switch (c) {
		case 'l':
			label=1;
			break;
		case 'm':
			make_mpeg=1;
			break;
		case 'v':
			verbose=1;
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
		default:
			cerr << "Usage: " << argv[0]
				 << "[-lmv -x <mag> -H <hmag> -V <vmag>] file1 file2 ..." 
				 << endl;
			exit(1);
		}
	}
	
	int nfiles=argc-optind;
	if(nfiles < 1) msg(ERROR,"File name required");
	if(nfiles > 1) label=1;
	char *const *argf=argv+optind;
		
	vmin=new double[nfiles];
	vmax=new double[nfiles];
	
	rgbdir=getenv("RGB_DIR");
	if(!rgbdir) rgbdir=".";
	
	for(int f=0; f < nfiles; f++) {
		char *fieldname=argf[f];
		ixstream fin(fieldname);
		if(!fin) msg(ERROR,"Cannot open input file %s",fieldname);
		fin >> nx >> ny >> nz;
		if(fin.eof()) msg(ERROR,"End of file during processing");
	
		int nxy=nx*ny;
		DynVector<float *> value;
		double vminf=DBL_MAX, vmaxf=-DBL_MAX;
		int n=0, l=0;
		while(1) {
			for(int k=0; k < nz; k++,l++) {
				value[l]=new float[nxy];
				for(int i=0; i < nxy; i++) {
					float v;
					fin >> v;
					if(fin.eof()) {
						msg(WARNING,"End of file during processing");
						break;
					}
					if(v < vminf) vminf=v;
					if(v > vmaxf) vmaxf=v;
					value[l][i]=v;
				}
				if(fin.eof()) break;
			}
			n++;
			int nx0,ny0,nz0;
			fin >> nx0 >> ny0 >> nz0;
			if(fin.eof()) break;
			if(nx0 != nx || ny0 != ny || nz0 != nz)
				msg(ERROR,"Inconsistent image size");
		}
	
		vmin[f]=vminf;
		vmax[f]=vmaxf;
		nset=n;
		xsize=mx*nx;
		int mpal=max(5,my);
		int msep=max(2,my);
		ysize=my*ny*nz+msep*nz+mpal;
		
		char *buf=new char[200+2*strlen(fieldname)];
// Delete old rgb files
		sprintf(buf,"rm %s/%s*.rgb > /dev/null 2>&1",rgbdir,fieldname);
		system(buf);

		char *oname=new char[20+strlen(fieldname)];
		double step=(vmaxf == vminf) ? 0.0 : 1.0/(vmaxf-vminf);
		step *= PaletteMax;
		l=0;
		for(n=0; n < nset; n++) {
			sprintf(oname,"%s/%s%04d.rgb",rgbdir,fieldname,n);
			ofstream fout(oname);
			if(!fout) msg(ERROR,"Cannot open output file %s",oname);
			for(int k=0; k < nz; k++,l++) {
				for(int j=0; j < ny; j++)  {
					for(int j2=0; j2 < my; j2++) {
						for(int i=0; i < nx; i++)  {
							int index=(step == 0.0) ? PaletteMax/2 : 
								(int) ((vmaxf-value[l][i+nx*j])*step+0.5);
							unsigned char 
								r=red[index], g=green[index], b=blue[index];
							for(int i2=0; i2 < mx; i2++) fout << r << g << b;
						}
					}
				}
				unsigned char black=0;
				for(int i=0; i < xsize*msep; i++) // Output separator
					fout << black << black << black;
			}
			for(int j2=0; j2 < mpal; j2++) { // Output palette
				for(int i=0; i < nx*mx; i++)  {
					int index;
					index=(int) (PaletteMax*(1.0-((double) i)/(nx*mx))+0.5);
					fout << red[index] << green[index] << blue[index];
				}
			}
			fout.close();
			if(!fout) msg(ERROR,"Cannot write to output file %s",oname);
		}
	}	
	
	if(!label) animate(nfiles,argf,nset-1,"yuv3",xsize,ysize);
	else {
		if(make_mpeg) montage(nfiles,argf,0,"miff");
		for(int n=0; n < nset; n++) 
			montage(nfiles,argf,n,make_mpeg ? "yuv3" : "miff");
		identify(nfiles,argf,0,"miff",xsize,ysize);
		
		if(make_mpeg) mpeg(nfiles,argf,nset-1,"mpg",xsize,ysize);
		else manimate(nfiles,argf,nset-1,"miff",xsize,ysize);
	}
}

void montage(int nfiles, char *const argf[], int n, char *const type)
{
	strstream buf;
	buf << "montage -size " << xsize << "x" << ysize << " -geometry "
		<< xsize << "x" << ysize << " -interlace none";
	for(int f=0; f < nfiles; f++) {
		char *fieldname=argf[f];
		buf << " -label \"" << setprecision(2) << vmin[f]
			<< "                       " << setprecision(2) << vmax[f] << "\\n"
			<< fieldname <<	"\" " << rgbdir << "/"
			<< fieldname << setfill('0') << setw(4) << n << ".rgb";
	}
	buf << " " << type << ":" << rgbdir << "/" << argf[0];
	if(type != "yuv3") buf << "." << type;
	buf << "." << n << ends;
	char *cmd=buf.str();
	if(verbose) cout << cmd << endl;
	system(cmd);
}

void identify(int, char *const argf[], int n, char *const type,
			  int& xsize, int& ysize)
{
	strstream buf;
	char *iname=".identify";
	buf << "identify " << rgbdir << "/" << argf[0] << "." << type << "." << n 
		<< " > " << iname  
		<< ends;
	char *cmd=buf.str();
	if(verbose) cout << cmd << endl;
	system(cmd);
	ifstream fin(iname);
	if(!fin) msg(ERROR,"Cannot open identify file %s",iname);
	char c='x';
	while(fin && c != ' ') fin.get(c);
	fin >> xsize;
	fin.get(c);
	if(c != 'x') msg(ERROR,"Parse error in reading combined image size");
	fin >> ysize;
}

void mpeg(int, char *const argf[], int n, char *const type,
		  int xsize, int ysize)
{
	strstream buf;
	buf << "mpeg -a 0 -b " << n << " -h " << xsize << " -v " << ysize
		<< " -PF " << rgbdir << "/" << argf[0] << "." << " -s " << argf[0]
		<< "." << type << ends;
	char *cmd=buf.str();
	if(verbose) cout << cmd << endl;
	system(cmd);
}

void animate(int, char *const argf[], int, char *const, int xsize, int ysize)
{
	strstream buf;
	buf << "animate -size " << xsize << "x" << ysize
		<< " -interlace none " << rgbdir << "/" << argf[0] << "*.rgb" << ends;
	char *cmd=buf.str();
	if(verbose) cout << cmd << endl;
	system(cmd);
}

void manimate(int, char *const argf[], int n, char *const type,
			  int xsize, int ysize)
{
	strstream buf;
	buf << "animate -scene 0-" << n << " -size " << xsize << "x" << ysize <<
		" " << type << ":" << rgbdir << "/" << argf[0] << "." << type
		<< ".%d" << ends;
	char *cmd=buf.str();
	if(verbose) cout << cmd << endl;
	system(cmd);
}
