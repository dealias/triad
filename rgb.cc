#include "xstream.h"
#include <iostream.h>
#include <limits.h>
#include <errno.h>
#include <strstream.h>
#include <iomanip.h>
#include <unistd.h>

#include "DynVector.h"
#include "rgb.h"

static float *vminf, *vmaxf;
static char *rgbdir;
static strstream rgbdirbuf;
static int xsize,ysize;

void cleanup(char *fieldname, char *type);
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

int verbose=0;
int floating_scale=0;
int byte=0;

extern "C" int getopt(int argc, char *const argv[], const char *optstring);

void usage(char *program)
{
	cerr << "Usage: " << program
		 << " [-bfghlmvz -x <mag> -H <hmag> -V <vmag>] file1 file2 ..."
		 << endl << endl;
}

int main(int argc, char *const argv[])
{
	int nx=1,ny=1,nz=1;
	int nset=0, mx=1, my=1;
	int c;
	int implicit=1;
	int gray=0;
	int label=0;
	int make_mpeg=0;
	int syntax=0;
	int zero=0;
	extern int optind;
	extern char *optarg;
	
#ifdef __GNUC__	
	optind=0;
#endif	
	while (1) {
		c = getopt(argc,argv,"bfghlmvzx:H:V:X:Y:Z:");
		if (c == -1) break;
		switch (c) {
		case 'b':
			byte=1;
			break;
		case 'f':
			floating_scale=1;
			break;
		case 'g':
			gray=1;
			break;
		case 'h':
			usage(argv[0]);
			cerr << "Options: " << endl;
			cerr << "-b\t\t unsigned byte (as opposed to float) input" << endl;
			cerr << "-f\t\t use a floating scale for each frame" << endl;
			cerr << "-g\t\t produce grey-scale output" << endl;
			cerr << "-h\t\t help" << endl;
			cerr << "-l\t\t label frames with file names and values" << endl;
			cerr << "-m\t\t generate .mpg mpeg file" << endl;
			cerr << "-v\t\t verbose output" << endl;
			cerr << "-z\t\t make color palette symmetric about zero" <<
				" (if possible)" << endl;
			cerr << "-x <mag>\t overall magnification factor" << endl;
			cerr << "-H <hmag>\t horizontal magnification factor" << endl;
			cerr << "-V <vmag>\t vertical magnification factor" << endl;
			cerr << "-X <xsize>\t use explicit horizontal size" << endl;
			cerr << "-Y <ysize>\t use explicit vertical size" << endl;
			cerr << "-Z <zsize>\t use explicit number of cross-sections/frame"
				 << endl;
			exit(0);
		case 'l':
			label=1;
			break;
		case 'm':
			make_mpeg=1;
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
		cerr << "Use '" << argv[0] << " -h' for a descriptions of options."
			 << endl;
		exit(1);
	}
	if(nfiles > 1) label=1;
	char *const *argf=argv+optind;
		
	vminf=new float[nfiles];
	vmaxf=new float[nfiles];
	
	rgbdir=getenv("RGB_DIR");
	if(!rgbdir) {
		rgbdirbuf << "/tmp/" << getenv("USER") << ends;
		rgbdir=rgbdirbuf.str();
	}
	
	char *const format=gray ? "gray" : "rgb";
	
	for(int f=0; f < nfiles; f++) {
		char *fieldname=argf[f];
		ixstream fin(fieldname);
		if(!fin) msg(ERROR,"Cannot open input file %s",fieldname);
		if(implicit) fin >> nx >> ny >> nz;
		if(fin.eof()) msg(ERROR,"End of file during processing");
	
		int nxy=nx*ny;
		DynVector<float *> value;
		DynVector<float> vminn,vmaxn;
		float gmin=FLT_MAX, gmax=-FLT_MAX; // Global min and max
		int n=0, l=0;
		while(1) {
			float vmin=FLT_MAX, vmax=-FLT_MAX;
			for(int k=0; k < nz; k++,l++) {
				float *valuel;
				value[l]=valuel=new float[nxy];
				for(int i=0; i < nxy; i++) {
					float v;
					if(byte) {
						unsigned char b;
						fin >> b;
						v=(float) b;
					} else fin >> v;
					if(fin.eof()) {
						if(implicit || i > 0) 
							msg(WARNING,"End of file during processing");
						break;
					}
					if(v < vmin) vmin=v;
					if(v > vmax) vmax=v;
					valuel[i]=v;
				}
				if(fin.eof()) break;
			}
			if(zero && vmin < 0 && vmax > 0) {
				vmax=max(-vmin,vmax);
				vmin=-vmax;
			}
			vminn[n]=vmin;
			vmaxn[n]=vmax;
			n++;
		
			if(vmin < gmin) gmin=vmin;
			if(vmax > gmax) gmax=vmax;
			if(implicit) {
				int nx0,ny0,nz0;
				fin >> nx0 >> ny0 >> nz0;
				if(fin.eof()) break;
				if(nx0 != nx || ny0 != ny || nz0 != nz)
					msg(ERROR,"Inconsistent image size");
			}
		}
	
		if(zero && gmin < 0 && gmax > 0) {
			gmax=max(-gmin,gmax);
			gmin=-gmax;
		}
			
		vminf[f]=gmin;
		vmaxf[f]=gmax;
		nset=n;
		xsize=mx*nx;
		int mpal=max(5,my);
		int msep=max(2,my);
		ysize=my*ny*nz+msep*nz+mpal;
		
		cleanup(fieldname,gray ? "*.gray" : "*.rgb");
		cleanup(fieldname,".*.*");
		
		int PaletteMax=gray ? 255 : ColorPaletteMax;
		l=0;
		double vmin=gmin;
		double vmax=gmax;
		double step=(vmax == vmin) ? 0.0 : PaletteMax/(vmax-vmin);
		
		for(n=0; n < nset; n++) {
			if(floating_scale) {
				vmin=vminn[n];
				vmax=vmaxn[n];
				step=(vmax == vmin) ? 0.0 : PaletteMax/(vmax-vmin);
			}
			
			strstream buf;
			buf << rgbdir << "/" << fieldname << setfill('0') << setw(4)
				<< n << "." << format << ends;
			char *oname=buf.str();
			ofstream fout(oname);
			if(!fout) msg(ERROR,"Cannot open output file %s",oname);
			
			for(int k=0; k < nz; k++,l++) {
				for(int j=0; j < ny; j++)  {
					for(int j2=0; j2 < my; j2++) {
						for(int i=0; i < nx; i++)  {
							int index=(step == 0.0) ? PaletteMax/2 : 
								(int) ((vmax-value[l][i+nx*j])*step+0.5);
							if(gray) {
								for(int i2=0; i2 < mx; i2++)
									fout << (unsigned char) index;
							} else {
								unsigned char
								r=red[index], g=green[index], b=blue[index];
								for(int i2=0; i2 < mx; i2++)
									fout << r << g << b;
							}
						}
					}
				}
				unsigned char black=0;
				for(int i=0; i < xsize*msep; i++) // Output separator
					if(gray) fout << black;
					else fout << black << black << black;
			}
			
			for(int j2=0; j2 < mpal; j2++) { // Output palette
				for(int i=0; i < nx*mx; i++)  {
					int index;
					index=(int) (PaletteMax*(1.0-((double) i)/(nx*mx))+0.5);
					if(gray) fout << (unsigned char) index;
					else fout << red[index] << green[index] << blue[index];
				}
			}
			
			fout.close();
			if(!fout) msg(ERROR,"Cannot write to output file %s",oname);
		}
	}	
	
	if(!label) animate(nfiles,argf,nset-1,format,xsize,ysize);
	else {
		if(make_mpeg) montage(nfiles,argf,0,format,"miff");
		for(int n=0; n < nset; n++) 
			montage(nfiles,argf,n,format,make_mpeg ? "yuv3" : "miff");
		identify(nfiles,argf,0,"miff",xsize,ysize);
		
		if(make_mpeg) mpeg(nfiles,argf,nset-1,"mpg",xsize,ysize);
		else manimate(nfiles,argf,nset-1,"miff",xsize,ysize);
	}
}

void cleanup(char *fieldname, char *type)
{
	strstream buf;
// Delete old files
	buf << "rm " << rgbdir << "/" << fieldname << type
		<< "> /dev/null 2>&1" << ends;
	char *cmd=buf.str();
	system(cmd);
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
	buf << "montage -size " << xsize << "x" << ysize << " -geometry "
		<< xsize << "x" << ysize << " -interlace none";
	for(int f=0; f < nfiles; f++) {
		char *fieldname=argf[f];
		buf << " -label \"";
		if(!floating_scale) 
			buf << setprecision(2) << vminf[f]
				<< separator << setprecision(2) << vmaxf[f] << "\\n";
		buf << fieldname << "\" " << rgbdir << "/"
			<< fieldname << setfill('0') << setw(4) << n << "." << format;
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
		<< " -interlace none " << rgbdir << "/" << argf[0] << "*."
		<< type << ends;
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
