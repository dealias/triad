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

template<class T>
void openfield(T& fin, char *fieldname, int& nx, int& ny, int& nz)
{
	fin.open(fieldname);
	if(!fin) msg(ERROR,"Cannot open input file %s",fieldname);
	if(implicit) {
		fin >> nx >> ny >> nz;
		if(fin.eof()) msg(ERROR,"End of file during processing");
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
			  double& vmin, double& vmax)
{
	vmin=DBL_MAX, vmax=-DBL_MAX;
	int nxy=nx*ny;
	for(int k=0; k < nz; k++) {
		float *valuek=value[k];
		for(int i=0; i < nxy; i++) {
			float v=get_value(fin);
			if(fin.eof()) {
				if(implicit || i > 0) 
					msg(WARNING,"End of file during processing");
				return EOF;
			}
			if(v < vmin) vmin=v;
			if(v > vmax) vmax=v;
			valuek[i]=v;
		}
	}
	if(zero && vmin < 0 && vmax > 0) {
		vmax=max(-vmin,vmax);
		vmin=-vmax;
	}
	if(implicit) {
		int nx0,ny0,nz0;
		fin >> nx0 >> ny0 >> nz0;
		if(fin.eof()) return 1;
		if(nx0 != nx || ny0 != ny || nz0 != nz)
			msg(ERROR,"Inconsistent image size");
	}
	return 0;
}

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
int implicit=1;
int zero=0;

extern "C" int getopt(int argc, char *const argv[], const char *optstring);

void usage(char *program)
{
	cerr << "Usage: " << program
		 << " [-bfghlmvz] [-x mag] [-H hmag] [-V vmag] [-B beg] [-E end]"
		 << endl 
		 << "           [-X xsize -Y ysize [-Z zsize]] file1 [file2 ...]"
		 << endl << endl;
}

int main(int argc, char *const argv[])
{
	int nx=1,ny=1,nz=1;
	int nset=0, mx=1, my=1;
	int n,begin=0, end=INT_MAX;
	int gray=0;
	int label=0;
	int make_mpeg=0;
	int syntax=0;
	extern int optind;
	extern char *optarg;
	
#ifdef __GNUC__	
	optind=0;
#endif	
	while (1) {
		char c = getopt(argc,argv,"bfghlmvzx:H:V:B:E:X:Y:Z:");
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
			cerr << "-b\t\t single-byte (unsigned char instead of float) input"
				 << endl;
			cerr << "-f\t\t use a floating scale for each frame" << endl;
			cerr << "-g\t\t produce grey-scale output" << endl;
			cerr << "-h\t\t help" << endl;
			cerr << "-l\t\t label frames with file names and values" << endl;
			cerr << "-m\t\t generate mpeg (.mpg) file" << endl;
			cerr << "-v\t\t verbose output" << endl;
			cerr << "-z\t\t make color palette symmetric about zero" <<
				" (if possible)" << endl;
			cerr << "-x mag\t\t overall magnification factor" << endl;
			cerr << "-H hmag\t\t horizontal magnification factor" << endl;
			cerr << "-V vmag\t\t vertical magnification factor" << endl;
			cerr << "-B beg\t\t first frame to process" << endl;
			cerr << "-E end\t\t last frame to process" << endl;
			cerr << "-X xsize\t explicit horizontal size" << endl;
			cerr << "-Y ysize\t explicit vertical size" << endl;
			cerr << "-Z zsize\t explicit number of cross-sections/frame"
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
		case 'B':
			begin=atoi(optarg);
			break;
		case 'E':
			end=atoi(optarg);
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
		cerr << "Type '" << argv[0] << " -h' for a descriptions of options."
			 << endl;
		exit(1);
	}
	if(nfiles > 1) label=1;
	char *const *argf=argv+optind;
		
	if(!floating_scale) {
		vminf=new double[nfiles];
		vmaxf=new double[nfiles];
	}
	
	rgbdir=getenv("RGB_DIR");
	if(!rgbdir) {
		rgbdirbuf << "/tmp/" << getenv("USER") << ends;
		rgbdir=rgbdirbuf.str();
	}
	
	char *const format=gray ? "gray" : "rgb";
	int PaletteMax=gray ? 255 : ColorPaletteMax;
		
	
	for(int f=0; f < nfiles; f++) {
		char *fieldname=argf[f];
		ifstream fin;
		ixstream xin;
		
		if(byte) openfield(fin,fieldname,nx,ny,nz);
		else openfield(xin,fieldname,nx,ny,nz);
	
		int mpal=max(5,my);
		int msep=max(2,my);
		xsize=mx*nx;
		ysize=my*ny*nz+msep*nz+mpal;
		
		float **value=new float* [nz];
		double gmin=DBL_MAX, gmax=-DBL_MAX; // Global min and max
		for(int k=0; k < nz; k++) value[k]=new float[nx*ny];
		
		cleanup(fieldname,gray ? "*.gray" : "*.rgb");
		cleanup(fieldname,".*.*");
		
		if(!floating_scale)	{
			n=0;
			int rc;
			do {
				double vmin,vmax;
				rc=byte ? readframe(fin,nx,ny,nz,value,vmin,vmax) :
					readframe(xin,nx,ny,nz,value,vmin,vmax);
				if(rc == EOF) break;
				if(vmin < gmin) gmin=vmin;
				if(vmax > gmax) gmax=vmax;
				n++;
			} while (rc == 0 && n < end);
			
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
		do {
			double vmin,vmax;
			rc=byte ? readframe(fin,nx,ny,nz,value,vmin,vmax) :
				readframe(xin,nx,ny,nz,value,vmin,vmax);
			if(rc == EOF) break;
			if(n < begin) continue;
			
			if(!floating_scale) {vmin=gmin;	vmax=gmax;}
			double step=(vmax == vmin) ? 0.0 : PaletteMax/(vmax-vmin);
			
			strstream buf;
			buf << rgbdir << "/" << fieldname << setfill('0') << setw(4)
				<< n << "." << format << ends;
			char *oname=buf.str();
			ofstream fout(oname);
			if(!fout) msg(ERROR,"Cannot open output file %s",oname);
			
			for(int k=0; k < nz; k++) {
				for(int j=0; j < ny; j++)  {
					for(int j2=0; j2 < my; j2++) {
						for(int i=0; i < nx; i++)  {
							int index=(step == 0.0) ? PaletteMax/2 : 
								(int) ((vmax-value[k][i+nx*j])*step+0.5);
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
			n++;
		} while (rc == 0 && n < end);
		nset=nset ? min(nset,n-begin) : n-begin;
	}
	
	if(label || make_mpeg) { 
		if(make_mpeg) montage(nfiles,argf,0,format,"miff");
		for(n=0; n < nset; n++) 
			montage(nfiles,argf,n,format,make_mpeg ? "yuv3" : "miff");
		identify(nfiles,argf,0,"miff",xsize,ysize);
		
		if(make_mpeg) mpeg(nfiles,argf,nset-1,"mpg",xsize,ysize);
		else manimate(nfiles,argf,nset-1,"miff",xsize,ysize);
	} else
		animate(nfiles,argf,nset-1,format,xsize,ysize);
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
