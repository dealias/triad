#include "xstream.h"
#include <iostream.h>
#include <limits.h>
#include "DynVector.h"
#include "rgb.h"
#include <errno.h>

void msg(char *text)
{
   cout << text << "." << endl;
   exit(FATAL);
}
	
int main(int argc, char *argv[])
{
	int nx,ny,nz;
	if(argc < 2) msg("File name required");
	ixstream fin(argv[1]);
	fin >> nx >> ny >> nz;
	int mx=(argc > 2) ? atoi(argv[2]) : 1;
	int my=(argc > 3) ? atoi(argv[3]) : mx;
	if(fin.eof()) msg("End of file during processing");
	
	int nxy=nx*ny;
	DynVector<float *> value;
	double vmin=DBL_MAX, vmax=-DBL_MAX;
	int n=0, l=0;
	while(1) {
		for(int k=0; k < nz; k++,l++) {
			value[l]=new float[nxy];
			for(int i=0; i < nxy; i++) {
				float v;
				fin >> v;
				if(fin.eof()) msg("End of file during processing");
				if(v < vmin) vmin=v;
				if(v > vmax) vmax=v;
				value[l][i]=v;
			}
		}
		n++;
		int nx0,ny0,nz0;
		fin >> nx0 >> ny0 >> nz0;
		if(fin.eof()) break;
		if(nx0 != nx || ny0 != ny || nz0 != nz) msg("Inconsistent image size");
	}
	
	int nset=n;
	int xsize=mx*nx;
	int ysize=my*((ny+1)*nz+1);
		
	char buf[200+2*strlen(argv[1])];
	sprintf(buf,"rm %s*.rgb >& /dev/null",argv[1]); // Delete old rgb files
	system(buf);

	char oname[20+strlen(argv[1])];
	double step=(vmax == vmin) ? 1.0 : 1.0/(vmax-vmin);
	l=0;
	for(n=0; n < nset; n++) {
		sprintf(oname,"%s%04d.rgb",argv[1],n);
		ofstream fout(oname);
		if(!fout) msg("Cannot open output file");
		for(int k=0; k < nz; k++,l++) {
			for(int j=0; j < ny; j++)  {
				for(int j2=0; j2 < my; j2++) {
					for(int i=0; i < nx; i++)  {
						int index=(int) (PaletteMax*(1.0-(
							value[l][i+nx*j]-vmin)*step)+0.5);
						unsigned char r=red[index], g=green[index],
							b=blue[index];
						for(int i2=0; i2 < mx; i2++) fout << r << g << b;
					}
				}
			}
			unsigned char black=0;
			for(int i=0; i < xsize*my; i++) // Output black separator
				fout << black << black << black;
		}
		for(int j2=0; j2 < my; j2++) { // Output palette
			for(int i=0; i < nx*mx; i++)  {
				int index=(int) (PaletteMax*(1.0-((double) i)/(nx*mx))+0.5);
				fout << red[index] << green[index] << blue[index];
			}
		}
		fout.close();
		if(!fout) msg("Cannot write to output file");
	}
	
	if(argc == 5)
		sprintf(buf,"convert -size %dx%d -interlace none %s*.rgb %s.mpg",
				xsize,ysize,argv[1],argv[1]);
	else
		sprintf(buf,"animate -size %dx%d -interlace none %s*.rgb",
				xsize,ysize,argv[1]);
	printf(buf);
	printf("\n");
	system(buf);
}
