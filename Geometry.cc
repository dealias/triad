#include "NWave.h"

NWave *GeometryProblem;
GeometryBase *Geometry;

Weight *weightBase;

int GeometryCompare(const void *a, const void *b)
{
	return GeometryProblem->GeometryTable->DefaultCompare(a,b);
}

int GeometryKeyCompare(const void *key, const void *p, const size_t n)
{
	return GeometryProblem->GeometryTable->DefaultKeyCompare(key,p,n);
}

static int formatted;

WeightIndex WeightN;

int get_weights(DynVector<Weight>& weight, int *Nweight, char *filename,
				char *filenamef)
{
	ifstream fin;
	int complete;
	int n=0;
	char *ufilename=filename;
	
	formatted=0;
	fin.open(filename);
	if(!fin) {
		errno=0;
		filename=filenamef;
		fin.open(filename);
		if(fin) formatted=1;
	}
	
	if(fin) {
		if(formatted) fin >> n;
		else fin.read((char *) &n,sizeof(int));
		if(fin.eof()) fin.close();
		if(n == 0) testlock();
	}
	
	complete=(n ? 1 : 0);
		
	if(fin) {
		cout << newl << "READING WEIGHT FACTORS FROM " << filename << "."
			 << endl;
		if(n) weight.Resize(n);
		if(formatted) {
			for(int i=0; i < n; i++) fin >> weight[i];
			save_weights(weight,n,ufilename);
		} else {
			if(complete) fin.read((char *) weight.Base(),n*sizeof(Weight));
			else {
				while(fin.read((char *) &weight[n],sizeof(Weight))) n++;
				if(n) cout << n << " WEIGHT FACTORS READ: ("
						   << weight[n-1].Index() << ")/(" << WeightN << ")."
						   << endl;
			}
		}
		fin.close();
		errno=0;
		if((formatted || complete) && !fin.good())
			msg(ERROR,"Error reading from weight file %s",filename);
	}
	*Nweight=n;
	return complete;
}

char *write_error="Error writing to weight file %s";

void save_weights(DynVector<Weight>& w, int n, char *filename)
{
	ofstream fout;
	
	fout.open(filename);
	fout.write((char *) &n,sizeof(int));
	fout.write((char *) w.Base(),n*sizeof(Weight));
	fout.close();
	if(!fout.good()) msg(ERROR,write_error,filename);
}

void save_formatted_weights(DynVector<Weight>& w, int n, char *filename)
{
	if(formatted) return;
	
	ofstream fout;
	fout.open(filename);
	fout.precision(FLT_DIG);
	if(!fout) msg(ERROR,"Weight file %s could not be opened",filename);
	fout << n << newl;
	for(int i=0; i < n; i++) fout << w[i] << newl;
	fout.close();
	if(!fout.good()) msg(ERROR,write_error,filename);
}	

int out_weights(ofstream& fout, Weight* w, int lastindex, int n)
{
	lock();
	fout.write((char *) (w+lastindex),(n-lastindex)*sizeof(Weight));
	fout.flush();
	unlock();
	return n;
}
