#include "options.h"
#include "NWave.h"

GeometryBase *Geometry;

Weight *weightBase;
WeightIndex WeightN;

const int maxbins=65536;

int get_weights(DynVector<Weight>& weight, int *Nweight, char *filename)
{
	int complete;
	int n=0;
	
	ixstream fin(filename);
	if(fin) {
		fin >> n;
		if(fin.eof()) {fin.close(); errno=0;}
		if(n == 0) testlock();
	}
	
	complete=(n ? 1 : 0);
		
	if(fin) {
		cout << newl << "READING WEIGHT FACTORS FROM " << filename << "."
			 << endl;
		if(complete) {
			weight.Resize(n);
			for(int i=0; i < n; i++) fin >> weight[i];
			if(fin.bad()) 
				msg(ERROR,"Error reading from weight file %s",filename);
		} else {
			Weight w;
			while(fin >> w, !fin.eof()) weight[n++]=w;
			if(n) cout << n << " WEIGHT FACTORS READ: ("
					   << weight[n-1].Index() << ")/(" << WeightN << ")."
					   << endl;
		}
	}
	*Nweight=n;
	return complete;
}

void save_weights(DynVector<Weight>& w, int n, char *filename)
{
	oxstream fout(filename);
	if(!fout) msg(ERROR,"Could not open weight file %s",filename);
	fout << n << newl;
	for(int i=0; i < n; i++) fout << w[i] << newl;
	fout.close();
	if(!fout) msg(ERROR,"Cannot write to weight file %s",filename);
}	

int out_weights(oxstream& fout, Weight *w, int lastindex, int n)
{
	lock();
	for(int i=lastindex; i < n; i++) fout << w[i] << newl;
	fout.flush();
	unlock();
	return n;
}
