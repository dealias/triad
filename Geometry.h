#ifndef __Geometry_h__
#define __Geometry_h__ 1

#include <errno.h>

#include "utils.h"
#include "Pair.h"
#include "Bin.h"

#define GEOMETRY(key) {new Entry<Partition<key>,GeometryBase>\
						   (#key,GeometryTable);}

class GeometryBase {
protected:	
	int Nmode;	// number of unreflected (explicitly evolved) bins
public:	
	virtual char *Name()=0;
	virtual int Create()=0;
	virtual Nu BinAveragedLinearity(int)=0;
	
	virtual Real Area(int)=0;
	virtual Real K(int)=0;
	virtual Real Th(int)=0;
	virtual Real Kx(int)=0;
	virtual Real Ky(int)=0;
};

extern GeometryBase *Geometry;
Compare_t GeometryCompare;
KeyCompare_t GeometryKeyCompare;

class Weight {
	unsigned int index;
	Mc value;
public:
	void Store(unsigned int index0, Mc value0) {
		index=index0;
		value=value0;
	}
	unsigned int Index() const {return index;}
	Mc Value() const {return value;}
};

template<class T>
class Partition : public GeometryBase {
	DynVector<Weight> weight;
	int Nweight;
	Bin<T> *bin; // pointer to table of bins
	int n; // total number of bins, including reflected bins
public:
	Partition() {}
	char *Name();
	char *WeightFileName();
	int Create();
	void MakeBins();
	void ListBins(ostream &os);
	Mc ComputeBinAverage(Bin<T> *k, Bin<T> *p, Bin<T> *q);
	Mc Ckpq(T, T, T);
	Mc FindWeight(int k, int p, int q);
	void GenerateWeights();
	void ComputeTriads();
	void ListTriads();
	Nu BinAveragedLinearity(int i);
	int Number() {return n;}
	
	Real Area(int k) {return bin[k].Area();}
	Real K(int k) {return bin[k].K();}
	Real Th(int k) {return bin[k].Th();}
	Real Kx(int k) {return bin[k].Kx();}
	Real Ky(int k) {return bin[k].Ky();}
	
	int pq(int p, int q) {return n*p-p*(p+1)/2+q;} // Index to element p <= q
	unsigned int WeightIndex(int k, int p, int q) {
		return n*k*(n-k+2)/2+k*(k-1)*(k-2)/6+(2*n-k-p-1)*(p-k)/2+q-k-(p+1)
			-(2*n-1)*(k+1)+k*(k+1)+n; // Index to element k < p < q
	}
	
	Mc Central(int k, int p, int q)	{
		if(p == q) return 0.0;
		int sign=1;
		// Ensure that p < q
		sort2(p,q,sign);
		return sign*Ckpq(bin[k].cen,bin[p].cen,bin[q].cen);
	}
};

template<class T>
int Partition<T>::Create()
{
	MakeBins();
	psibuffer=new Var[reality ? 2*Nmode : Nmode];
	
	if(verbose > 2) {
		cout.precision(3);
		ListBins(cout);
		cout.precision(REAL_DIG);
	}
	
	ComputeTriads();
	if(verbose > 2) ListTriads();
	return Nmode;
}

template<class T>
void Partition<T>::GenerateWeights() {
	Mc binaverage;
	int k,p,q;
	unsigned int kpq,nkpq,lastkpq=0;
	double realtime,lasttime=time(NULL);
	double interval=15.0;
	
	Nweight=0;
	nkpq=WeightIndex(Nmode,Nmode,n);
	
	for(k=0; k < Nmode; k++) {	// Loop for k <= p <= q
		for(p=k+1; p < Nmode; p++) {
			for(q=p+1; q < n; q++) {
				if((coangular(&bin[k],&bin[p]) ||
					coangular(&bin[p],&bin[q]) ||
					coangular(&bin[q],&bin[k]))) binaverage=0.0;
				else binaverage=ComputeBinAverage(&bin[k],&bin[p],&bin[q]);
				
				if(binaverage) {
					if(verbose > 3) cout << k << "," << p << "," << q <<
										": " << binaverage << endl;
					kpq=WeightIndex(k,p,q);
					if(kpq > lastkpq || Nweight==0) {
						weight[Nweight++].Store(kpq,binaverage);
					} else msg(ERROR,"Index is not strictly increasing");
					lastkpq=kpq;
					
					realtime=time(NULL);
					if(realtime-lasttime > interval || verbose > 2) {
						lasttime=realtime;
						cout << 100*(kpq+1)/nkpq << 
							"% of weight factors generated (" << kpq+1 <<
								"/" <<	nkpq << ")." << endl;
					}
				}
			}
		}
	}
	if(Nweight >= 100) cout << endl;
}

template<class T>
inline Mc Partition<T>::FindWeight(int k, int p, int q) {
	// Exploit the full antisymmetry of the weight factors: reorder so that
	// $k\prime\le p\prime\le q\prime$.
	unsigned int kpq;
	int l,u,i,cmp,sign,conjflag=0;
	
	if(k==p || p==q || q==k) return 0.0;
	sign=1;
	
	// Ensure that k <= p <= q
	sort2(p,q,sign);
	sort2(k,p,sign);
	sort2(p,q,sign);

	if(p >= Nmode) {int K=k; k=p-Nmode; p=q-Nmode; q=K+Nmode; conjflag=1;}
	
	kpq=WeightIndex(k,p,q);
	l=0;
	u=Nweight;
	while(l < u) {
		i=(l+u)/2;
		cmp=kpq-weight[i].Index();
		if(cmp == 0) {
			Mc value=weight[i].Value();
			if(conjflag) value=conj(value);
			return sign*value;
		}
		if (cmp < 0) u=i;
		else if (cmp > 0) l=i+1;
	}
	return 0.0;	// no match found.
}

template<class T>
void Partition<T>::ComputeTriads() {
	Mc nkpq;
	Var **pqindex,**index;
	Real sym,denom;
	int k,p,q;
	ifstream fin;
	ofstream fout;
	char *filename;
	
	Ntriad=Npair=0;
	
	fin.open(filename=WeightFileName());
	if(fin) {
		fin.read((char *) &Nweight,sizeof(int));
		if(fin.eof()) fin.close();
	}
	if(fin) {
		(void) weight[Nweight-1];
		fin.read((char *) weight.Base(),Nweight*sizeof(Weight));
		if(!fin.good()) msg(ERROR,"Error reading from weight file %s",
							filename);
		fin.close();
	} else {
		errno=0;
		fout.open(filename);
		if(!fout) msg(ERROR,"Weight file %s could not be opened",filename);
		GenerateWeights();
		fout.write((char *) &Nweight,sizeof(int));
		fout.write((char *) weight.Base(),Nweight*sizeof(Weight));
		if(!fout.good()) msg(ERROR,"Error writing to weight file %s",
							 filename);
		fout.close();
	}
	
	pqindex=new Var*[pq(n,n)];
	pair=new Pair[pq(n,n)];
	for(p=0; p < n; p++) for(q=p; q < n; q++) pqindex[pq(p,q)]=NULL;
	
	for(k=0; k < Nmode; k++) {
		denom=twopi2*bin[k].Area();
		for(p=0; p < n; p++) {
			for(q=p; q < n; q++) {
				nkpq=FindWeight(k,p,q);
				nkpq *= Central(k,p,q);
				
				if(nkpq != 0.0)	{
					index=pqindex+pq(p,q);
					if(!*index) *index=pair[Npair++].Store(&psibuffer[p],
														   &psibuffer[q]);
					sym=(p==q) ? 0.5 : 1.0;
					nkpq *= sym/denom;
					
					triad[Ntriad++].Store(*index,nkpq);
				}
			}
		}
		triad[Ntriad++].Store(NULL,0.0);
	}
	
	triadBase=triad.Base();
	delete [] pqindex;
	weight.~DynVector();

	cout << endl << Npair << " WAVENUMBER PAIRS ALLOCATED." << endl <<
		Ntriad-Nmode << " WAVENUMBER TRIADS ALLOCATED." << endl;
}

template<class T>
void Partition<T>::ListTriads() {
	int j;
	cout << endl << Npair << " Pairs:" << endl;
	for(j=0; j < Npair; j++) { 
		cout << "(" << pair[j].p << "," << pair[j].q << ")" << endl;
	}
	cout << endl << Ntriad-Nmode << " Triads:" << endl;
	for(j=0; j < Ntriad; j++) {
		if(triad[j].pq)	cout << triad[j].Mkpq << endl;
	}
}

template<class T>
void Partition<T>::ListBins(ostream &os)
{
	os << "         " << Name() << " Bin Geometry:" << endl;
#if _CRAY	
	for(int i=0; i < n; i++) os << "[" << bin[i].min << "\t" << bin[i].cen <<
								 "\t" << bin[i].max << "]" << endl;
#else	
	for(int i=0; i < n; i++) os << bin[i];
#endif	
	cout << endl;
}

#endif
