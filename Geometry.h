#ifndef __Geometry_h__
#define __Geometry_h__ 1

#include <errno.h>

#include "kernel.h"
#include "utils.h"
#include "Pair.h"
#include "Bin.h"

#define GEOMETRY(key) {new Entry<Partition<key>,GeometryBase>\
						   (#key,GeometryTable);}

class GeometryBase {
protected:	
	int Nmode; // number of unreflected (explicitly evolved) bins
	int n; // total number of bins, including reflected bins
public:	
	virtual char *Name()=0;
	virtual int Create()=0;
	int TotalNumber() {return n;}
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

inline istream& operator >> (istream& s, Weight& y) {
	unsigned int index;
	Mc value;
	s >> index >> value;
	y.Store(index,value);
	return s;
}

inline ostream& operator << (ostream& s, const Weight& y) {
	s << y.Index() << endl << y.Value();
	return s;
}

template<class T>
class Partition : public GeometryBase {
	DynVector<Weight> weight;
	int Nweight;
	Bin<T> *bin; // pointer to table of bins
public:
	Partition() {}
	char *Name();
	char *WeightFileName(char *suffix);
	int Create();
	void MakeBins();
	void ListBins(ostream &os);
	Mc ComputeBinAverage(Bin<T> *k, Bin<T> *p, Bin<T> *q);
	inline Mc Ckpq(T, T, T);
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
};

template<class T>
int Partition<T>::Create()
{
	MakeBins();
	psibuffer=new Var[n];
	psibufferR=(reality ? psibuffer+Nmode : psibuffer);
	psibufferStop=psibuffer+n;

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
	
	for(k=0; k < Nmode; k++) {	// Loop for k < p < q
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
	if(k==p || p==q || q==k) return 0.0;
	
	int sign=1, conjflag=0, k0=k, dummy;
	
	// Exploit the full antisymmetry of the weight factors: reorder so that
	// k <= p <= q.

	sort2(p,q,dummy);
	int p0=p, q0=q;
	
	sort2(k,p,sign);
	sort2(p,q,sign);

	if(p >= Nmode) {int K=k; k=p-Nmode; p=q-Nmode; q=K+Nmode; conjflag=1;}
	
	unsigned int kpq=WeightIndex(k,p,q);
	int l=0;
	int u=Nweight;
	while(l < u) {
		int i=(l+u)/2;
		int cmp=kpq-weight[i].Index();
		if(cmp < 0) u=i;
		else {
			if(cmp > 0) l=i+1;
			else {
				Mc value=weight[i].Value();
				if(conjflag) value=conj(value);
				return sign*value*Ckpq(bin[k0].cen,bin[p0].cen,bin[q0].cen);
			}
		}
	}
	return 0.0;	// no match found.
}

template<class T>
void Partition<T>::ComputeTriads() {
	Mc nkpq;
	Real norm;
	int k,p,q;
	ifstream fin;
	ofstream fout;
	char *filename;
	int i,formatted=0;
	
	Ntriad=0;
	
	filename=WeightFileName("");
	fin.open(filename);
	if(!fin) {
		filename=WeightFileName("f");
		fin.open(filename);
		formatted=1;
	}
	
	if(fin) {
		if(formatted) fin >> Nweight;
		else fin.read((char *) &Nweight,sizeof(int));
		if(fin.eof() || !Nweight) fin.close();
	}
	
	if(fin) {
		(void) weight[Nweight-1];
		if(formatted) for(i=0; i < Nweight; i++) fin >> weight[i];
		else fin.read((char *) weight.Base(),Nweight*sizeof(Weight));
		fin.close();
		if(!fin.good()) msg(ERROR,"Error reading from weight file %s",
							filename);
		if(!formatted && output) {
			filename=WeightFileName("f");
			fout.open(filename);
			fout.precision(REAL_DIG);
			if(!fout) msg(ERROR,"Weight file %s could not be opened",filename);
			fout << Nweight << endl;
			for(i=0; i < Nweight; i++) fout << weight[i] << endl;
			fout.close();
			if(!fout.good())
				msg(ERROR,"Error writing to weight file %s",filename);
		}
	}
	
	if(!fin || formatted) {
		errno=0;
		filename=WeightFileName("");
		fout.open(filename);
		if(!fout) msg(ERROR,"Weight file %s could not be opened",filename);
		if(!fin) GenerateWeights();
		fout.write((char *) &Nweight,sizeof(int));
		fout.write((char *) weight.Base(),Nweight*sizeof(Weight));
		fout.close();
		if(!fout.good()) msg(ERROR,"Error writing to weight file %s",filename);
	}
	
	int npq=reality ? Nmode*(3*Nmode+1)/2 : pq(n,n);
	pqbuffer=new Var[npq];
	pqIndex=new Var*[n];
	triadLimits=new TriadLimits[Nmode];
	int *ntriad=new int[Nmode];
	
	int NmodeR=psibufferR-psibuffer;
	Var *pq=pqbuffer;
	for(p=0; p < n; p++) {
			int m=max(NmodeR,p);
			pqIndex[p]=pq-m;
			if(n > m) pq += n-m;
	}
	
	triad.Resize(Nmode*n);
	for(k=0; k < Nmode; k++) {
		norm=1.0/(twopi2*Area(k));
		pq=pqbuffer;
		for(p=0; p < n; p++) {
			for(q=max(NmodeR,p); q < n; q++, pq++) {
				nkpq=FindWeight(k,p,q);
				if(nkpq != 0.0)	{
					if(p==q) nkpq *= 0.5;
					triad[Ntriad++].Store(pq,nkpq*norm);
				}
			}
		}
		ntriad[k]=Ntriad;
	}
	triad.Resize(Ntriad);
	
	triadLimits[0].start=triad.Base();
	for(k=0; k < Nmode-1; k++) {
		triadLimits[k+1].start=triadLimits[k].stop=triad.Base()+ntriad[k];
	}
	triadLimits[Nmode-1].stop=triad.Base()+Ntriad;

	delete [] ntriad;
	weight.~DynVector();

	cout << endl << Ntriad-Nmode << " WAVENUMBER TRIADS ALLOCATED." << endl;
}

template<class T>
void Partition<T>::ListTriads() {
	int j;
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
