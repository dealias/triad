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

typedef unsigned int Index_t;

class Weight {
	Index_t index;
	McWeight value;
public:
	void Store(Index_t index0, Mc value0) {
		index=index0;
		value=value0;
	}
	Index_t Index() const {return index;}
	Mc Value() const {return value;}
};

inline istream& operator >> (istream& s, Weight& y) {
	Index_t index;
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
class Hash {
	int n;
	T first, last;
	double factor,constant;
	int *table;
public:
	virtual inline int hash(T value) {
		return (int) (value*factor+constant);
	}
	virtual inline int hash_verify(T value) {
		int h=hash(value);
		if(h < 0 || h >= n) 
			msg(ERROR,"Hash for %d is outside the interval [0,%d]",value,n-1);
		return h;
	}
	
	Hash(int nhash, int nvalue, T value(int)) {
		n=nhash; first=value(0); last=value(nvalue-1);
		factor=((double) n-1)/(last-first); constant=0.5-first*factor;
		table=new int[n+1];
		int j=0;
		for(int i=0; i < n; i++) {
			while (hash_verify(value(j)) < i) j++;
			if(j >= nvalue)
				msg(ERROR,"Hash table entry %d is out of bounds",j);
			table[i]=j;
		}
		table[n]=nvalue;
	}
	
	inline int Table(int hash) {return table[hash];}
};
	
extern Weight *weightBase;
inline Index_t HashWeightIndex(int j)
{
	return weightBase[j].Index();
}

template<class T>
class Partition : public GeometryBase {
	DynVector<Weight> weight;
	Hash<Index_t> *hash;
	int Nweight,Nhash;
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
	Index_t WeightIndex(int k, int p, int q) {
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
	int k,p,q,first;
	int lastindex=0;
	Index_t kpq,nkpq,previous,lastkpq=0;
	double realtime,lasttime=time(NULL);
	double interval=15.0;
	ofstream fout;
	
	if(Nweight) {
		previous=weight[Nweight-1].Index();
		first=0;
	} else {
		previous=0;
		first=1;
	}
	
	char *filename=WeightFileName("");
	fout.open(filename);
	if(!fout) msg(ERROR,"Weight file %s could not be opened",filename);
	fout.write((char *) &lastindex,sizeof(int));
	
	cout << endl << "GENERATING WEIGHT FACTORS." << endl;
	nkpq=WeightIndex(Nmode,Nmode,n);
	
	for(k=0; k < Nmode; k++) {	// Loop for k < p < q
		for(p=k+1; p < Nmode; p++) {
			for(q=p+1; q < n; q++) {
				kpq=WeightIndex(k,p,q);
				if(kpq > previous || first) {
					if((coangular(&bin[k],&bin[p]) ||
						coangular(&bin[p],&bin[q]) ||
						coangular(&bin[q],&bin[k]))) binaverage=0.0;
					else binaverage=ComputeBinAverage(&bin[k],&bin[p],&bin[q]);
				
					if(binaverage) {
						if(verbose > 3) cout << k << "," << p << "," << q <<
											": " << binaverage << endl;
						if(kpq > lastkpq || first) {
							weight[Nweight++].Store(kpq,binaverage);
						} else msg(ERROR,"Index is not strictly increasing");
						lastkpq=kpq; first=0;
					}
					realtime=time(NULL);
					if(realtime-lasttime > interval || verbose > 2) {
						lasttime=realtime;
						cout << 100*(kpq+1)/nkpq << 
							"% of weight factors generated (" << kpq+1 <<
								"/" <<	nkpq << ")." << endl;
						lock();
						fout.write((char *) (weight.Base()+lastindex),
								   Nweight*sizeof(Weight));
						fout.flush();
						unlock();
						lastindex=Nweight;
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
	
	int dummy, sign=1, conjflag=0, k0=k;
	
	// Exploit the full antisymmetry of the weight factors: reorder so that
	// k <= p <= q.

	sort2(p,q,dummy);
	int p0=p, q0=q;
	
	sort2(k,p,sign);
	sort2(p,q,sign);

	if(p >= Nmode) {int K=k; k=p-Nmode; p=q-Nmode; q=K+Nmode; conjflag=1;}
	
	Index_t kpq=WeightIndex(k,p,q);
	int h=hash->hash(kpq);
	if(h < 0 || h >= Nhash) return 0.0; // no match found.
	int l=hash->Table(h);
	int u=hash->Table(h+1);
	
	while(l < u) {
		int i=(l+u)/2;
		int cmp=kpq-weight[i].Index();
		if(cmp == 0) {
			Mc value=weight[i].Value();
			if(conjflag) value=conj(value);
			return sign*value*Ckpq(bin[k0].cen,bin[p0].cen,bin[q0].cen);
		}
		if(cmp < 0) u=i;
		else if(cmp > 0) l=i+1;
	}
	return 0.0;	// no match found.
}

int get_weights(DynVector<Weight>& weight, int *Nweight, char *filename,
				char *filenamef);
void save_weights(DynVector<Weight>& w, int n, char *filename);
void save_formatted_weights(DynVector<Weight>& w, int n, char *filename);

template<class T>
void Partition<T>::ComputeTriads() {
	Mc nkpq;
	Real norm;
	int k,p,q;
	char *filename=WeightFileName(""),*filenamef=WeightFileName("f");
	
	if(!get_weights(weight,&Nweight,filename,filenamef)) {
		GenerateWeights();
		save_weights(weight,Nweight,filename);
	}
	if(output) save_formatted_weights(weight,Nweight,filenamef);

	Ntriad=0;
	
	int npq=reality ? Nmode*(3*Nmode+1)/2 : pq(n,n);
	pqbuffer=new Var[npq];
	pqIndex=new Var*[n];
	qStart=new int[n];
	triadLimits=new TriadLimits[Nmode];
	int *ntriad=new int[Nmode];
	
	int NmodeR=psibufferR-psibuffer;
	Var *pq=pqbuffer;
	for(p=0; p < n; p++) {
			int m=max(NmodeR,p);
			pqIndex[p]=pq-m;
			if(n > m) pq += n-m;
	}
	for(p=0; p < n; p++) qStart[p]=max(p,NmodeR);
	
	weightBase=weight.Base();

	Nhash=2*Nweight;
	hash=new Hash<Index_t>(Nhash,Nweight,HashWeightIndex);
	cout << "HASH TABLE CONSTRUCTED." << endl;

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
	
	delete hash;
	
	triadLimits[0].start=triad.Base();
	for(k=0; k < Nmode-1; k++) {
		triadLimits[k+1].start=triadLimits[k].stop=triad.Base()+ntriad[k];
	}
	triadLimits[Nmode-1].stop=triad.Base()+Ntriad;

	delete [] ntriad;
	weight.~DynVector();

	cout << Ntriad-Nmode << " WAVENUMBER TRIADS ALLOCATED." << endl;
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
