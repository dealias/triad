#ifndef __Partition_h__
#define __Partition_h__ 1

#include <errno.h>

#include "Geometry.h"
#include "DynVector.h"
#include "Cartesian.h"

#define PARTITION(key,discrete) {\
new Entry<Partition<key,discrete>,GeometryBase>(#key,GeometryTable);}

typedef unsigned int Index_t;

class Triad {
public:
	Var *pq;
	Mc Mkpq;
	void Store(Var *pqindex, Mc value) {pq=pqindex; Mkpq=value;}
};

struct TriadLimits {
	Triad *start;
	Triad *stop;
};

extern int Ntriad;

extern Var *pqbuffer;
extern Var **pqIndex;
extern int *qStart;
extern DynVector<Triad> triad;
extern TriadLimits *triadLimits;

template<class D>
class Discrete
{
public:
	D value;
	Real weight;
	void Store(const D& value0, const Real weight0) {
		value=value0; weight=weight0;
	}
};

template<class D>
ostream& operator << (ostream& os, const Discrete<D>& y) {
	os << y.value << ": " << y.weight;
	return os;
}

template<class T, class D>
class Bin {
public:
	T min,max,cen;
	int nmode;
	Real area;
	DynVector<Discrete<D> > mode;
	
	Bin() {nmode=0; area=0.0;}
	T Delta() {return (max-min);}
	Real Area();
	
	Real K() {return cen.K();}
	Real K2() {return cen.K2();}
	Real Th() {return cen.Th();}
	Real X() {return cen.X();}
	Real Y() {return cen.Y();}
	
	Real InBin(const D& m) {return InInterval(m,min,max);}
	void Count(const D& m) {
		const Real weight=InBin(m); 
		if(weight) {
			area += weight;
			mode[nmode++].Store(m,weight);
		}
	}
	void MakeModes();
};

template<class T, class D>
ostream& operator << (ostream& os, const Bin<T,D>& y) {
	os << "[" << y.min << "\t" << y.cen << "\t" << y.max << "]";
	if(discrete) os << ": " << y.area;
	os << endl;
	return os;
}

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
		if(nvalue < 1) return;
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

template<class T, class D>
class Partition : public GeometryBase {
	DynVector<Weight> weight;
	Hash<Index_t> *hash;
	int Nweight,Nhash;
	Index_t nmax;
	Bin<T,D> *bin; // pointer to table of bins
public:
	char *Name();
	int ValidApproximation(char *s) {return strcmp(s,"SR")==0;}
	char *WeightFileName(char *suffix);
	void MakeBins();
	void List(ostream &os);
	Mc ComputeBinAverage(Bin<T,D> *k, Bin<T,D> *p, Bin<T,D> *q);
	Mc FindWeight(int k, int p, int q);

	void GenerateWeights();
	void Initialize();
	void ListTriads();
	
	Real Area(int k) {return bin[k].Area();}
	Real K(int k) {return bin[k].K();}
	Real K2(int k) {return bin[k].K2();}
	Real Th(int k) {return bin[k].Th();}
	Real X(int k) {return bin[k].X();}
	Real Y(int k) {return bin[k].Y();}
	
// Factor which converts |y|^2 to energy in various normalizations:
	Real Normalization(int);
	
	Nu Linearity(int i);
	inline Mc Ckpq(T&, T&, T&);
	
	int pq(int p, int q) {return n*p-p*(p+1)/2+q;} // Index to element p <= q
	Index_t WeightIndex(int k, int p, int q) {
		return n*k*(n-k+2)/2+k*(k-1)*(k-2)/6+(2*n-k-p-1)*(p-k)/2+q-k-(p+1)
			-(2*n-1)*(k+1)+k*(k+1)+n; // Index to element k < p < q
	}
};


int out_weights(ofstream& fout, Weight* w, int lastindex, int n);
	
template<class T, class D>
void Partition<T,D>::GenerateWeights() {
	Mc binaverage;
	int k,p,q,first;
	int lastindex=0;
	Index_t kpq,previous,lastkpq=0;
	double realtime,lasttime=time(NULL);
	double interval=15.0;
	ofstream fout;
	
	char *filename=WeightFileName("");
	fout.open(filename);
	if(!fout) msg(ERROR,"Weight file %s could not be opened",filename);
	fout.write((char *) &lastindex,sizeof(int));
	
	if(Nweight) {
		lastindex=out_weights(fout,weight.Base(),lastindex,Nweight);
		previous=weight[Nweight-1].Index();
		first=0;
	} else {
		previous=0;
		first=1;
	}
	
	cout << endl << "GENERATING WEIGHT FACTORS." << endl;
	
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
						cout << 100*(kpq+1)/nmax << 
							"% of weight factors generated (" << kpq+1 <<
							"/" <<	nmax << ")." << endl;
						lastindex=out_weights(fout,weight.Base(),lastindex,
											  Nweight);
					}
				}
			}
		}
	}
	if(Nweight >= 100) cout << endl;
}




template<class T, class D>
inline Mc Partition<T,D>::FindWeight(int k, int p, int q) {
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

int get_weights(DynVector<Weight>& weight, int nmax, int *Nweight,
				char *filename,	char *filenamef);
void save_weights(DynVector<Weight>& w, int n, char *filename);
void save_formatted_weights(DynVector<Weight>& w, int n, char *filename);

template<class T, class D>
void Partition<T,D>::Initialize() {
	Mc nkpq;
	Real norm;
	int k,p,q;
	
	psibuffer=new Var[n];
	psibufferR=(reality ? psibuffer+Nmode : psibuffer);
	psibufferStop=psibuffer+n;
		
	char *filename=WeightFileName(""), *filenamef=WeightFileName("f");
	
	nmax=WeightIndex(Nmode,Nmode,n);
	if(!get_weights(weight,nmax,&Nweight,filename,filenamef)) {
		GenerateWeights();
		save_weights(weight,Nweight,filename);
	}
	if(output) save_formatted_weights(weight,Nweight,filenamef);
	
	weightBase=weight.Base();

	Nhash=2*Nweight;
	hash=new Hash<Index_t>(Nhash,Nweight,HashWeightIndex);
	cout << "HASH TABLE CONSTRUCTED." << endl;
	
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
	
	triad.Resize(Nmode*n);
	for(k=0; k < Nmode; k++) {
		norm=Area(k);
		if(norm != 0.0) norm=1.0/norm;
		pq=pqbuffer;
		for(p=0; p < n; p++) {
			for(q=max(NmodeR,p); q < n; q++, pq++) {
				nkpq=FindWeight(k,p,q);
				if(nkpq != 0.0)	{
					if(p==q) nkpq *= 0.5;
					if(norm == 0.0) msg(ERROR,"Invalid weight factor");
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

	cout << Ntriad << " WAVENUMBER TRIADS ALLOCATED." << endl;
	if(verbose > 2) ListTriads();
}

template<class T, class D>
void Partition<T,D>::ListTriads() {
	int j;
	cout << endl << Ntriad << " Triads:" << endl;
	for(j=0; j < Ntriad; j++) {
		if(triad[j].pq)	cout << triad[j].Mkpq << endl;
	}
}

template<class T, class D>
void Partition<T,D>::List(ostream &os)
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
