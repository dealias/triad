#ifndef __Partition_h__
#define __Partition_h__ 1

#include <errno.h>

#include "Geometry.h"
#include "DynVector.h"

#define PARTITION(key,discrete) {\
new Entry<Partition<key,discrete>,GeometryBase>(#key,GeometryTable);}

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
extern int discrete;
extern int movie;
extern int truefield;

extern Var *pqbuffer;
extern Var **pqIndex;
extern int *qStart;
extern DynVector<Triad> triad;
extern TriadLimits *triadLimits;

template<class T, class D>
INLINE int InInterval(const D& m, const T& a, const T& b);

template<class T, class D>
class Bin : public Mode {
public:
	T min,max,cen;
	int nmode;
	DynVector<D> mode;
	
	Bin() {nmode=0;}
	T Delta() {return max-min;}
	Real Area();
	
	Real X() const {return cen.X();}
	Real Y() const {return cen.Y();}
	Real K2() const {return cen.K2();}
	Real K() const {return cen.K();}
	Real Th() const {return cen.Th();}
	
	int InBin(const D& m) {return InInterval(m,min,max);}
	void Count(const D& m) {if(InBin(m)) mode[nmode++]=m;}
	
	void MakeModes();
	
	void WriteModes(oxstream &s) {
		s << nmode;
		for(int i=0; i < nmode; i++) s << mode[i];
	}
	
	void ReadModes(ixstream &s) {
		s >> nmode;
		for(int i=0; i < nmode; i++) {
			s >> mode[i];
		}
	}
};

template<class T, class D>
INLINE ostream& operator << (ostream& os, const Bin<T,D>& y) {
	os << "[" << y.min << "\t" << y.cen << "\t" << y.max << "]";
	return os;
}

typedef unsigned short int Index_t;

class WeightIndex {
	Index_t k,p,q;
public:
	WeightIndex() {}
	WeightIndex(Index_t k0, Index_t p0, Index_t q0) : k(k0), p(p0), q(q0) {}
	inline operator double();
	friend inline int operator > (const WeightIndex& a, const WeightIndex& b);
	friend inline int operator < (const WeightIndex& a, const WeightIndex& b);
	friend inline int operator == (const WeightIndex& a, const WeightIndex& b);
	friend inline ixstream& operator >> (ixstream& s, WeightIndex& y);
	friend inline oxstream& operator << (oxstream& s, const WeightIndex& y);
	friend inline ostream& operator << (ostream& s, const WeightIndex& y);
};

extern WeightIndex WeightN;

inline int operator > (const WeightIndex& a, const WeightIndex& b)
{
	return a.k > b.k || (a.k == b.k && a.p > b.p) ||
		(a.k == b.k && a.p == b.p && a.q > b.q);
}
	
inline int operator < (const WeightIndex& a, const WeightIndex& b)
{
	return a.k < b.k || (a.k == b.k && a.p < b.p) ||
		(a.k == b.k && a.p == b.p && a.q < b.q);
}
	
inline int operator == (const WeightIndex& a, const WeightIndex& b)
{
	return (a.k == b.k && a.p == b.p && a.q == b.q);
}
	
inline WeightIndex::operator double()
{
	return ((double)q)/WeightN.k/WeightN.p/WeightN.q+
		((double)p)/WeightN.k/WeightN.p+((double)k)/WeightN.k;
}
	
extern const int maxbins;

inline ixstream& operator >> (ixstream& s, WeightIndex& y)
{
	unsigned int kp;
	s >> kp >> y.q;
	y.k=kp/maxbins;
	y.p=kp-maxbins*y.k;
    return s;
}

inline oxstream& operator << (oxstream& s, const WeightIndex& y)
{
	unsigned int kp=maxbins*y.k+y.p; 
	s << kp << y.q;
    return s;
}

inline ostream& operator << (ostream& s, const WeightIndex& y)
{
	s << y.k << " " << y.p << " " << y.q;
    return s;
}

class Weight {
	WeightIndex index;
	McWeight value;
public:
	void Store(WeightIndex index0, Mc value0) {
		index=index0;
		value=value0;
	}
	WeightIndex Index() const {return index;}
	Mc Value() const {return value;}
};

inline ixstream& operator >> (ixstream& s, Weight& y) {
	WeightIndex index;
	Mc value;
	s >> index >> value;
	y.Store(index,value);
	return s;
}

inline oxstream& operator << (oxstream& s, const Weight& y) {
	s << y.Index() << newl << y.Value();
	return s;
}

template<class T>
class Hash {
	int n;
	double factor,constant;
	int *table;
public:
	virtual inline int hash(T value) {
		return (int) (((double) value)*factor+constant);
	}
	virtual inline int hash_verify(T value) {
		int h=hash(value);
		if(h < 0 || h >= n) 
			msg(ERROR,"Hash for %d is outside the interval [0,%d]",value,n-1);
		return h;
	}
	
	Hash(int nhash, int nvalue, T value(int)) {
		if(nvalue < 1) return;
		n=nhash; 
		T first=value(0);
		factor=(n-1)/(((double) value(nvalue-1))-(double) (first));
		constant=0.5-((double) first)*factor;
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
inline WeightIndex HashWeightIndex(int j)
{
	return weightBase[j].Index();
}

template<class T, class D>
class Partition : public GeometryBase {
	DynVector<Weight> weight;
	Hash<WeightIndex> *hash;
	int Nweight,Nhash;
	Bin<T,D> *bin; // pointer to table of bins
public:
	char *Name();
	int Valid(char *s) {return strcmp(s,"SR")==0;}
	char *WeightFileName();
	void MakeBins();
	INLINE void List(ostream &os);
	Mc ComputeBinAverage(Bin<T,D> *k, Bin<T,D> *p, Bin<T,D> *q);
	Mc FindWeight(int k, int p, int q);

	INLINE void out_weights(oxstream& fout, int lastindex);
	INLINE void save_weights();
	INLINE int get_weights();
		
	INLINE void GenerateWeights();
	INLINE void Initialize();
	INLINE void ListTriads(ostream &os);
	
	Mode& ModeOf(int k) {return bin[k];}
	Real Area(int k) {return bin[k].Area();}
	
// Factor which converts |y|^2 to energy in various normalizations:
	Real Normalization(int);
	
	Nu Linear(int);
	Real Forcing(int);
	
	INLINE Mc Ckpq(T&, T&, T&);
	
	int pq(int p, int q) {return n*p-p*(p+1)/2+q;} // Index to element p <= q
};


template<class T, class D>
INLINE int coangular(Bin<T,D> *k, Bin<T,D> *p);

template<class T, class D>
INLINE void Partition<T,D>::GenerateWeights() {
	Mc binaverage;
	Index_t k,p,q;
	int i,first;
	int lastindex=0;
	WeightIndex kpq,previous,lastkpq=WeightIndex(0,0,0);
	double realtime,lasttime=time(NULL);
	double interval=15.0;
	oxstream fout;
	
	if(discrete && !Nweight) for (i=0; i < n; i++) bin[i].MakeModes();
	
	char *filename=WeightFileName();
	fout.open(filename);
	if(!fout) msg(ERROR,"Weight file %s could not be opened",filename);
	fout << lastindex;
	if(discrete) for(i=0; i < n; i++) bin[i].WriteModes(fout);
	
	if(Nweight) {
		out_weights(fout,lastindex);
		lastindex=Nweight;
		previous=weight[Nweight-1].Index();
		first=0;
	} else {
		previous=WeightIndex(0,0,0);
		first=1;
	}
	
	cout << newl << "GENERATING WEIGHT FACTORS." << endl;
	if(n > USHRT_MAX || n >= maxbins)
		msg(ERROR, "Internal limit on number of bins exceeded");	
	
	for(k=0; k < Nmode; k++) {	// Loop for k < p < q
		for(p=k+1; p < Nmode; p++) {
			for(q=p+1; q < n; q++) {
				kpq=WeightIndex(k,p,q);
				if(kpq > previous || first) {
					realtime=time(NULL);
					if(realtime-lasttime > interval || verbose > 2) {
						lasttime=realtime;
						cout << (int) (100.0*((double) kpq)+0.5) << 
							"% of weight factors generated: (" << kpq <<
							")/(" <<	WeightN << ")." << endl;
						out_weights(fout,lastindex);
						lastindex=Nweight;
					}
					if(!discrete && (coangular(&bin[k],&bin[p]) ||
									 coangular(&bin[p],&bin[q]) ||
									 coangular(&bin[q],&bin[k])))
						binaverage=0.0;
					else
						binaverage=ComputeBinAverage(&bin[k],&bin[p],&bin[q]);
				
					if(binaverage) {
						if(verbose > 3) cout << kpq << ": " << binaverage <<
											endl;
						if(kpq > lastkpq || first) {
							weight[Nweight++].Store(kpq,binaverage);
						} else msg(ERROR,"Index is not strictly increasing");
						lastkpq=kpq; first=0;
					}
				}
			}
		}
	}
	cout << endl;
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
	
	WeightIndex kpq=WeightIndex(k,p,q);
	int h=hash->hash(kpq);
	if(h < 0 || h >= Nhash) return 0.0; // no match found.
	int l=hash->Table(h);
	int u=hash->Table(h+1);
	
	while(l < u) {
		int i=(l+u)/2;
		Weight *w=weight+i;
		if(kpq == w->Index()) {
			Mc value=w->Value();
			if(conjflag) value=conj(value);
			return sign*value*Ckpq(bin[k0].cen,bin[p0].cen,bin[q0].cen);
		} else {
			if(kpq < w->Index()) u=i;
			else l=i+1;
		}
	}
	return 0.0;	// no match found.
}

template<class T, class D>
INLINE int Partition<T,D>::get_weights()
{
	int complete=0;
	int i,N=0;
	
	char *filename=WeightFileName();
	ixstream fin(filename);
	if(fin) {
		fin >> N;
		if(!fin.eof()) {
			if(N) complete=1;
			else testlock();
		
			cout << newl << "READING WEIGHT FACTORS FROM " << filename << "."
				 << endl;
		
			if(discrete) {
				for(i=0; i < n; i++) bin[i].ReadModes(fin);
				if(fin.bad())
					msg(ERROR,"Error reading from weight file %s",filename);
			}
		
			if(complete) {
				weight.Resize(N);
				for(i=0; i < N; i++) fin >> weight[i];
				if(fin.bad()) 
					msg(ERROR,"Error reading from weight file %s",filename);
			} else {
				Weight w;
				while(fin >> w, !fin.eof()) weight[N++]=w;
				if(N) cout << N << " WEIGHT FACTORS READ: ("
						   << weight[N-1].Index() << ")/(" << WeightN << ")."
						   << endl;
			}
		}
	}
	errno=0;

	Nweight=N;
	return complete;
}

template<class T, class D>
INLINE void Partition<T,D>::save_weights()
{
	int i;
	char *filename=WeightFileName();
	oxstream fout(filename);
	if(!fout) msg(ERROR,"Could not open weight file %s",filename);
	fout << Nweight << newl;
	if(discrete) for(i=0; i < n; i++) bin[i].WriteModes(fout);
	for(i=0; i < Nweight; i++) fout << weight[i] << newl;
	fout.close();
	if(!fout) msg(ERROR,"Cannot write to weight file %s",filename);
}	

template<class T, class D>
INLINE void Partition<T,D>::out_weights(oxstream& fout, int lastindex)
{
	lock();
	for(int i=lastindex; i < Nweight; i++) fout << weight[i] << newl;
	fout.flush();
	unlock();
}

template<class T, class D>
INLINE void Partition<T,D>::Initialize() {
	Mc nkpq;
	Real norm;
	int k,p,q;
	
	psibuffer=new Var[n];
	psibufferR=(reality ? psibuffer+Nmode : psibuffer);
	psibufferStop=psibuffer+n;
		
	WeightN=WeightIndex(Nmode,Nmode,n);
	if(!get_weights()) {
		GenerateWeights();
		save_weights();
	}
	
	if(discrete && movie && truefield) {
		int i;
		Nevolved=(Ny-1)/2*Nx+(Nx-1)/2;
		int Ndiscrete=Nx*Ny-1;
		ModeBin=new int[Ndiscrete];
		const int Unassigned=-1;
		for(i=0; i < Ndiscrete; i++) ModeBin[i]=Unassigned;
		
		for (i=0; i < n; i++) {
			Bin<T,D> *p=&bin[i];
			int nmode=p->nmode;
			D *mode=p->mode.Base();
			for(int m=0; m < nmode; m++) {
				int index=mode[m].ModeIndex();
				if(index >= Ndiscrete) 
					msg(ERROR, "Number of modes is greater than Nx*Ny-1");
				if(ModeBin[index] != Unassigned) 
					msg(ERROR,
						"Discrete mode (%d) is already assigned to bin %d",
						m,ModeBin[index]);
				ModeBin[index]=i;
			}
		}
	}
	
	weightBase=weight.Base();

	Nhash=2*Nweight;
	hash=new Hash<WeightIndex>(Nhash,Nweight,HashWeightIndex);
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
	if(verbose > 2) ListTriads(cout);
}

template<class T, class D>
INLINE void Partition<T,D>::ListTriads(ostream &os) {
	int j;
	os << newl << Ntriad << " Triads:" << endl;
	for(j=0; j < Ntriad; j++) {
		if(triad[j].pq)	os << triad[j].Mkpq << newl;
	}
	os << flush;
}

template<class T, class D>
INLINE void Partition<T,D>::List(ostream &os)
{
	os << "         " << Name() << " Bin Geometry:" << endl;
	for(int i=0; i < n; i++) {
#if _CFRONT	// Work around Cfront bug
	os << "[" << bin[i].min << "\t" << bin[i].cen << "\t" << bin[i].max
	   << "]" << newl;
#else	
	os << bin[i] << newl;
#endif	
	}
	os << endl;
}

#endif
