#ifndef __Basis_h__
#define __Basis_h__ 1

#include "Geometry.h"
#include "DynVector.h"

#define BASIS(key) {new Entry<Basis<key>,GeometryBase> (#key,GeometryTable);}

class Chain {
public:
	int k,p; // mode indices
	int stop; // column index where chain stops
	Var *psiq;
	void Store(int k0, int p0, int stop0, Var *q0) {
		k=k0; p=p0; stop=stop0; psiq=q0;
	}
};

extern int Nchainp,Nchainn,Ntriad;
extern DynVector<Chain> chainp,chainn;
extern Chain *chainpBase,*chainnBase;
extern Real *kinv2;

template<class T>
class Basis : public GeometryBase {
	T *mode; // pointer to table of modes
	T low; // lower limits of grid
	T high; // upper limits of grid
public:
	char *Name();
	int ValidApproximation(char *s) {
		return (strcmp(s,"NONE")==0 || strcmp(s,"PS")==0);
	}

	void MakeBins();
	void List(ostream &);
	void ListTriads();
	void ListChains(Chain *chain, int Nchain, char *type);
	void ComputeTriads();
	
	int InGrid(T &);
	Real Area(int) {return 1.0;}
	
	Real K(int k) {return mode[k].K();}
	Real K2(int k) {return mode[k].K2();}
	Real Th(int k) {return mode[k].Th();}
	Real Kx(int k) {return mode[k].Kx();}
	Real Ky(int k) {return mode[k].Ky();}
	
// Factor which converts |y|^2 to energy in various normalizations:
	Real Normalization(int);
	
	Nu Linearity(int);
	Mc Mkpq(T& k, T& p, T& q);
};

template<class T>
inline void Basis<T>::ListChains(Chain *chain, int Nchain, char *type) {
	cout << endl << Nchain << " " << type << " Chains:" << endl;
	Chain *c,*chainStop=chain+Nchain;
	for(c=chain; c < chainStop; c++) {
		cout << "k=" << mode[c->k] << " p=" << mode[c->p] << " stop=" <<
			c->stop << " q=" << mode[c->psiq-psibuffer] << endl;
	}
}

template<class T>
void Basis<T>::ListTriads() {
	ListChains(chainp.Base(),Nchainp,"Positive");
	ListChains(chainn.Base(),Nchainn,"Negative");
}

template<class T>
void Basis<T>::List(ostream &os)
{
	os << "         " << Name() << " Mode Geometry:" << endl;
	for(int i=0; i < n; i++) os << mode[i] << endl;
	cout << endl;
}

inline void StoreChain(int k, int p, int stop, Var *q, int sign)
{
	if(sign >= 0) chainp[Nchainp++].Store(k,p,stop,q);
	else chainn[Nchainn++].Store(k,p,stop,q);
}


template<class T>
void Basis<T>::ComputeTriads()
{
	int k,p,q,sign=0,newchain;
	int ck,cp,cq;
	int lastp=-1;
	T mq;
	
	kinv2=new Real[Nmode];
	for(k=0; k < Nmode; k++) kinv2[k]=1.0/mode[k].K2();

	if(pseudospectral) return;
	
	chainp.Resize(8*Nmode);
	chainn.Resize(4*Nmode);
	
	Ntriad=Nchainp=Nchainn=0;
	ck=cp=cq=q=0;
	
	for(k=0; k < Nmode; k++) {
		newchain=1;
		for(p=0; p < n; p++) {
			mq = -(mode[k]+mode[p]);
			if(mode[p].Row() != mode[cp].Row() || mq.Row() != mode[cq].Row())
			   newchain=1;
			if(InGrid(mq) && !(newchain && Mkpq(mode[k],mode[p],mq) == 0.0)) {
				if(!newchain) {
				    if(sign == 0) {
						if(q < n-1 && mode[q+1] == mq) sign=1;
						else if(q > 0 && mode[q-1] == mq) sign=-1;
					}
					q += sign;
				}
				if(newchain || mode[q] != mq) {
					for(q=0; q < n && mode[q] != mq; q++);
					if(q == n) msg(ERROR, "Invalid beat mode computed");
					if(lastp >= 0) StoreChain(ck,cp,mode[lastp].Column()+1,
											  psibuffer+cq,sign);
					ck=k; cp=p; cq=q;
					sign=0; newchain=0;
				}
				Ntriad++;
				lastp=p;
			}
			else newchain=1;
		}
	}
	if(lastp >= 0) StoreChain(ck,cp,mode[lastp].Column()+1,psibuffer+cq,sign);

	chainp.Resize(Nchainp);
	chainpBase=chainp.Base();
	
	chainn.Resize(Nchainn);
	chainnBase=chainn.Base();

	cout << Nchainp+Nchainn << " WAVENUMBER CHAINS ALLOCATED." << endl;
	cout << Ntriad << " WAVENUMBER TRIADS ALLOCATED." << endl;
}

Nu LinearityAt(int i);

template <class T>
Nu Basis<T>::Linearity(int i)
{
	return LinearityAt(i);
}

#endif
