#ifndef __Basis_h__
#define __Basis_h__ 1

#include "Geometry.h"
#include "DynVector.h"

#define BASIS(key) {new Entry<Basis<key>,GeometryBase> (#key,GeometryTable);}

class Chain{
public:
	int k;
	Var *pstart,*pstop;
	Var *qstart;
	Mc *Mkpq;
	void Store(int k0, Var *pstart0, Var *pstop0, Var *qstart0, Mc *Mkpq0) {
		k=k0; pstart=pstart0; pstop=pstop0; qstart=qstart0; Mkpq=Mkpq0;
		if(pstop < pstart)
			msg(ERROR,"Negative length chain [%x,%x] encountered",
				pstart,pstop);
	}
};

extern int Nchainp,Nchainn;
extern DynVector<Chain> chainp,chainn;
extern Chain *chainpBase,*chainnBase;

template<class T>
class Basis : public GeometryBase {
	T *mode; // pointer to table of modes
	T period; // periodicity
public:
	char *Name();
	char *Approximation() {return "NONE";}
	void MakeBins();
	void List(ostream &);
	void ListTriads();
	void ComputeTriads();
	
	void SetPeriod(T period0) {period=period0;}
	int InGrid(T &);
	
	Nu Linearity(int);
	Mc Ckpq(T& k, T& p, T& q);

	Real Area(int) {return 1.0;}
	Real K(int k) {return mode[k].K();}
	Real Th(int k) {return mode[k].Th();}
	Real Kx(int k) {return mode[k].Kx();}
	Real Ky(int k) {return mode[k].Ky();}
};

inline void ListChains(Chain *chain, int Nchain, char *type) {
	cout << endl << Nchain << " " << type << " Chains:" << endl;
	Chain *c,*chainStop=chain+Nchain;
	for(c=chain; c < chainStop; c++) {
		cout << "k=" << c->k << " p1=" << c->pstart-psibuffer << " p2=" <<
			c->pstop-psibuffer << " q=" << c->qstart-psibuffer << " Mkpq=" <<
			*(c->Mkpq) << endl;
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

inline void StoreChain(int k, Var *pstart, Var *pstop, Var *qstart, Mc *mstart,
					   int sign)
{
	if(sign >= 0) chainp[Nchainp++].Store(k,pstart,pstop,qstart,mstart);
	else chainn[Nchainn++].Store(k,pstart,pstop,qstart,mstart);
}


template<class T>
void Basis<T>::ComputeTriads()
{
	Mc *m,*Mkpq,*mstart;
	Var *pstart=psibuffer,*pstop=psibuffer,*qstart=psibuffer;
	int k,p,q,sign,newchain;
	int kstart=0, lastp=-1;
	T mq;
	mstart=m=Mkpq=new Mc[Nmode*n];
	chainp.Resize(8*Nmode);
	chainn.Resize(4*Nmode);
	Nchainp=Nchainn=0;
	q=0;

	for(k=0; k < Nmode; k++) {
		newchain=1;
		for(p=0; p < n; p++) {
			mq = -(mode[k]+mode[p]);
			if(InGrid(mq)) {
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
					if(lastp >= 0) {
						pstop=psibuffer+lastp+1;
						StoreChain(kstart,pstart,pstop,qstart,mstart,sign);
					}
					kstart=k; pstart=psibuffer+p; qstart=psibuffer+q; mstart=m;
					sign=0; newchain=0;
				}
				(*m++)=Ckpq(mode[k],mode[p],mode[q]);
				lastp=p;
			}
			else newchain=1;
		}
	}
	pstop=psibuffer+lastp+1;
	StoreChain(kstart,pstart,pstop,qstart,mstart,sign);

	chainp.Resize(Nchainp);
	chainpBase=chainp.Base();
	
	chainn.Resize(Nchainn);
	chainnBase=chainn.Base();

	cout << Nchainp+Nchainn << " WAVENUMBER CHAINS ALLOCATED." << endl;
	cout << m-Mkpq << " WAVENUMBER TRIADS ALLOCATED." << endl;
}

Nu LinearityAt(int i);

template <class T>
Nu Basis<T>::Linearity(int i)
{
	return LinearityAt(i);
}

#endif
