#ifndef __Param_h__
#define __Param_h__ 1

#include <iostream.h>
#include <stdlib.h>
#include <string.h>

#include "kernel.h"

extern int NParam;

inline void VocabularyBase::ParamAdd(ParamBase *p)
{
	ParamList[NParam++]=p;
}

#define VOCAB(var,min,max) Vocab(&var, #var, min, max, 1, 1)
#define VOCAB_NODUMP(var,min,max) Vocab(&var, #var, min, max, 1, 0)

#define VOCAB_ARRAY(var) \
Vocab(var, #var, *var-*var, *var-*var, (int) (sizeof(var)/sizeof(*var)),1)
	
template<class T>
class Param : public ParamBase {
	char *name;
	int nvar;
	T *var;
	T min;
	T max;
	int dump;
public:
	Param(T *address, char *s, int n, T min0, T max0, int dump0) {
		ParamBase *param;
		int match_type;

		param=Vocabulary->Locate(s,&match_type);
		if(match_type==2) {
			cout << beep << "Duplicate vocabulary entry: " << s << endl;
			return;
		}
		
		name=s; nvar=n; var=address; min=min0; max=max0; dump=dump0;
		Vocabulary->ParamAdd(this);
		Vocabulary->Sort();
	}

	void Set(T x) {int i; for(i=0; i<nvar; i++) var[i]=x;}
	void Set(T *x) {int i; for(i=0; i<nvar; i++) var[i]=x[i];}
	
	void SetStr(char *);
	
	void Retrieve(T *p) {*p=*var;}
	char *Name() {return name;}
	
	void Display(ostream& os) {
		os << name << " = ";
		for(int i=0; i < nvar-1; i++) os << var[i] << ", ";
		os << var[nvar-1] << endl;
	}
	
	void Output(ostream& os) {
		if(!dump) return; // Don't dump control parameters
		os << name << "=";
		for(int i=0; i < nvar-1; i++) os << var[i] << ",";
		os << var[nvar-1] << endl;
	}
	
	void GraphicsOutput(ostream& os) {
		if(!dump) return; // Don't dump control parameters
		if(nvar==1) os << "define " << name << " \"" << var[0] << "\"";
		else {
			os << "set " << name << "={";
			for(int i=0; i < nvar; i++) os << var[i] << " ";
			os << "}";
		}
		os << endl;
	}
	
	inline int InRange(T);
	
	void get_values(char *optarg, T (*rtn)(const char *))
	{
		int i;
		char *ptr,*optarg0;
		T value;

		i=0;
		do {
			if(!(optarg0=strchr(optarg,')'))) optarg0=optarg;
			if((ptr=strchr(optarg0,','))) *ptr=0;
			else ptr=optarg;
			if(i >= nvar) {
				cout << beep << "Warning: excess initializer " << optarg <<
					" to " << name << " discarded." << endl;
			}
			else {
				value=(*rtn)(optarg);

				if(InRange(value)) var[i++]=value;
				else {
					cout << beep << "Value " << value << " for " <<
						name << " is invalid: Limits are " <<
						min << " to " << max << "." << endl;
				}
			}	
		} while ((ptr==optarg) ? 0 : (optarg=ptr+1));
	}
	
};

template<class T>
inline void Vocab(T *var, char *s, T min, T max, int n, int dump)
{
	new Param<T>(var,s,n,min,max,dump);
}

char *atos(const char *s);

inline void Param<char *>::SetStr(char *s) {get_values(s,&atos);}
inline void Param<double>::SetStr(char *s) {get_values(s,&atof);}
inline void Param<int>::SetStr(char *s) {get_values(s,&atoi);}
inline void Param<Complex>::SetStr(char *s) {get_values(s,&atoc);}

inline int Param<double>::InRange(double value)
{
	return (min == max) || (value >= min && value <= max);
}

inline int Param<int>::InRange(int value)
{
	return (min == max) || (value >= min && value <= max);
}

inline int Param<char *>::InRange(char *) {return 1;}
inline int Param<Complex>::InRange(Complex) {return 1;}

#endif
