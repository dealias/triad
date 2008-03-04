#ifndef __Param_h__
#define __Param_h__ 1

#include <cstdlib>

#include "kernel.h"
#include "Integrator.h"

extern unsigned int NParam;
extern int param_warn;

class ParamBase {
  int nvar;
 public:
  ParamBase() {}
  virtual ~ParamBase() {}
  virtual void Display(ostream& os)=0;
  virtual void Help(ostream& os)=0;
  virtual void GraphicsOutput(ostream& os, int define=0)=0;
  virtual void Output(ostream& os)=0;
  virtual void SetStr(const char *)=0;		// Set from string	
  virtual const char *Name()=0;
  virtual int Dump()=0;
};

class VocabularyBase {
 protected:
  DynVector<ParamBase *> ParamList;
 public:	
  VocabularyBase();
  virtual ~VocabularyBase() {}
  Table<ProblemBase> *ProblemTable;
  Table<IntegratorBase> *IntegratorTable;
  Table<RK> *RKIntegratorTable;
	
  ParamBase *Locate(const char *key, int *match_type);
  void ParamAdd(ParamBase *p);
  void Parse(char *s);
  void Assign(const char *s, int warn=1);
  void Sort();
  void List(ostream& os);
  void Dump(ostream& os);
  void GraphicsDump(ostream& os, int define);
  virtual const char *Name()=0;
  virtual const char *Abbrev()=0;
  virtual const char *Directory() {return "";}

  ProblemBase *NewProblem(const char *& key) {
    ProblemBase *p=ProblemTable->Locate(key);
    p->SetAbbrev(key);
    return p;
  }
	
  RK *NewRKIntegrator(const char *& key) {
    char *key2=strdup(key);
    undashify(key,key2);
    const char *key0=key2;
    RK *p=RKIntegratorTable->Locate(key0);
    p->SetAbbrev(key0);
    return p;
  }  

  IntegratorBase *NewIntegrator(const char *& key) {
    char *key2=strdup(key);
    undashify(key,key2);
    const char *key0=key2;
    IntegratorBase *p=IntegratorTable->Locate(key0);
    p->SetAbbrev(key0);
    return p;
  }
	
  virtual const char *FileName(const char* delimiter="", 
			       const char *suffix="");
};

extern VocabularyBase *Vocabulary;

#define METHOD(key) (void) new Entry<key,ProblemBase> (#key,ProblemTable);

#define PLURAL(x) ((x)==1 ? "" : "s")

template<class T>
inline void open_output(T& fout, const char *delimiter, const char *suffix,
			int append)
{
  const char *filename=Vocabulary->FileName(delimiter,suffix);
  if(append) fout.open(filename,fout.app); // Append to end of output file.
  else fout.open(filename);
  if(!fout) msg(ERROR,"Output file %s could not be opened",filename);
  fout.precision(digits);
  errno=0;
}

template<class T>
inline void open_output(T& fout, const char *delimiter, const char *suffix)
{
  open_output(fout,delimiter,suffix,restart);
}	

inline void VocabularyBase::ParamAdd(ParamBase *p)
{
  ParamList[NParam++]=p;
}

#define VOCAB(var,min,max,help) Vocab(&var, #var, min, max, help, 1, 1)
#define VOCAB_NODUMP(var,min,max,help) Vocab(&var, #var, min, max, help, 1, 0)
#define VOCAB_OBSOLETE(var,min,max,help) Vocab(&var, #var, min, max, help, 1,2)

#define VOCAB_ARRAY(var,help) \
Vocab(var, #var, *var-*var, *var-*var, help, \
(int) (sizeof(var)/sizeof(*var)), 1)
	
template<class T>
class Param : public ParamBase {
  const char *name;
  int nvar;
  T *var;
  T min;
  T max;
  const char *help;
  int dump;
 public:
  Param(T *address, const char *s, int n, T min0, T max0, const char *help0,
	int dump0) {
    name=s; nvar=n; var=address; min=min0; max=max0; help=help0; dump=dump0;
    Vocabulary->ParamAdd(this);
  }

  virtual ~Param() {}
  void Set(T x) {int i; for(i=0; i < nvar; i++) var[i]=x;}
  void Set(T *x) {int i; for(i=0; i < nvar; i++) var[i]=x[i];}
	
  void SetStr(const char *);
	
  void Retrieve(T *p) {*p=*var;}
  const char *Name() {return name;}
	
  void Display(ostream& os) {
    if(dump == 2) return; // Don't display obsolete parameters
    os << name << "=";
    for(int i=0; i < nvar-1; i++) os << var[i] << ", ";
    os << var[nvar-1] << endl;
  }
	
  void Help(ostream& os) {
    os << help << endl;
  }
	
  int Dump() {return dump == 1 && nvar == 1;}
  
  void Output(ostream& os) {
    if(dump != 1) return; // Don't dump control parameters
    os << name << "=";
    for(int i=0; i < nvar-1; i++) os << var[i] << ",";
    os << var[nvar-1] << endl;
  }
	
  void Out(ostream& os, const char *type, int define=0, const char *delim="") 
  {
    if(dump != 1) return; // Don't dump control parameters
    switch(define) {
    case 0:
      os << type;
      if(nvar > 1) os << "[]";
      os << " " << name << ";" << endl;
      break;
      
    case 1: 
      os << name << "=";
      if(nvar == 1) {
	os << delim << var[0] << delim << ";";
      }
      else {
	os << "new " << type << "[]{" 
	   << delim << var[0] << delim;
	for(int i=1; i < nvar; i++) os << "," << delim << var[i] << delim;
	os << "};";
      }
      os << endl;
      break;
    }
  }


  void GraphicsOutput(ostream& os, int define=0);
  
  inline int InRange(T);
	
  inline void get_values(const char *arg, T (*rtn)(const char *));
};

template<class T>
inline void Vocab(T *var, const char *s, T min, T max, const char *help, int n,
		  int dump)
{
  (void) new Param<T>(var,s,n,min,max,help,dump);
}

inline void Vocab(unsigned int *var, const char *s, int min, int max,
		  const char *help, int n, int dump)
{
  (void) new Param<unsigned int>(var,s,n,(unsigned int) min,(unsigned int) max,
				 help,dump);
}

inline const char *atosc(const char *s)
{
  return atos(s);
}

template<>
inline void Param<char *>::SetStr(const char *s) {get_values(s,&atos);}
template<>
inline void Param<const char *>::SetStr(const char *s) {get_values(s,&atosc);}
template<>
inline void Param<double>::SetStr(const char *s) {get_values(s,&atof);}
template<>
inline void Param<int>::SetStr(const char *s) {get_values(s,&atoi);}
template<>
inline void Param<unsigned int>::SetStr(const char *s) {get_values(s,&atou);}
template<>
inline void Param<Complex>::SetStr(const char *s) {get_values(s,&atoc);}

template<>
inline int Param<double>::InRange(double value)
{
  return (min == max) || (value >= min && value <= max);
}

template<>
inline int Param<int>::InRange(int value)
{
  return (min == max) || (value >= min && value <= max);
}

template<>
inline int Param<unsigned int>::InRange(unsigned int value)
{
  return (min == max) || (value >= min && value <= max);
}

template<>
inline int Param<const char *>::InRange(const char *) {return 1;}
template<>
inline int Param<char *>::InRange(char *) {return 1;}
template<>
inline int Param<Complex>::InRange(Complex) {return 1;}

template<class T>
inline void Param<T>::get_values(const char *arg, T (*rtn)(const char *))
{
  int i=0;
  char *ptr;
  T value;
  char *optarg=strdup(arg);

  do {
    const char *optarg0=strchr(optarg,')');
    if(!optarg0) optarg0=optarg;
    ptr=strchr(optarg0,',');
    if(ptr) *ptr=0;
    else ptr=optarg;
    if(i >= nvar && param_warn) 
      msg(OVERRIDE_GLOBAL,"Excess initializer \"%s\" to %s",
	  optarg,name);
    else {
      value=(*rtn)(optarg);

      if(InRange(value)) var[i++]=value;
      else if(param_warn) {
	ostringstream buf;
	buf << "Value \"" << value << "\" for "
	    << name << " is invalid. Limits are "
	    << min << " to " << max << ends;
	msg(OVERRIDE_GLOBAL,"%s",buf.str().c_str());
      }
    }	
  } while ((ptr == optarg) ? 0 : (optarg=ptr+1));
}
	
template<>
inline void Param<unsigned int>::GraphicsOutput(ostream& os, int define) {
  Out(os,"int",define);
}

template<>
inline void Param<int>::GraphicsOutput(ostream& os, int define) {
  Out(os,"int",define);
}

template<>
inline void Param<double>::GraphicsOutput(ostream& os, int define) {
  Out(os,"real",define);
}

template<>
inline void Param<Complex>::GraphicsOutput(ostream& os, int define) {
  Out(os,"pair",define);
}

template<>
inline void Param<const char *>::GraphicsOutput(ostream& os, int define) {
  Out(os,"string",define,"\"");
}

#endif
