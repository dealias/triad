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
  virtual void Display(ostream& os)=0;
  virtual void Help(ostream& os)=0;
  virtual void GraphicsOutput(ostream& os)=0;
  virtual void Output(ostream& os)=0;
  virtual void SetStr(const char *)=0;		// Set from string	
  virtual const char *Name()=0;
};

class VocabularyBase {
 protected:
  DynVector<ParamBase *> ParamList;
 public:	
  VocabularyBase();
  virtual ~VocabularyBase() {}
  Table<ProblemBase> *ProblemTable;
  Table<IntegratorBase> *IntegratorTable;
	
  ParamBase *Locate(const char *key, int *match_type);
  void ParamAdd(ParamBase *p);
  void Parse(char *s);
  void Assign(const char *s, int warn=1);
  void Sort();
  void List(ostream& os);
  void Dump(ostream& os);
  void GraphicsDump(ostream& os);
  virtual const char *Name()=0;
  virtual const char *Abbrev()=0;
  virtual const char *Directory() {return "";}

  ProblemBase *NewProblem(const char *& key) {
    ProblemBase *p=ProblemTable->Locate(key);
    p->SetAbbrev(key);
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
    os << name << "=";
    for(int i=0; i < nvar-1; i++) os << var[i] << ", ";
    os << var[nvar-1] << endl;
  }
	
  void Help(ostream& os) {
    os << help << endl;
  }
	
  void Output(ostream& os) {
    if(!dump) return; // Don't dump control parameters
    os << name << "=";
    for(int i=0; i < nvar-1; i++) os << var[i] << ",";
    os << var[nvar-1] << endl;
  }
	
  void Out(ostream& os, const char *type, const char *delim="") 
  {
    if(!dump) return; // Don't dump control parameters
    if(strcmp(name,"dynamic")==0) // An asy keyword; don't need anyway.
      return;
    if(nvar == 1) {
      os << type << " " << name << "=" << delim << var[0] << delim << ";";
    }
    else {
      os << type << "[] " << name << "={" << delim << var[0] << delim;
      for(int i=1; i < nvar; i++) os << "," << delim << var[i] << delim;
      os << "};";
    }
    os << endl;
  }

#if 0 // For sm 
  void GraphicsOutput(ostream& os) {
    if(!dump) return; // Don't dump control parameters
    if(nvar == 1) os << "define " << name << " \"" << var[0] << "\"";
    else {
      os << "set " << name << "={";
      for(int i=0; i < nvar; i++) os << var[i] << " ";
      os << "}";
    }
    os << endl;
  }
#endif  
	
  void GraphicsOutput(ostream& os);
  
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

inline void Param<char *>::SetStr(const char *s) {get_values(s,&atos);}
inline void Param<const char *>::SetStr(const char *s) {get_values(s,&atosc);}
inline void Param<double>::SetStr(const char *s) {get_values(s,&atof);}
inline void Param<int>::SetStr(const char *s) {get_values(s,&atoi);}
inline void Param<unsigned int>::SetStr(const char *s) {get_values(s,&atou);}
inline void Param<Complex>::SetStr(const char *s) {get_values(s,&atoc);}

inline int Param<double>::InRange(double value)
{
  return (min == max) || (value >= min && value <= max);
}

inline int Param<int>::InRange(int value)
{
  return (min == max) || (value >= min && value <= max);
}

inline int Param<unsigned int>::InRange(unsigned int value)
{
  return (min == max) || (value >= min && value <= max);
}

inline int Param<const char *>::InRange(const char *) {return 1;}
inline int Param<char *>::InRange(char *) {return 1;}
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
	
inline void Param<unsigned int>::GraphicsOutput(ostream& os) {Out(os,"int");}
inline void Param<int>::GraphicsOutput(ostream& os) {Out(os,"int");}
inline void Param<double>::GraphicsOutput(ostream& os) {Out(os,"real");}
inline void Param<Complex>::GraphicsOutput(ostream& os) {Out(os,"pair");}
inline void Param<const char *>::GraphicsOutput(ostream& os) {
  Out(os,"string","\"");
}

#endif
