#include "options.h"
#include "kernel.h"
#include "Param.h"

unsigned int NParam=0;
int param_warn=1;

static int ParamCompare(const void *a, const void *b);
static int ParamKeyCompare(const void*key, const void *p, const size_t n);

void VocabularyBase::Sort()
{
  qsort(ParamList,NParam,sizeof(ParamBase *),ParamCompare);
	
  for(unsigned int i=0; i < NParam-1; i++)
    if(ParamCompare(&ParamList[i],&ParamList[i+1]) == 0)
      msg(ERROR_GLOBAL,
	  "Duplicate vocabulary entry: %s",ParamList[i]->Name());
}

ParamBase *VocabularyBase::Locate(const char *key, int *match_type) 
{
  ParamBase **ptr;
  ptr=(ParamBase **) bsearch2(key,ParamList,NParam,sizeof(ParamBase *),
			      ParamKeyCompare,match_type);
  if(ptr) return *ptr;
  else return NULL;
}

void VocabularyBase::List(ostream& os)
{
  for(unsigned int i=0; i < NParam; i++) ParamList[i]->Display(os);
}

void VocabularyBase::Dump(ostream& os)
{
  for(unsigned int i=0; i < NParam; i++) ParamList[i]->Output(os);
}

void VocabularyBase::GraphicsDump(ostream& os)
{
  for(unsigned int i=0; i < NParam; i++) ParamList[i]->GraphicsOutput(os);
}

void VocabularyBase::Parse(char *s)
{
  char *command=strtok(s," \n\t");
  if(command)  {
    Assign(command);
    while((command=strtok(NULL," \n\t"))) {
      Assign(command);
    }
    errno=0;
  }
}	

void VocabularyBase::Assign(const char *key, int warn)
{
  char *buffer=new char[strlen(key)+1];
  char *ptr;
  strcpy(buffer,key);
  int match_type;
  ParamBase *param;
	
  param_warn=warn;
	
  ptr=strchr(buffer,'=');
  if(!ptr) msg(ERROR_GLOBAL,"Invalid assignment: %s",buffer);
  *ptr=0;
	
  param=Locate(buffer,&match_type);
  if(check_match(match_type,"command",buffer,warn)) param->SetStr(++ptr);
  delete [] buffer;
}


static int ParamCompare(const void *a, const void *b)
{
  return strcmp((*(ParamBase **)a)->Name(),(*(ParamBase **)b)->Name());
}

static int ParamKeyCompare(const void* key, const void *p, const size_t n)
{
  return strcmpn((char *) key,(*(ParamBase **)p)->Name(),n);
}
