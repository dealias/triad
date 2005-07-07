#include "options.h"
#include "kernel.h"

#include <sys/stat.h>
#include <unistd.h>

extern VocabularyBase *Vocabulary;
extern double polltime;

void poll()
{
  struct stat buf;
  const char *stopfile=Vocabulary->FileName(dirsep,"STOP");
  
  if(stat(stopfile,&buf) == 0) {
    msg(WARNING_GLOBAL, "Exit forced by presence of file %s",stopfile);
    unlink(stopfile);
    exit(CONTINUE);
  }
}
