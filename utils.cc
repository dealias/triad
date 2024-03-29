#include "utils.h"

#include <cerrno>
#include <string>
#include <ctype.h>
#include <time.h>

const char* run="";
ExitCode exit_signal=COMPLETE;

const double pi=PI;
const double twopi=2.0*pi;
const double twopi2=twopi*twopi;

char *upcase(const char *s, char *s2)
{
  char *p,*pstop;
  size_t n;

  n = strlen(s);
  if(s2==NULL) s2=new char[n+1];
  pstop=s2+n;
	
  for(p=s2; p < pstop; p++) *p = (char) toupper(*s++);
  *p='\0';

  return s2;
}

char *downcase(const char *s, char *s2)
{
  char *p,*pstop;
  size_t n;

  n = strlen(s);
  if(s2==NULL) s2=new char[n+1];
  pstop=s2+n;
	
  for(p=s2; p < pstop; p++) *p = (char) tolower(*s++);
  *p='\0';

  return s2;
}

char *dashify(const char *s, char *s2)
{
  return convert(s, '_', '-', s2);
}

char *undashify(const char *s, char *s2)
{
  return convert(s, '-', '_', s2);
}

char *convert(const char *s, char from, char to, char *s2)
{
  char *p;
  if(s2==NULL) {s2=strdup(s);}
  p=s2;
  while((p=strchr((char *)p,from))) {*p++=to;}
  return s2;
}

int RealCompare(const void *a, const void *b)
{
  return isgn(*(Real *) a - *(Real *)b);
}

// Generalized strncmp utility.
int strcmpn(const char *s1, const char *s2, size_t n)
{
  if(n==0) return strcmp(s1,s2);
  else return strncmp(s1,s2,n);
}

// Generalized strncasecmp utility.
int strcasecmpn(const char *s1, const char *s2, size_t n)
{
  if(n==0) return strcasecmp(s1,s2);
  else return strncasecmp(s1,s2,n);
}

// Like bsearch but also set flag match_type to:
//
//	-1 if no match is found,
//       0 if match is ambiguous,
//       1 if match is an unambiguous abbreviation,
//       2 if match is exact.
//
void *bsearch2(register const void *key,
	       register const void *base,
	       size_t nmemb,
	       register size_t size,
	       int (*compar)(const void *, const void *, const size_t),
	       int *match_type)
{
  register size_t l, u, idx, n;
  register const void *p, *p0;
  register int comparison;

  l = 0;
  u = nmemb;
  *match_type=-1;
	
  while (l < u) {
    idx = (l + u) / 2;
    p = (void *) (((const char *) base) + (idx * size));
    comparison = (*compar)((char *) key,p,0);
    if(comparison == 0) {
      *match_type=2;	// Exact match.
      return (void *) p;
    }
			
    n=strlen((char *) key);
			
    if((*compar)(key,p,n) == 0) {
      // Strings agree to n characters; an abbreviation has been found.
      *match_type=1;
				
      // Check for unambiguous match:
      // If the key also matches the previous entry to n characters,
      // it is ambiguous. 
      if(idx > 0) {
	p0 = (void *) (((const char *) base) + ((idx-1) * size));
	if((*compar)((char *) key,p0,n) == 0) *match_type=0; 
      }
				
      // If the key also matches the next entry to n characters,
      // it is ambiguous.
      if(idx < nmemb-1) {
	p0 = (void *) (((const char *) base) + ((idx+1) * size));
	if((*compar)((char *) key,p0,n) == 0) *match_type=0; 
      }
				
      if(*match_type) return (void *) p;
    }

    if(comparison < 0) u=idx;
    else l=idx+1;
  }

  // No match found.
  return (void *) base;
}

int check_match(int match_type, const char *object, const char *key, int warn)
{
  switch(match_type) {
  case -1:
    if(warn) msg(warn == 1 ? OVERRIDE_ : RETRY_,"",0,
		 "Unknown %s: %s",object,key);
    return 0;
  case 0:
    if(warn) msg(warn == 1 ? OVERRIDE_ : RETRY_,"",0,
		 "Ambiguous %s: %s",object,key);
    return 0;
  }
  return 1;
}

size_t atou(const char *s)
{
  return atoi(s);
}
	
Complex atoc(const char *s)
{
  double re,im;
  const char *ptr;
  const char sep=',';

  while(isspace(*s)) s++;
  if(*s != '(') return atof(s);
  if(!(ptr=strchr(s,sep))) return 0.0;
  re=atof(++s);
  if(!(strrchr(s,')'))) return 0.0;
  im=atof(++ptr);
  return Complex(re,im);
}

char *atos(const char *s)
{
  const int n=100;
  char *s2;
	
  s2=new char[n+1];
  strncpy(s2,s,n);
  return s2;
}

double drand_gauss()
{
  double factor,r2,v1,v2;
  static int flag=0;
  static double save;
			  
  if (flag) {
    flag=0;
    return save;
  } else {
    flag=1;
    do {
      v1=2.0*drand()-1.0;
      v2=2.0*drand()-1.0;
      r2=v1*v1+v2*v2;
    } while (r2 >= 1.0 || r2 == 0.0);
    factor=sqrt(-2.0*log(r2)/r2);
    save=v2*factor;
    return v1*factor;
  }
}

Complex crand_gauss()
{
  double r2,v1,v2;
  do {
    v1=2.0*drand()-1.0;
    v2=2.0*drand()-1.0;
    r2=v1*v1+v2*v2;
  } while (r2 >= 1.0 || r2 == 0.0);
  return Complex(v1,v2)*sqrt(-log(r2)/r2);
}

char *output_filename(char *basename, char *suffix)
{
  char *filename;
  suffix=downcase(suffix);
  filename=new char[strlen(basename)+strlen(suffix)+1];
  strcpy(filename,basename);
  strcat(filename,suffix);
  return filename;
}

Complex *out_base;
Real out_re(size_t i) {return out_base[i].re;}
Real out_im(size_t i) {return out_base[i].im;}

