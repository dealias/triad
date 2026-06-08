// g++ -DHAVE_LIBTIRPC -I/usr/include/tirpc xstreamDemo.cc xstream.cc -ltirpc

#include <iostream>
#include "xstream.h"

using namespace std;
using namespace xdr;

int main()
{
  oxstream fout("data");
  if(fout)
    for(size_t i=0; i < 10; i++) fout << 20+i << endl;
  fout.close();

  ixstream fin("data");
  size_t i;
  if(fin)
    while(fin >> i, !fin.eof()) cout << i << endl;

  return 0;
}
