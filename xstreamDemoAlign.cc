// g++ -DHAVE_LIBTIRPC -I/usr/include/tirpc xstreamDemoAlign.cc xstream.cc -ltirp

#include <cassert>
#include "xstream.h"

using namespace std;
using namespace xdr;

int main()
{
//  oxstream fout("data");
  memoxstream fout;

//  size_t n=16000004;
  size_t n=16000001;
  size_t sumIn=0;
  if(fout) {
#if 1
    for(size_t i=0; i < n; i++) {
      unsigned char x=i & 0xff;
      sumIn += x;
      fout << (xbyte) x;
    }
#endif
    for(size_t i=0; i < n; i++)
      fout << (double) i << endl;
  }
  fout.close();

//  ixstream fin("data");

  std::vector<unsigned char> data=fout.createCopyOfCurrentData();
  memixstream fin(data);

  xbyte x;
  size_t i=0;
  size_t sumOut=0;
  if(fin) {
#if 1
    for(size_t i=0; i < n; i++) {
      fin >> x;
      sumOut += (unsigned char) x;
    }
#endif
    double v;
    for(size_t i=0; i < n; i++)
      fin >> v;
  }
  assert(sumIn == sumOut);

  return 0;
}
