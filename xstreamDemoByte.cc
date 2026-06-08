// g++ -DHAVE_LIBTIRPC -I/usr/include/tirpc xstreamDemoByte.cc xstream.cc -ltirp

#include <cassert>
#include "xstream.h"

using namespace std;
using namespace xdr;

int main()
{
//  oxstream fout("data");
  memoxstream fout;

  size_t n=5294967296;
  n=10;
  size_t sumIn=0;
  if(fout) {
    for(size_t i=0; i < n; i++) {
      unsigned char x=i & 0xff;
      sumIn += x;
      fout << (xbyte) x;
    }
  }
  fout.close();

//  ixstream fin("data");

  std::vector<unsigned char> data=fout.createCopyOfCurrentData();
  memixstream fin(data);

  xbyte x;
  size_t i=0;
  size_t sumOut=0;
  if(fin) {
    while(fin >> x, !fin.eof()) {
      sumOut += (unsigned char) x;
      ++i;
    }
  }
  assert(sumIn == sumOut);

  return 0;
}
