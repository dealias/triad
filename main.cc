#include "iostream.h"
#include "types.h"

void fft1(Complex data[], unsigned int log2n, int isign);
	
main()
{
	const unsigned int log2n=4;
	unsigned long n=1 << log2n;
	int i;
	Complex data[n];
	
	for(i=0; i < n; i++) {data[i]=Complex(i,i*i);}
	
	fft1(data, log2n, 1);
	for(i=0; i < n; i++) {cout << data[i] << endl;}
	
	for(i=0; i < n; i++) {data[i]=Complex(i,i*i);}
	cout << endl;
	
	fft1(data, 1, 1);
	for(i=0; i < n; i++) {cout << data[i] << endl;}
}
