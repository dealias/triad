#ifndef __Spline_h__
#define __Spline_h__ 1

#include "Matrix.h"
#include "Tridiagonal.h"

namespace Array {
  
// Do a binary search of an ordered array to find an interval containing key.
// Return the index corresponding to the left-hand endpoint of the matching
// interval, or -1 if key is less than the first element.
template<class X>
unsigned int bintsearch(X key, unsigned int n, const typename array1<X>::opt x)
{
  if(key < x[0]) return -1;
  if(key >= x[n-1]) return n-1;
  
  unsigned int l=0;
  unsigned int u=n-1;
	
  while (l < u) {
    unsigned int i=(l+u)/2;
    if(x[i] <= key && key < x[i+1]) return i;
    if(key < x[i]) u=i;
    else l=i+1;
  }
  exit(1); // Statement not reachable
}

// Compute a natural cubic spline (zero second derivatives at the endpoints)
// of n data points (x_i,y_i) where the x_i are listed in ascending order.
// The constructor computes the second derivatives at each point;
// interpolated values can then be computed with Interpolate.
  
template<class Y, class X>
class CubicSpline {
protected:
  static const double sixth=1.0/6.0;
  static typename array1<X>::opt work;
  static typename array1<Y>::opt y2;
  static unsigned int size;
public:
  CubicSpline() {size=0;}
  CubicSpline(unsigned int n,
	      const typename array1<X>::opt x,
	      const typename array1<Y>::opt y) {
    if(n > size) {
      Reallocate(y2,n);
      Reallocate(work,n-3);
      size=n;
    }
    y2[0]=y2[n-1]=0.0;
    
    X lastdx=x[1]-x[0];
    Y lastdy=y[1]-y[0];
    if(n > 2) {
      X dx=x[2]-x[1];
      Y dy=y[2]-y[1];
      X temp=3.0/(dx+lastdx);
      y2[1]=(dy/dx-lastdy/lastdx)*temp;
      work[0]=-sixth*dx*temp;
      lastdx=dx;
      lastdy=dy;
	
      for(unsigned int i=1; i < n-3; i++) {
	X dx=x[i+2]-x[i+1];
	Y dy=y[i+2]-y[i+1];
	X temp=3.0/((dx+lastdx)+0.5*lastdx*work[i-1]);
	y2[i+1]=(dy/dx-lastdy/lastdx-sixth*lastdx*y2[i])*temp;
	work[i]=-sixth*dx*temp;
	lastdx=dx;
	lastdy=dy;
      }
      
      dx=x[n-1]-x[n-2];
      dy=y[n-1]-y[n-2];
      temp=3.0/((dx+lastdx)+0.5*lastdx*work[n-4]);
      y2[n-2]=(dy/dx-lastdy/lastdx-sixth*lastdx*y2[n-3])*temp;

      for(int i=(int)n-4; i >= 0; i--) y2[i+1] += work[i]*y2[i+2];
    }
    return;
  }
  
  Y Interpolate(unsigned int n,
		const typename array1<X>::opt x,
		const typename array1<Y>::opt y, X x0, unsigned int i) {
    if(x0 == x[n-1]) return y[n-1];
    if(i < 0 || i >= n) ArrayExit("Index out of range");
    X D=x[i+1]-x[i];
    X B=(x0-x[i])/D;
    X A=1.0-B;
    D *= sixth*D;
    X C=(A*A*A-A)*D;
    D *=B*B*B-B;
    
    return A*y[i]+B*y[i+1]+C*y2[i]+D*y2[i+1];
  }
  
  Y Interpolate(unsigned int n,
		const typename array1<X>::opt x,
		const typename array1<Y>::opt y, X x0) {
    return Interpolate(n,x,y,x0,bintsearch(x0,n,x));
  }
};

template<class Y, class X>
unsigned int CubicSpline<Y,X>::size=0;

template<class Y, class X>
typename array1<X>::opt CubicSpline<Y,X>::work;

template<class Y, class X>
typename array1<Y>::opt CubicSpline<Y,X>::y2;

}

#endif
