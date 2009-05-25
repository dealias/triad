#ifndef __Spline_h__
#define __Spline_h__ 1

#include "Array.h"

namespace Array {
  
// Do a binary search of an ordered array of n elements to find an interval
// containing a given key. Return the index corresponding to the left-hand
// endpoint of the matching interval, n-1 if the key is greater than or
// equal to the last element, or n if the key is less than the first element.
template<class X>
unsigned int bintsearch(X key, unsigned int n, const typename array1<X>::opt x)
{
  if(n == 0) return 0;
  if(key < x[0]) return n;
  unsigned int u=n-1;
  if(key >= x[u]) return u;
  unsigned int l=0;
	
  while (l < u) {
    unsigned int i=(l+u)/2;
    if(x[i] <= key && key < x[i+1]) return i;
    if(key < x[i]) u=i;
    else l=i+1;
  }
  ArrayExit("Statement not reachable");
  return 0;
}

// Linearly interpolation of data points (x_i,y_i), where the x_i are listed
// in ascending order. Values outside the available data range are linearly
// extrapolated using the first derivative at the nearest endpoint. 
  
template<class X, class Y>
Y LinearInterpolate(unsigned int n,
		    const typename array1<X>::opt x,
		    const typename array1<Y>::opt y, X x0, unsigned int i) 
{
  if(n == 0) ArrayExit("Zero data points in linear interpolation");
  if(n == 1) return y[0];
  if(i >= n) {
    X dx=x[1]-x[0];
    return y[0]+(y[1]-y[0])/dx*(x0-x[0]);
  }
  if(i == n-1) {
    X dx=x[n-1]-x[n-2];
    return y[n-1]+(y[n-1]-y[n-2])/dx*(x0-x[n-1]);
  }

  X D=x[i+1]-x[i];
  X B=(x0-x[i])/D;
  X A=1.0-B;
  return A*y[i]+B*y[i+1];
}
    
template<class X, class Y>
Y LinearInterpolate(unsigned int n,
		    const typename array1<X>::opt x,
		    const typename array1<Y>::opt y, X x0) 
{
    return LinearInterpolate<X,Y>(n,x,y,x0,bintsearch(x0,n,x));
}

// Compute a natural cubic spline (zero second derivatives at the endpoints)
// of n data points (x_i,y_i), where the x_i are listed in ascending order.
// The constructor computes the second derivatives at each point;
// interpolated values can then be computed with Interpolate. Values
// outside the available data range are linearly extrapolated using
// the first derivative at the nearest endpoint.
  
template<class X, class Y>
class CubicSpline {
protected:
  static double sixth;
  static typename array1<X>::opt work;
  static typename array1<Y>::opt y2;
  static unsigned int size;
public:
  CubicSpline() {size=0;}
  CubicSpline(unsigned int n,
	      const typename array1<X>::opt x,
	      const typename array1<Y>::opt y) {
    if(n < 2) return;
    if(n > size) {
      Reallocate(y2,n);
      if(n > 3) Reallocate(work,n-3);
      size=n;
    }
    y2[0]=y2[n-1]=0.0;
    
    if(n > 2) {
      X lastdx=x[1]-x[0];
      Y lastdy=y[1]-y[0];
      X dx=x[2]-x[1];
      Y dy=y[2]-y[1];
      X temp=3.0/(dx+lastdx);
      y2[1]=(dy/dx-lastdy/lastdx)*temp;
      if(n > 3) work[0]=-sixth*dx*temp;
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
      
      if(n > 3) {
        dx=x[n-1]-x[n-2];
        dy=y[n-1]-y[n-2];
        y2[n-2]=(dy/dx-lastdy/lastdx-sixth*lastdx*y2[n-3])*
          (3.0/((dx+lastdx)+0.5*lastdx*work[n-4]));
      
        for(int i=(int)n-4; i >= 0; i--) y2[i+1] += work[i]*y2[i+2];
      }
    }
  }
  
  Y Interpolate(unsigned int n,
		const typename array1<X>::opt x,
		const typename array1<Y>::opt y, X x0, unsigned int i) {
    if(n == 0) ArrayExit("Zero data points in cubic spline");
    if(n == 1) return y[0];
    if(i >= n) {
      X dx=x[1]-x[0];
      return y[0]+((y[1]-y[0])/dx-sixth*dx*y2[1])*(x0-x[0]);
    }
    if(i == n-1) {
      X dx=x[n-1]-x[n-2];
      return y[n-1]+((y[n-1]-y[n-2])/dx+sixth*dx*y2[n-2])*(x0-x[n-1]);
    }
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

template<class X, class Y>
double CubicSpline<X,Y>::sixth=1.0/6.0;
  
template<class X, class Y>
unsigned int CubicSpline<X,Y>::size=0;

template<class X, class Y>
typename array1<X>::opt CubicSpline<X,Y>::work;

template<class X, class Y>
typename array1<Y>::opt CubicSpline<X,Y>::y2;

}

#endif
