#include <math.h>
#include "Spline.h"

using namespace Array;

typedef array1<double> vector;

int main()
{
  unsigned int n=7;
  vector x(n),y(n);
  
  x[0]=0.0;
  x[1]=0.5;
  x[2]=1.0;
  x[3]=2.0;
  x[4]=2.5;
  x[5]=5.0;
  x[6]=8.0;
  
  y[0]=0.0;
  y[1]=0.5;
  y[2]=1.0;
  y[3]=-1.0;
  y[4]=2.0;
  y[5]=1.0;
  y[6]=1.0;
  
  cout << x << endl;
  cout << y << endl;
  CubicSpline<double,double> S(n,x,y);
  int m=100;
  double delta=(x[n-1]-x[0])/m;
  for(int i=-m/2; i <= 3*m/2; i++) {
    double xi=x[0]+i*delta;
    cout << xi << " " << S.Interpolate(n,x,y,xi) << endl;
  }
    
  return 0;
}
