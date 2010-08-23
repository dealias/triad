import graph;

size(200,150,IgnoreAspect);

string run=getstring("integrator");
file in=input("test/"+run).line();
real[][] a=in.dimension(0,0);
a=transpose(a);

real[] x=a[0];
real[] y=a[1];

draw(graph(x,y),red);

xaxis("$x$",BottomTop,LeftTicks);
yaxis("$y$",LeftRight,RightTicks);
