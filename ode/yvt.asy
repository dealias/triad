include getparam;

pen p=linewidth(1);

size(200,150,IgnoreAspect);

while(nextrun()) {
  file fin=input("test/"+run).line();
  real[][] a=fin.dimension(0,0);
  a=transpose(a);
  real[] t=a[0];
  real[] y=a[1];
  draw(graph(t,y),p+Pen(n),texify(run),marker(scale(2.0)*unitcircle,Pen(n),Fill));
}

xaxis("$t$",BottomTop,LeftTicks);
yaxis("$y$",LeftRight,RightTicks);

attach(legend(),point(plain.E),20plain.E);
