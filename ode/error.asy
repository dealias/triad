include getparam;

pen p=linewidth(1);

real[] t,y,exact;

file fin=input("test/exact").line();
real[][] a=fin.dimension(0,0);
a=transpose(a);
t=a[0]; exact=a[1];

//real[] exact=exp(-t)*cos(t);

colorPen[2]=heavygreen;

while(nextrun()) {
  file fin=input("test/"+run).line();
  a=fin.dimension(0,0);
  a=transpose(a);
  write(a[1]-exact);
  y=abs(a[1]-exact);
  draw(graph(t,y,t <= 5),p+Pen(n),texify(run));
}

xaxis("$t$",BottomTop,LeftTicks);
yaxis("error",LeftRight,RightTicks);

attach(legend(),point(plain.E),20plain.E);
