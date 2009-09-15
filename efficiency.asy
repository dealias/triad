import graph;
include getparam;
include averages;
size(200,150,IgnoreAspect);

while(nextrun()) {
  real[][] a;
  real[] data;
  file fin=input(run+"/stat",comment="").line();
  
  string[] names=fin.word();
    
  int field=find(names == "t");
  a=fin.word(false).dimension(0,0);
  a=transpose(a);
  data=a[field];

  real[] stops;
  for (int i=1; i < a[0].length; ++i) {
    if (a[0][i] == 1) {
      stops.push(a[5][i-1]);
    }
  }

  draw(graph(a[5],data),p+Pen(n),texify(run));

  //  for (int i=0; i < stops.length; ++i)
  //    xequals(stops[i],Pen(n)+dashed);
  
  real[] data0=copy(data);
  data0.delete(0);
  if (a[5].length>1)
    write("slope="+(string) (data[data.length-1]/a[5][a[5].length-1]));
  else
    write("insufficient progress.");
}

if(n > 1) attach(legend(),point(E),20E);

xaxis("CPU",BottomTop,LeftTicks);
yaxis("$t$",LeftRight,RightTicks(trailingzero));
