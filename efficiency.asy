import graph;
include getparam;
include averages;
size(200,150,IgnoreAspect);

real[][] a;
real[] data;

while(nextrun()) { 
  file fin=line(input(run+"/stat",comment=""));
  
  string[] names=word(fin);
    
  int field=find(names == "t");
  a=dimension(word(fin,false),0,0);
  a=transpose(a);
  data=a[field];

  draw(graph(a[5],data),p+Pen(n),texify(run));
  
  real[] data0=copy(data);
  data0.delete(0);
  write("slope="+(string) (data[data.length-1]/a[5][a[5].length-1]));
}

if(n > 1) attach(legend(),point(E),20E);

xaxis("CPU",BottomTop,LeftTicks);
yaxis("t",LeftRight,RightTicks(trailingzero));
