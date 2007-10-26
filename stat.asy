import graph;
include getparam;
include averages;
size(200,150,IgnoreAspect);

string fieldname=getstring("field","dt");

if(fieldname == "dt")
  scale(Linear,Log);

while(nextrun()) {
  file fin=line(input(run+"/stat",comment=""));
  
  string[] names=word(fin);
  int field=find(names == fieldname);
  if(field < 0) abort("No such field: "+fieldname);
  real[][] a=dimension(word(fin,false),0,0);
  a=transpose(a);
  real[] data=a[field];

  draw(graph(a[0],data),p+Pen(n),texify(run));
  
  write("max=",max(data));
  write("min=",min(data));
  write("mean=",sum(data)/data.length);
}

if(n > 1) attach(legend(),point(E),20E);

xaxis("it",BottomTop,LeftTicks);
yaxis(texify(fieldname),LeftRight,RightTicks(trailingzero));
