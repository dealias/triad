import graph;
include getparam;
include averages;
size(200,150,IgnoreAspect);

string fieldname=getstring("field","dt");

if(fieldname == "dt")
  scale(Linear,Log);


while(nextrun()) {
  real[][] a;
  real[] data;
  file fin=input(run+"/stat",comment="").line();
  
  string[] names=fin.word();
  int field=find(names == fieldname);
  if(field < 0) abort("No such field: "+fieldname);
  a=fin.word(false).dimension(0,0);
  a=transpose(a);
  data=a[field];

  real recount=0;
  real[] stops;
  for (int i=1; i < a[0].length; ++i) {
    if (a[0][i] == 1) {
      recount += a[0][i-1];
      stops.push(a[0][i-1]);
    }
    a[0][i] += recount;
  }
	
  draw(graph(a[0],data),p+Pen(n),texify(run));
  for (int i=1; i < stops.length; ++i)
    xequals(stops[i],Pen(n)+dashed);

  real[] data0=copy(data);
  data0.delete(0);
  write(run+":");
  write(" max=",max(data));
  write(" min=",min(data));
  write(" min0=",min(data0));
  write(" mean=",sum(data)/data.length);
  write();
}

if(n > 1) attach(legend(),point(E),20E);

xaxis("it",BottomTop,LeftTicks);
yaxis(texify(fieldname),LeftRight,RightTicks(trailingzero));
