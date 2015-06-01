
circleIntersectArea<-function (R,r,d)
{
  a<-r*r*acos(x=((d*d+r*r-R*R)/(2*d*r)))
  b<-R*R*acos(x=((d*d+R*R-r*r)/(2*d*R)))
  c<-0.5*sqrt(x=(-d+r+R)*(d+r-R)*(d-r+R)*(d+r+R))
  return (a+b-c);
}

R=1
r=0.9263
d=seq(from=0.55,to=0.6,by=0.0005)

A<-lapply(d,circleIntersectArea,R=R,r=r)