fisherZTransform<-function(rho)
{
  return (0.5*log((1+rho)/(1-rho),base=exp(1)));
}

#are the degrees of freedome here correct?
compareCorrelations<-function(rho1,rho2,n1,n2,stats=FALSE)
{
  #transform correlation 1
  z1 <- fisherZTransform(rho1);
  #transform correlation 2
  z2 <- fisherZTransform(rho2);
  #estimate variance of z1-z2
  theta <- sqrt(x=(1/(n1-3))+(1/(n2-3)))
  
  #translate to Z distribution
  z = (rho1-rho2)/theta
  
  #calculate p value using t test
  p = 2*pt(-abs(z),df=n1+n2-6)
  
  if(stats)
  {
    return (list(z1=z1, z2=z2, theta=theta, z=z, p=p));
    
  }else {
    return (p);
  }
}