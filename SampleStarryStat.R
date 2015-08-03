genRandomModuleSet<-function(U,M){
  m<-matrix(data = 0,nrow = length(M), ncol = max(M));
  for(i in 1:length(M)){
      m[i,1:M[i]]=sample(U,M[i],replace = FALSE);
  }
  return(m);
}

permutationTest<-function(s,intersectedUniverseSize,lengths1,lengths2,iterations=10000,drawHistogram=FALSE){
  
  totalGenes<-intersectedUniverseSize;
  sampleStarryStat<-vector(length = iterations);
  U<-seq(1:totalGenes);
  for(iteration in 1:iterations)
  {
    clusts1<-genRandomModuleSet(U,lengths1);
    clusts2<-genRandomModuleSet(U,lengths2);
    
    countsMat <- matrix(nrow = dim(clusts1)[1], ncol= dim(clusts2)[1]);
    row.names(countsMat)<-row.names(clusts1)
    colnames(countsMat)<-row.names(clusts2)
    countsMat[,]<-0;
    gMat <- matrix(nrow = dim(clusts1)[1], ncol= dim(clusts2)[1]);
    row.names(gMat)<-row.names(clusts1)
    colnames(gMat)<-row.names(clusts2)
    gMat[,]<-0;
    
    for(i in 1:dim(clusts1)[1])
    {
      for(j in 1:dim(clusts2)[1])
      {
        countsMat[i,j] <- length(intersect(clusts1[i,][clusts1[i,]!=0],clusts2[j,][clusts2[j,]!=0]))
        r<- lengths1[i];
        g<- countsMat[i,j];
        b<- lengths2[j];
        x<- r + b - g;
        gMat[i,j] <- g/x;
      }
    }
    sampleStarryStat[iteration]<-sum(gMat)/length(lengths1)/length(lengths2);
  }

  su<-summary(sampleStarryStat)
  p <- length(sampleStarryStat[sampleStarryStat>s])/length(sampleStarryStat);
  print(p);
  out <- list(pval<-p, summary<-su);
  if(drawHistogram)
  {
    h<-hist(sampleStarryStat);
    lh<-hist(log(sampleStarryStat))
    out <- list(pval<-out$pval, 
                summary<-out$summary, hist=h, lhist=lh);
  }
  return(out);
}

# universeSize<-5948;
# M <- c(76,49,42,26,26,14,13,12,11,11,11,11,11,10,9,8,8,8,7,7,7);
# N <- c(76,49,42,26,26,14,13,12,11,11,11,11,11,10,9,8,8,8,7,7,7);
# #N <- c(51,48,46,34,34,31,26,24,21,17,16,15,15,13,13,12,12,12,12,11,11,11,11,10,10,10,10,10,9,9,9,9,9,9,8,8,8,8,8,8,8,8,8,8,8,8,8,8);
# s = 0.8;
# 
# results<-permutationTest(s,universeSize,lengths1 = M,lengths2=N,drawHistogram = TRUE);
# plot(x=results$hist.mids,y=results$hist.counts,type = "b");
# plot(x=results$lhist.mids,y=results$lhist.counts,type = "b");
