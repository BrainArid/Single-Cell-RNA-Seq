genRandomModuleSet<-function(U,M){
  m<-matrix(data = 0,nrow = length(M), ncol = max(M));
  for(i in 1:length(M)){
      m[i,1:M[i]]=sample(U,M[i],replace = FALSE);
  }

  return(m);
}

permutationTest<-function(s,universeSize,M,N,iterations=10000,histogram=FALSE)
  U <- seq(1,5948)
  M <- c(76,49,42,26,26,14,13,12,11,11,11,11,11,10,9,8,8,8,7,7,7);
  N <- c(76,49,42,26,26,14,13,12,11,11,11,11,11,10,9,8,8,8,7,7,7);
  #N <- c(51,48,46,34,34,31,26,24,21,17,16,15,15,13,13,12,12,12,12,11,11,11,11,10,10,10,10,10,9,9,9,9,9,9,8,8,8,8,8,8,8,8,8,8,8,8,8,8);
  totalGenes<-length(U);
  sampleStarryStat<-vector(length = iterations);
  s = 0.8;
  
  for(iteration in 1:iterations)
  {
    clusts1<-genRandomModuleSet(U,M);
    clusts2<-genRandomModuleSet(U,N);
    
    lengths1<-vector(length = dim(clusts1)[1])
    for(i in 1:dim(clusts1)[1])
    {
      lengths1[i]<-length(clusts1[i,][clusts1[i,]!=0])
    }
    
    lengths2<-vector(length = dim(clusts2)[1])
    for(i in 1:dim(clusts2)[1])
    {
      lengths2[i]<-length(clusts2[i,][clusts2[i,]!=0])
    }
    
    #sort clusters' modules by module length
    order1<- order(lengths1,decreasing=TRUE)
    order2<- order(lengths2,decreasing=TRUE)
    clusts1 <- clusts1[order1, ]
    clusts2 <- clusts2[order2, ]
    lengths1<- lengths1[order1]
    lengths2<- lengths2[order2]
    
    countsMat <- matrix(nrow = dim(clusts1)[1], ncol= dim(clusts2)[1]);
    row.names(countsMat)<-row.names(clusts1)
    colnames(countsMat)<-row.names(clusts2)
    countsMat[,]<-0;
    #rMat <- matrix(nrow = dim(clusts1)[1], ncol= dim(clusts2)[1]);
    #row.names(rMat)<-row.names(clusts1)
    #colnames(rMat)<-row.names(clusts2)
    #rMat[,]<-0;
    gMat <- matrix(nrow = dim(clusts1)[1], ncol= dim(clusts2)[1]);
    row.names(gMat)<-row.names(clusts1)
    colnames(gMat)<-row.names(clusts2)
    gMat[,]<-0;
    #bMat <- matrix(nrow = dim(clusts1)[1], ncol= dim(clusts2)[1]);
    #row.names(bMat)<-row.names(clusts1)
    #colnames(bMat)<-row.names(clusts2)
    #bMat[,]<-0;
    
    #if(fisherExact)
    #{
    #  feMat <- matrix(nrow = dim(clusts1)[1], ncol= dim(clusts2)[1],dimnames=list(row.names(clusts1),row.names(clusts2)));
    #  row.names(feMat)<-row.names(clusts1)
    #  colnames(feMat)<-row.names(clusts2)
    #  feMat[,]<-0;
    #}
    
    #print("Processing Overlabs...")
    for(i in 1:dim(clusts1)[1])
    {
      for(j in 1:dim(clusts2)[1])
      {
    #    print(paste0("Cross comparing row ", i, " of ", dim(clusts1)[1]," and column ", j, " of ", dim(clusts2)[1]));
        countsMat[i,j] <- length(intersect(clusts1[i,][clusts1[i,]!=0],clusts2[j,][clusts2[j,]!=0]))
        r<- lengths1[i];
        g<- countsMat[i,j];
        b<- lengths2[j];
        x<- r + b - g;
        #rMat[i,j] <- r/x;
        gMat[i,j] <- g/x;
        #bMat[i,j] <- b/x;
        #if(fisherExact)
        #{
        #  feMat[i,j] <- fisher.test(matrix(data=c(g,r-g,b-g,totalGenes-(r+b-g)),nrow=2,ncol=2))$p.value
        #}
      }
    }
    
    sampleStarryStat[iteration]<-sum(gMat);
  }
  hist(sampleStarryStat);
  summary(sampleStarryStat)
  p <- length(sampleStarryStat[sampleStarryStat>0.8])/length(sampleStarryStat);
}