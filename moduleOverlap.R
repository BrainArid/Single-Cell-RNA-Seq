source("SampleStarryStat.R");

#input is expected to be module per row otherwise specify clusts#Columned=TRUE
#totalGenes=22143 comes from GRCh37/hg19 annotation
moduleOverlap <- function(clustsFile1, clusts1Columned=FALSE, clusts1Names=FALSE, 
                          clustsFile2, clusts2Columned=FALSE, clusts2Names=FALSE, 
                          outDir, threshold, clustHeatMap=FALSE, fisherExact=FALSE, 
                          histogram=FALSE, geneUniverse1, geneUniverse2)
{
  #setwd(dir);#setwd("/Users/Brian/Documents/Research/microArray v RNA Seq/BRCA/")
  print("opening files");
  clusts1 <- read.table(file=clustsFile1,sep = "\t",stringsAsFactors = FALSE,fill = TRUE)
  clusts2 <- read.table(file=clustsFile2,sep = "\t",stringsAsFactors = FALSE,fill = TRUE)
  print("files opened");
  
  #rotate if necessary
  if(clusts1Columned)
    clusts1<-t(clusts1)
  if(clusts2Columned)
    clusts2<-t(clusts2)
  #remove first column if module names found
  if(clusts1Names)
  {
    row.names(clusts1)<-factor(clusts1[,1])
    clusts1<- clusts1[,-1] 
  }else
  {
    n<- paste0("module", row.names(clusts1))
    row.names(clusts1)<- factor(n,levels=n)
    rm(n);
  }
  if(clusts2Names)
  {
    row.names(clusts2)<-factor(clusts2[,1])
    clusts2<- clusts2[,-1]
  } else
  {
    n<- paste0("module", row.names(clusts2))
    row.names(clusts2)<- factor(n,levels=n)
    rm(n);
  }
  #remove any columns containing NAs
  clusts1 <- clusts1[,colSums(is.na(clusts1))==0]
  clusts2 <- clusts2[,colSums(is.na(clusts2))==0]
  
  #get intersection of universes
  geneUniverseIntersect<-intersect(x=geneUniverse1,y=geneUniverse2);
  totalGenes<-length(geneUniverseIntersect);
  
  #intersect each clust with universe Intersection
  
  lengths1<-vector(length = dim(clusts1)[1])
  for(i in 1:dim(clusts1)[1])
  {
    temp <- intersect(x=clusts1[i,][clusts1[i,]!=""],
                             y=geneUniverseIntersect);
    clusts1[i,]<-"";
    clusts1[i,1:length(temp)]<-temp;
    lengths1[i]<-length(clusts1[i,][clusts1[i,]!=""]);
  }
  
  lengths2<-vector(length = dim(clusts2)[1])
  for(i in 1:dim(clusts2)[1])
  {
    temp <- intersect(x=clusts2[i,][clusts2[i,]!=""],
                      y=geneUniverseIntersect);
    clusts2[i,]<-"";
    clusts2[i,1:length(temp)]<-temp;
    lengths2[i]<-length(clusts2[i,][clusts2[i,]!=""]);
  }
  remove(temp);
  
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
  rMat <- matrix(nrow = dim(clusts1)[1], ncol= dim(clusts2)[1]);
  row.names(rMat)<-row.names(clusts1)
  colnames(rMat)<-row.names(clusts2)
  rMat[,]<-0;
  gMat <- matrix(nrow = dim(clusts1)[1], ncol= dim(clusts2)[1]);
  row.names(gMat)<-row.names(clusts1)
  colnames(gMat)<-row.names(clusts2)
  gMat[,]<-0;
  bMat <- matrix(nrow = dim(clusts1)[1], ncol= dim(clusts2)[1]);
  row.names(bMat)<-row.names(clusts1)
  colnames(bMat)<-row.names(clusts2)
  bMat[,]<-0;
  if(fisherExact)
  {
    feMat <- matrix(nrow = dim(clusts1)[1], ncol= dim(clusts2)[1],dimnames=list(row.names(clusts1),row.names(clusts2)));
    row.names(feMat)<-row.names(clusts1)
    colnames(feMat)<-row.names(clusts2)
    feMat[,]<-0;
  }
  
  print("Processing Overlaps...")
  for(i in 1:dim(clusts1)[1])
  {
    for(j in 1:dim(clusts2)[1])
    {
      print(paste0("Cross comparing row ", i, " of ", dim(clusts1)[1]," and column ", j, " of ", dim(clusts2)[1]));
      countsMat[i,j] <- length(intersect(clusts1[i,][clusts1[i,]!=""],clusts2[j,][clusts2[j,]!=""]))
      r<- lengths1[i];
      g<- countsMat[i,j];
      b<- lengths2[j];
      x<- r + b - g;
      rMat[i,j] <- r/x;
      gMat[i,j] <- g/x;
      bMat[i,j] <- b/x;
      if(fisherExact)
      {
        feMat[i,j] <- fisher.test(matrix(data=c(g,r-g,b-g,totalGenes-(r+b-g)),nrow=2,ncol=2))$p.value
      }
    }
  }
  
  #visualize
  print("Drawing pretty pictures...")
  library("ggplot2")
  library("reshape")
  library("plyr")
  library("scales")

  rMat.m <- melt(rMat)
  rMat.m$X1 <- factor(rMat.m$X1, levels = row.names(clusts1))
  rMat.m$X2 <- factor(rMat.m$X2, levels = row.names(clusts2))

  gMat.m <- melt(gMat)
  gMat.m$X1 <- factor(gMat.m$X1, levels = row.names(clusts1))
  gMat.m$X2 <- factor(gMat.m$X2, levels = row.names(clusts2))

  #if(clustHeatMap)
  #{
  #  p <- ggplot(gMat.m, aes(y=X1, x=X2))+ 
  #    geom_tile(aes(y=X1,x=X2,fill = rgb(red = 0,green=value,blue=0)), colour = "black")+scale_fill_identity()+
  #    ylab(label=basename(clustsFile1))+
  #    xlab(label=basename(clustsFile2))+
  #    theme(axis.text.x=element_text(angle=90));
  #  ggsave(filename = paste0(outDir, basename(clustsFile1), "_VS_", basename(clustsFile2), "_INTERSECT.png"),plot = p,width=3,height=3)
  #}
   if(fisherExact)
   {
      feMat.m <- melt(round(feMat,5))
      feMat.m$X1 <- factor(feMat.m$X1, levels = row.names(clusts1))
      feMat.m$X2 <- factor(feMat.m$X2, levels = row.names(clusts2))
#       p <- ggplot(feMat.m, aes(y=X1, x=X2))+
#         geom_tile(aes(y=X1,x=X2,fill = rgb(red = 0,green=1-value,blue=0)), colour = "black")+scale_fill_identity()+
#         ylab(label=basename(clustsFile1))+
#         xlab(label=basename(clustsFile2))+
#         theme(axis.text.x=element_text(angle=90));
#       ggsave(filename = paste0(outDir, basename(clustsFile1), "_VS_", basename(clustsFile2), "_FISHER.png"),plot = p,width=3,height=3)
   }
  
  bMat.m <- melt(bMat)
  bMat.m$X1 <- factor(bMat.m$X1, levels = row.names(clusts1))
  bMat.m$X2 <- factor(bMat.m$X2, levels = row.names(clusts2))
  maxSize<-1;
  colMat.m <- data.frame(X1=rMat.m$X1,X2=rMat.m$X2,r=rMat.m$value,g=gMat.m$value,b=bMat.m$value,fe=feMat.m$value,feSize=(feMat.m$value<threshold)*maxSize/2)
  colMat.m$X1 <- factor(colMat.m$X1, levels = row.names(clusts1))
  colMat.m$X2 <- factor(colMat.m$X2, levels = row.names(clusts2))
  if(clustHeatMap)
  {
    p <- ggplot(colMat.m, aes(y=X1, x=X2))+
      geom_tile(aes(y=X1,x=X2,fill = rgb(red = r,green=g,blue=b),width=maxSize,height=maxSize), colour = "black")+
      scale_fill_identity()+
      ylab(label=basename(clustsFile1))+
      xlab(label=basename(clustsFile2))+
      theme(axis.text=element_text(size=9),axis.text.x=element_text(angle=-90,hjust=0,vjust=0.5));
    if(fisherExact)
    {
      p <- p + geom_tile(aes(y=X1,x=X2,fill = rgb(red = 1,green=1,blue=1),width=feSize,height=feSize), colour = "black");
    }
    ggsave(filename = paste0(outDir, basename(clustsFile1), "_VS_", basename(clustsFile2), "_COMPOSITE.png"),plot = p,width=7,height=7)
  }

#   if(histogram)
#   {
#     
#     p <- ggplot(colMat.m[colMat.m$g>0.25,], aes(x=g))+ geom_histogram(binwidth=0.02)+
#       xlab(label="Intersection Size over Union Size")+
#       ylab(label="module count");
#     ggsave(filename = paste0(outDir, basename(clustsFile1), "_VS_", basename(clustsFile2), "_HIST.png"),plot = p,width=3,height=3)
#     
#   }
  

  #permutation test for p-value
  print("Permutation testing to estimate p-value...")
  s<-sum(gMat)/length(lengths1)/length(lengths2);
  permResults<-permutationTest(s,totalGenes,lengths1,lengths2,iterations=10000, drawHistogram = FALSE);
  print("...done.");

  #list white tiles:
  if(dim(colMat.m[colMat.m$fe<threshold,])[1]>0)
  {
    print(colMat.m[colMat.m$fe<threshold,]);
    
    x <- colMat.m[colMat.m$fe<threshold,]$X1
    y <- colMat.m[colMat.m$fe<threshold,]$X2
    
    commonModules <- list();
    commonModules$profiles <- list();
    commonModules$clusts1 <- list();
    commonModules$clusts2 <- list();
    
    relativeRank <- function(L1, L2)
    {      
      l1 <- length(L1);
      l2 <- length(L2);
      m1<-match(L1,L2);
      m2<-match(L2,L1);
      
      return(list(r1=m1-seq(from = 1,to = l1,1),r2=m2-seq(from = 1,to = l2,1)));
    }
    
    print(paste0("Selecting modules for output that meet overlap threshold of ", threshold))
    for(i in 1:length(x))
    {
      print(paste0("Getting overlap ", i, " of ", length(x)) )
      commonModules$profiles[i] <- list(colMat.m[colMat.m$X1==x[i] & colMat.m$X2==y[i],c("X1","X2","r","g","b", "fe")]);
      #if(fisherExact)
      #{
      #  commonModules$profiles[[i]]$fe <- feMat.m[feMat.m$X1==x[i] & feMat.m$X2==y[i],c("value")]
      #}
      commonModules$clusts1[i] <- list(clusts1[x[i],clusts1[x[i],]!=""]);
      commonModules$clusts2[i] <- list(clusts2[y[i],clusts2[y[i],]!=""]);
      if(enrich)
      {
        print(paste0("Calculating enrichment ", i, " of ", length(x)) )
        commonModules$clusts1GO[i] <- list(topGOEnrich(genesOfInterest = commonModules$clusts1[[i]],geneID2GO = geneID2GO));
        commonModules$clusts2GO[i] <- list(topGOEnrich(genesOfInterest = commonModules$clusts2[[i]],geneID2GO = geneID2GO));
      
        relativeClustsLocs <- relativeRank(L1 <- commonModules$clusts1GO[[i]]$BP[,1], L2 <- commonModules$clusts2GO[[i]]$BP[,1])
        commonModules$clusts1GO[[i]]$BP <- cbind(commonModules$clusts1GO[[i]]$BP, relativeRank=relativeClustsLocs$r1);
        commonModules$clusts2GO[[i]]$BP <- cbind(commonModules$clusts2GO[[i]]$BP, relativeRank=relativeClustsLocs$r2);
        relativeClustsLocs <- relativeRank(L1 <- commonModules$clusts1GO[[i]]$MF[,1], L2 <- commonModules$clusts2GO[[i]]$MF[,1])
        commonModules$clusts1GO[[i]]$MF <- cbind(commonModules$clusts1GO[[i]]$MF, relativeRank=relativeClustsLocs$r1);
        commonModules$clusts2GO[[i]]$MF <- cbind(commonModules$clusts2GO[[i]]$MF, relativeRank=relativeClustsLocs$r2);
        relativeClustsLocs <- relativeRank(L1 <- commonModules$clusts1GO[[i]]$CC[,1], L2 <- commonModules$clusts2GO[[i]]$CC[,1])
        commonModules$clusts1GO[[i]]$CC <- cbind(commonModules$clusts1GO[[i]]$CC, relativeRank=relativeClustsLocs$r1);
        commonModules$clusts2GO[[i]]$CC <- cbind(commonModules$clusts2GO[[i]]$CC, relativeRank=relativeClustsLocs$r2);
      }
    }
    
    #and thier profiles:
    gSortedOrder <- order(matrix(data=unlist(commonModules$profiles),nrow=length(commonModules$profiles),ncol=6,byrow=TRUE)[,6],decreasing=TRUE)
    commonModules <- list(profiles=commonModules$profiles[gSortedOrder],
                          clusts1=commonModules$clusts1[gSortedOrder],
                          clusts2=commonModules$clusts2[gSortedOrder],
                          stat=s,
                          pval=permResults$pval);
    return(commonModules);
  }
  return(NULL);
}

print("Reading in command line arguments.");
args <- commandArgs(trailingOnly = TRUE);
print(paste0("commandArgs: ",args));

if(length(args) > 0)
{
  #Parse arguments (we expec the form --argName=argValue)
  parseArgs <- function (x) 
  {
    s<- unlist(strsplit(sub("^--","",x), "="));
    return(list(V1=s[1],V2=paste(s[-1],collapse = "=")))
  }
  argsDF <- as.data.frame(do.call("rbind", lapply(X = args,FUN = parseArgs)));
  args <- as.list(argsDF$V2)
  names(args) <- argsDF$V1
  rm(argsDF)
} else
{
  args <- list();
}

print(paste0("commandArgs: ",args));

#initialize arguments if 
initializeBooleanArg <- function(arg, default){
  if(is.null(arg))
  {
    arg <- default;
  } else if(is.character(arg))
  {
    arg <- as.logical(arg);
  }
  return(arg);
}

initializeStringArg <- function(arg, default){
  if(is.null(arg))
  {
    arg <- default;
  } else if(!is.character(arg))
  {
    arg <- as.character(arg);
  }
  return(arg);
}

initializeFloatArg <- function(arg, default){
  if(is.null(arg))
  {
    arg <- default;
  } else if(!is.numeric(arg))
  {
    arg <- as.numeric(arg);
  }
  return(arg);
}

initializeIntArg <- function(arg, default){
  if(is.null(arg))
  {
    arg <- default;
  } else if(!is.integer(arg))
  {
    arg <- as.integer(arg);
  }
  return(arg);
}

args$dir1 <- initializeStringArg(arg=args$dir1, default="../Data/codensedModules/frequencyNetwork_0.3/");
args$dir2 <- initializeStringArg(arg=args$dir2, default="../Data/codensedModules/GSE48865_norm_TPM_spearman_filt_0.0002/");
args$clustsFile1 <- initializeStringArg(arg=args$clustsFile1, default="GSE57872StarryMap.txt");
args$clusts1Columned <- initializeBooleanArg(arg=args$clusts1Columned, default=TRUE);
args$clusts1Names <- initializeBooleanArg(arg=args$clusts1Names, default=TRUE);
args$clustsFile2 <- initializeStringArg(arg=args$clustsFile2, default="GSE48865StarryMap.txt");
args$clusts2Columned <- initializeBooleanArg(arg=args$clusts2Columned, default=TRUE);
args$clusts2Names <- initializeBooleanArg(arg=args$clusts2Names, default=TRUE);
args$outDir <- initializeStringArg(arg=args$outDir, default="../../Figures");
args$threshold <- initializeFloatArg(arg=args$threshold, default=0.001);
args$clustHeatMap <- initializeBooleanArg(arg=args$clustHeatMap, default=TRUE);
args$fisherExact <- initializeBooleanArg(arg=args$enrich, default=TRUE);
args$histogram <- initializeBooleanArg(arg=args$histogram, default=TRUE);
args$geneUniverseFile1 <- initializeStringArg(arg=args$geneUniverseFile1, default="../Data/coexpressionNetworks/geneOrder.txt");
args$geneUniverseFile2 <- initializeStringArg(arg=args$geneUniverseFile2, default="../Data/coexpressionNetworks/GSE48865_geneOrder.txt");

if(DEBUG)
{
  dir1<-args$dir1
  dir2<-args$dir2
  clustsFile1<-paste0(args$dir1, args$clustsFile1)
  clusts1Columned<-args$clusts1Columned
  clusts1Names<-args$clusts1Names
  clustsFile2<-paste0(args$dir2, args$clustsFile2)
  clusts2Columned<-args$clusts2Columned
  clusts2Names<-args$clusts2Names
  outDir<-args$outDir
  threshold<-args$threshold
  clustHeatMap<-args$clustHeatMap
  fisherExact<-args$fisherExact
  histogram<-args$histogram
}

geneUniverse1 <- as.character(read.table(file=args$geneUniverseFile1,skip=1,header=FALSE)[,2]);
geneUniverse2 <- as.character(read.table(file=args$geneUniverseFile2,skip=1,header=FALSE)[,2]);

moduleFiles1 <- list.files(pattern="G=5_E=3_D=0\\.[6]_Q=0\\.4_B=0\\.4_S=80_C=0\\.6\\.modulesFO\\.hgnc\\.david\\.module$",path=args$dir1)
moduleFiles2 <- list.files(pattern="G=5_E=1_D=0\\.[4-8]_Q=0\\.4_B=0\\.4_S=80_C=0\\.6\\.modulesFO\\.hgnc\\.david\\.module$",path=args$dir2)
moduleFile1<-"../Data/codensedModules/frequencyNetwork_0.3/G=5_E=3_D=0.6_Q=0.4_B=0.4_S=80_C=0.6.modulesFO.hgnc.david.module"
moduleFile2<-"../Data/codensedModules/GSE48865_norm_TPM_spearman_filt_0.0002/G=5_E=1_D=0.6_Q=0.4_B=0.4_S=80_C=0.6.modulesFO.hgnc.david.module"

for(moduleFile1 in moduleFiles1)
{
  for(moduleFile2 in moduleFiles2)
  {
    commonModules <-moduleOverlap(paste0(args$dir1, moduleFile1), args$clusts1Columned, args$clusts1Names, paste0(args$dir2, moduleFile2), args$clusts2Columned, args$clusts2Names, args$outDir, args$threshold, args$clustHeatMap, args$histogram, geneUniverse1, geneUniverse2);
  }
}

#...
commonModules <-moduleOverlap(paste0(args$dir1, args$clustsFile1),args$clusts1Columned, args$clusts1Names, paste0(args$dir2, args$clustsFile2), args$clusts2Columned, args$clusts2Names, args$outDir, args$threshold, args$clustHeatMap, args$histogram, args$enrich);

fileName<-paste0(args$outDir, basename(args$clustsFile1), "_VS_", basename(args$clustsFile2), "_MATCHING_MODULES.csv");
if(is.null(commonModules))
{
  write("No consensus.",fileName);
} else 
{
  unlistAndWrite <- function(x, index, file, append=TRUE, ncolumns=1000, sep=";")
  {
      write(unlist(x[[index]]),file,ncolumns=ncolumns,append = TRUE, sep=sep);
  }
  
  if(args$enrich)
  {
    GOCharts1 <- commonModules$clusts1GO
    GOCharts2 <- commonModules$clusts2GO
    commonModules$clusts1GO <- NULL;
    commonModules$clusts2GO <- NULL;
  }
  
  write(paste0(basename(args$clustsFile1), ";", basename(args$clustsFile2)),fileName,append=FALSE);#FALSE to clear file contents
  write(paste0("number of matching modules:,",length(commonModules[[1]])),fileName,append=TRUE);
  write("X1 index;X2 index;r;g;b;fe",fileName,append=TRUE);

    for(i in 1:length(commonModules[[1]]))
    {
      lapply(commonModules, unlistAndWrite, i, fileName);
      write(x = paste("\n",args$clustsFile1,";;;;;;;;;;",args$clustsFile2,"\n"),fileName,append = TRUE);
      if(args$enrich)
      {
        write(x = "BIOLOGICAL PROCESS",fileName,append = TRUE);
        write.table(x=cbind(GOCharts1[[i]]$BP,x=rep(x="",times = length(GOCharts1[[i]]$BP[,1])),GOCharts2[[i]]$BP),file = fileName,append = TRUE,sep = ";",quote = FALSE,row.names = FALSE,col.names = TRUE);
        write(x = "MOLECULAR FUNCTION",fileName,append = TRUE);
        write.table(x=cbind(GOCharts1[[i]]$MF,x=rep(x="",times = length(GOCharts1[[i]]$MF[,1])),GOCharts2[[i]]$MF),file = fileName,append = TRUE,sep = ";",quote = FALSE,row.names = FALSE,col.names = TRUE);
        write(x = "CELL CYCLE",fileName,append = TRUE);
        write.table(x=cbind(GOCharts1[[i]]$CC,x=rep(x="",times = length(GOCharts1[[i]]$CC[,1])),GOCharts2[[i]]$CC),file = fileName,append = TRUE,sep = ";",quote = FALSE,row.names = FALSE,col.names = TRUE);
        write(x = "\n",fileName,append = TRUE);
      }
  }
}
