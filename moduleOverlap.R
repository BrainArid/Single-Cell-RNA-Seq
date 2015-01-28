moduleOverlap <- function(dir, clustsFile1, clustsFile2, outDir, threshold, clustHeatMap)
{
  #setwd(dir);#setwd("/Users/Brian/Documents/Research/microArray v RNA Seq/BRCA/")
  clusts1 <- read.table(file=paste0(dir,clustsFile1),sep = ",",stringsAsFactors = FALSE,fill = TRUE)
  clusts2 <- read.table(file=paste0(dir,clustsFile2),sep = ",",stringsAsFactors = FALSE,fill = TRUE)
  #remove any columns containing NAs
  clusts1 <- clusts1[,colSums(is.na(clusts1))==0]
  clusts2 <- clusts2[,colSums(is.na(clusts2))==0]
  countsMat <- matrix(nrow = dim(clusts1)[1], ncol= dim(clusts2)[1]);
  countsMat[,]<-0;
  rMat <- matrix(nrow = dim(clusts1)[1], ncol= dim(clusts2)[1]);
  rMat[,]<-0;
  gMat <- matrix(nrow = dim(clusts1)[1], ncol= dim(clusts2)[1]);
  gMat[,]<-0;
  bMat <- matrix(nrow = dim(clusts1)[1], ncol= dim(clusts2)[1]);
  bMat[,]<-0;
  
  lengths1<-vector(length = dim(clusts1)[1])
  for(i in 1:dim(clusts1)[1])
  {
    lengths1[i]<-length(clusts1[i,][clusts1[i,]!=""])
  }
  
  lengths2<-vector(length = dim(clusts2)[1])
  for(i in 1:dim(clusts2)[1])
  {
    lengths2[i]<-length(clusts2[i,][clusts2[i,]!=""])
  }
  
  for(i in 1:dim(clusts1)[1])
  {
    for(j in 1:dim(clusts2)[1])
    {
      print(paste0("Cross comparing microarray row ", i, " of ", dim(clusts1)[1]," and column ", j, " of ", dim(clusts2)[1]));
      countsMat[i,j] <- length(intersect(clusts1[i,][clusts1[i,]!=""],clusts2[j,][clusts2[j,]!=""]))
      r<- lengths1[i];
      g<- countsMat[i,j];
      b<- lengths2[j];
      x<- r + b - g;
      rMat[i,j] <- r/x;
      gMat[i,j] <- g/x;
      bMat[i,j] <- b/x;
    }
  }
  
  #visualize
  library("ggplot2")
  library("reshape")
  library("plyr")
  library("scales")
  countsMat.m <- melt(countsMat)
  countsMat.m <- countsMat.m[with(countsMat.m, order(X1, X2)), ]
  countRange <- range(countsMat);
  #countsMat.m <- ddply(countsMat.m, c("value"), transform, rescale = rescale(value,to = c(0,1),from=countRange))
  #p <- ggplot(countsMat.m, aes(y=X1, x=X2)) + geom_tile(aes(y=X1, x=X2, fill = rescale), colour = "black") + scale_fill_gradient(low = "black", high = "white")
  
  rMat.m <- melt(rMat)
  rMat.m <- rMat.m[with(rMat.m, order(X1, X2)), ]
  #countRange <- range(rMat);
  #rMat.m <- ddply(rMat.m, c("value"), transform, rescale = rescale(value,to = c(0,1),from=countRange))
  #p <- ggplot(rMat.m, aes(y=X1, x=X2)) + geom_tile(aes(y=X1,x=X2,fill = rgb(red = value,green=0,blue=0)), colour = "black")+scale_fill_identity();
  
  gMat.m <- melt(gMat)
  gMat.m <- gMat.m[with(gMat.m, order(X1, X2)), ]
  #countRange <- range(gMat);
  #gMat.m <- ddply(gMat.m, c("value"), transform, rescale = rescale(value,to = c(0,1),from=countRange))
  if(clustHeatMap)
  {
    p <- ggplot(gMat.m, aes(y=X1, x=X2)) + geom_tile(aes(y=X1,x=X2,fill = rgb(red = 0,green=value,blue=0)), colour = "black")+scale_fill_identity();
    ggsave(filename = paste0(outDir, basename(clustsFile1), "_VS_", basename(clustsFile2), "_INTERSECT.png"),plot = p)
  }
  
  bMat.m <- melt(bMat)
  bMat.m <- bMat.m[with(bMat.m, order(X1, X2)), ]
  #countRange <- range(bMat);
  #bMat.m <- ddply(bMat.m, c("value"), transform, rescale = rescale(value,to = c(0,1),from=countRange))
  #p <- ggplot(bMat.m, aes(y=X1, x=X2)) + geom_tile(aes(y=X1,x=X2,fill = rgb(red = 0,green=0,blue=value)), colour = "black")+scale_fill_identity();
  
  colMat.m <- data.frame(X1=rMat.m$X1,X2=rMat.m$X2,r=rMat.m$value,g=gMat.m$value,b=bMat.m$value)
  if(clustHeatMap)
  {
    p <- ggplot(colMat.m, aes(y=X1, x=X2)) + geom_tile(aes(y=X1,x=X2,fill = rgb(red = r,green=g,blue=b)), colour = "black")+scale_fill_identity();
    ggsave(filename = paste0(outDir, basename(clustsFile1), "_VS_", basename(clustsFile2), "_COMPOSITE.png"),plot = p)
  }
  
  source("topGOEnrich.R");
  geneID2GO <- topGo_get_geneID2GO();
  
  #list white tiles:
  if(dim(colMat.m[colMat.m$r>threshold & colMat.m$g>threshold & colMat.m$b>threshold,])[1]>0)
  {
    print(colMat.m[colMat.m$r>threshold & colMat.m$g>threshold & colMat.m$b>threshold,]);
    
    x <- colMat.m[colMat.m$r>threshold & colMat.m$g>threshold & colMat.m$b>threshold,]$X1
    y <- colMat.m[colMat.m$r>threshold & colMat.m$g>threshold & colMat.m$b>threshold,]$X2
    
    commonModules <- list();
    commonModules$profiles <- list();
    commonModules$clusts1 <- list();
    commonModules$clusts2 <- list();
    commonModules$clusts1GO <- list();
    commonModules$clusts2GO <- list();
    
    relativeRank <- function(L1, L2)
    {      
      l1 <- length(L1);
      l2 <- length(L2);
      m1<-match(L1,L2);
      m2<-match(L2,L1);
      
      return(list(r1=m1-seq(from = 1,to = l1,1),r2=m2-seq(from = 1,to = l2,1)));
    }
    
    for(i in 1:length(x))
    {
      commonModules$profiles[i] <- list(colMat.m[colMat.m$X1==x[i] & colMat.m$X2==y[i],c("X1","X2","r","g","b")]);
      commonModules$clusts1[i] <- list(clusts1[x[i],clusts1[x[i],]!=""]);
      commonModules$clusts2[i] <- list(clusts2[y[i],clusts2[y[i],]!=""]);
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
    
    #and thier profiles:
    gSortedOrder <- order(matrix(data=unlist(commonModules$profiles),nrow=length(commonModules$profiles),ncol=5,byrow=TRUE)[,4],decreasing=TRUE)
    commonModules <- list(profiles=commonModules$profiles[gSortedOrder],
                          clusts1=commonModules$clusts1[gSortedOrder],
                          clusts2=commonModules$clusts2[gSortedOrder],
                          clusts1GO=commonModules$clusts1GO[gSortedOrder],
                          clusts2GO=commonModules$clusts2GO[gSortedOrder])
    return(commonModules);
  }
  return(NULL)
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

args$dir <- initializeStringArg(arg=args$dir, default="./");
args$clustsFile1 <- initializeStringArg(arg=args$clustsFile1, default="ma_pearson_allGenes_int.txtg=0.90.modules");
args$clustsFile2 <- initializeStringArg(arg=args$clustsFile2, default="rs_DESeq_spearman_allGenes_int.txtg=0.60.modules");
args$outDir <- initializeStringArg(arg=args$outDir, default="out/");
args$threshold <- initializeFloatArg(arg=args$threshold, default=0.70);
args$clustHeatMap <- initializeBooleanArg(arg=args$clustHeatMap, default=FALSE);

commonModules <-moduleOverlap(args$dir, args$clustsFile1, args$clustsFile2, args$outDir, args$threshold, args$clustHeatMap);

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
  
  GOCharts1 <- commonModules$clusts1GO
  GOCharts2 <- commonModules$clusts2GO
  commonModules$clusts1GO <- NULL;
  commonModules$clusts2GO <- NULL;
  
  write(paste0(basename(args$clustsFile1), ";", basename(args$clustsFile2)),fileName,append=FALSE);#FALSE to clear file contents
  write(paste0("number of matching modules:,",length(commonModules[[1]])),fileName,append=TRUE);
  write("X1 index;X2 index;r;g;b",fileName,append=TRUE);
  for(i in 1:length(commonModules[[1]]))
  {
    lapply(commonModules, unlistAndWrite, i, fileName);
    write(x = paste("\n",args$clustsFile1,";;;;;;;;;;",args$clustsFile2,"\n"),fileName,append = TRUE);
    write(x = "BIOLOGICAL PROCESS",fileName,append = TRUE);
    write.table(x=cbind(GOCharts1[[i]]$BP,x=rep(x="",times = length(GOCharts1[[i]]$BP[,1])),GOCharts2[[i]]$BP),file = fileName,append = TRUE,sep = ";",quote = FALSE,row.names = FALSE,col.names = TRUE);
    write(x = "MOLECULAR FUNCTION",fileName,append = TRUE);
    write.table(x=cbind(GOCharts1[[i]]$MF,x=rep(x="",times = length(GOCharts1[[i]]$MF[,1])),GOCharts2[[i]]$MF),file = fileName,append = TRUE,sep = ";",quote = FALSE,row.names = FALSE,col.names = TRUE);
    write(x = "CELL CYCLE",fileName,append = TRUE);
    write.table(x=cbind(GOCharts1[[i]]$CC,x=rep(x="",times = length(GOCharts1[[i]]$CC[,1])),GOCharts2[[i]]$CC),file = fileName,append = TRUE,sep = ";",quote = FALSE,row.names = FALSE,col.names = TRUE);
    write(x = "\n",fileName,append = TRUE);

  }
}
