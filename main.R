#Single Cell Coexpression Network Project main script
#By: Brian Arand
#November 2014

workingDirectories <- c("/home/barand/Single_Cell_RNA_Seq/","C:/Users/Student/My Research/Single-Cell-RNA-Seq");
for(wd in workingDirectories)
  if(file.exists(wd)){expr=setwd(wd)}

rm(workingDirectories,wd);

print("Library Paths:");
.libPaths();

source("http://bioconductor.org/biocLite.R");

#read in arguments
print("Reading in command line arguments.");
args <- commandArgs(trailingOnly = TRUE);
print(paste0("commandArgs: ",args));

if(length(args) > 0)
{
  #Parse arguments (we expec the form --argName=argValue)
  parseArgs <- function (x) strsplit(sub("^--","",x), "=");
  argsDF <- as.data.frame(do.call("rbind", parseArgs(args)));
  args <- as.character(argsDF$V2)
  names(args) <- argsDF$V1
  rm(argsDF);
}
args<- as.list(args);

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

args$diffCoexFlag <- initializeBooleanArg(arg=args$diffCoexFlag, default=TRUE);
args$diffExprsFlag <- initializeBooleanArg(arg=args$diffExprsFlag, default=FALSE);
args$normFlagRPKM <- initializeBooleanArg(arg=args$normFlagRPKM, default=TRUE);
args$normFlagUbi <- initializeBooleanArg(arg=args$normFlagUbi, default=TRUE);
args$normFlagDESeq <- initializeBooleanArg(arg=args$normFlagDESeq, default=TRUE);
args$normFlagQuant <- initializeBooleanArg(arg=args$normFlagQuant, default=TRUE);
args$QCFlag <- initializeBooleanArg(arg=args$QCFlag, default=FALSE);
args$dataFromRDS <- initializeBooleanArg(arg=args$dataFromRDS, default=FALSE);
args$saveNormalizationRDS <- initializeBooleanArg(arg=args$saveNormalizationRDS, default=FALSE);


#metaData <- cbind(metaData, control=(substr(x=metaData[,5],start=14,stop=16)=='11'));#addControl bool
#metaData[,6] <- gsub(pattern="-",replacement=".", x=metaData[,6]);#replace '-' with '.' to make mapping easier later

if(args$dataFromRDS)
{
  print("Reading normalization.RDS");
  Data <- readRDS("normalization.RDS");
} else
{
  #import data
  print("Importing data files.");
  dataFile <- "Data/GSE57872_GBM_data_matrix.txt/GSE57872_GBM_data_matrix.txt"
  metaFile <- "Data/GSE57872_series_matrix.txt/GSE57872_series_matrix.txt"
  #Metadata
  metaData <- read.table(file=metaFile,header=TRUE,sep="\t",skip=35,fill=TRUE,blank.lines.skip=TRUE);
  metaData <- metaData[,metaData[18,]!="Please note that this sample did not pass the quality control filtering (described in the data processing field), thus was excluded from further data processing"]
  metaData <- metaData[,-1];
  
  Data <- list();
  Data$all <- read.table(file=dataFile);
  indecies <- NULL;
  metaNames <- sapply(X=names(metaData), FUN=function(x){ return(paste0(strsplit(x=paste0(strsplit(x=x,split="\\.")[[1]][c(-1,-2)],collapse='-'),split="_")[[1]][-1],collapse="_"));});
  for(i in 1:dim(Data$all)[2])
  {
    sampleName <- names(Data$all)[i];
    #do some mapping because metadata authors were dumb
    if(grepl("^MGH264_...$",sampleName))
    {
      sampleName<- paste0("MGH26-2_",strsplit(sampleName,"_")[[1]][2]);
    } else if(grepl("CSC$",sampleName))
    {
      sampleName<- paste0(substr(sampleName,1,nchar(sampleName)-3),"GSC");
    } else if(grepl("FCS$",sampleName))
    {
      sampleName<- paste0(substr(sampleName,1,nchar(sampleName)-3),"DGC")
    }
    for(j in 1:(length(metaNames)+1))
    {
      if(j==length(metaNames)+1)
      {
        j=-1;
        break;
      }
      if(sampleName==metaNames[[j]])
        break;
    }
    indecies <- c(indecies, j);
  }
  Data$all <- Data$all[,indecies!=-1];
  indecies <- indecies[indecies!=-1];
  Data$meta <- metaData[,indecies];
  
  #particion by primary cancer sample
  Data$MGH26 <- as.matrix(Data$all[,grepl(x=colnames(Data$meta), pattern="^Single.cell.mRNA.seq_MGH26")])
  #Data$MGH262 <- as.matrix(Data$all[,grepl(x=colnames(Data$meta), pattern="^Single.cell.mRNA.seq_MGH26")])
  Data$MGH28 <- as.matrix(Data$all[,Data$meta[7,]=="Single cell mRNA-seq_MGH28"])
  Data$MGH29 <- as.matrix(Data$all[,Data$meta[7,]=="Single cell mRNA-seq_MGH29"])
  Data$MGH30 <- as.matrix(Data$all[,Data$meta[7,]=="Single cell mRNA-seq_MGH30"])
  Data$MGH31 <- as.matrix(Data$all[,Data$meta[7,]=="Single cell mRNA-seq_MGH31"])
  Data$CSC6 <- as.matrix(Data$all[,Data$meta[7,]=="Single cell mRNA-seq_CSC6"])
  Data$CSC8 <- as.matrix(Data$all[,Data$meta[7,]=="Single cell mRNA-seq_CSC8"])
  Data$Population <- as.matrix(Data$all[,grepl(x=lapply(Data$meta[7,], as.character), pattern="Population mRNA-seq_MGH..")])[,c(3,6,7,8,11)]
  Data$Samples <- as.matrix(cbind(Data$MGH26,Data$MGH28,Data$MGH29,Data$MGH30,Data$MGH31))
  Data$Average <-cbind(apply(Data$MGH26,MARGIN=1,FUN=mean),apply(Data$MGH28,MARGIN=1,FUN=mean),apply(Data$MGH29,MARGIN=1,FUN=mean),apply(Data$MGH30,MARGIN=1,FUN=mean),apply(Data$MGH31,MARGIN=1,FUN=mean));
  
  #hierarchical clustering
  Data$dist <- dist(x=t(Data$Samples),method="euclidian")
  Data$hclust <- hclust(d=Data$dist)
  plot(Data$hclust,cex=0.8)
  
  # vector of colors labelColors = c('red', 'blue', 'darkgreen', 'darkgrey',
  # 'purple')
  # using dendrogram objects
  Data$dendogram = as.dendrogram(Data$hclust)
  labelColors = c("RED", "BLUE", "GREEN", "YELLOW" , "PURPLE", "ORANGE", "PINK")
  # cut dendrogram in 4 clusters
  clusMember = cutree(Data$hclust, 7)
  # function to get color labels
  colLab <- function(n) {
    if (is.leaf(n)) {
      a <- attributes(n)
      labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
      attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
    }
    n
  }
  # using dendrapply
  clusDendro = dendrapply(Data$dendogram, colLab)
  # make plot
  plot(clusDendro, main = "Cool Dendrogram")#, type = "triangle")
  
  #partician by clustering like the original paper.
  
  rm(indecies,metaData,i,j,metaNames, sampleName)

if(args$saveNormalizationRDS)
{
  saveRDS(Data, file="normalization.RDS")
}

if(args$diffExprsFlag)
{
  print("Gene Differential Expression:");
  
  #microarray differential expression
  #may not be working... produces rediculously small p-values
  #source("http://bioconductor.org/biocLite.R");
  biocLite("limma");
  library("limma");
  
  condition=c(rep(x="control",times=Data$conCount), rep(x="cancer", times=Data$canCount));
  #combn <- factor(paste(pData(phenoData)[,1], pData(phenoData)[,2], sep = "_"));
  design <- model.matrix(~condition);# describe model to be fit
  
  fit <- lmFit(Data$ma, design);# fit each probeset to model
  efit <- eBayes(fit);# empirical Bayes adjustment
  write.csv(topTable(efit, coef=2, number=length(efit$p.value)), file=paste0("DiffExpression MicroArray.csv"),quote=FALSE);
  
  #significance comparison
  maPRank <- rank(efit$p.value[,2]);
  rsPRank <- rank(res$padj);
  names(rsPRank) <- row.names(res);
  maPRank <- maPRank[sort(names(maPRank))];
  rsPRank <- rsPRank[sort(names(rsPRank))];
  plot2Groups(GroupA=maPRank, GroupB=rsPRank, xlab="MicroArray Rank", ylab="RNASeq Rank", main="DEG Rank comparison by p-value", file="Comp_DEG_pVal_rank_across_tech.png");
  Rs_PRank <- cor(x=maPRank, y=rsPRank, method="spearman");
  rm(maPRank);
  rm(rsPRank);
  
  #log foldchange comparison
  maFC <- topTable(efit, coef=2, number=length(efit$p.value))[,1];
  names(maFC) <- row.names(topTable(efit, coef=2, number=length(efit$p.value)));
  maFCRank <- rank(maFC);
  rsFC <- res$log2FoldChange;
  names(rsFC) <- row.names(res);
  rsFCRank <- rank(rsFC);
  maFCRank <- maFCRank[sort(names(maFCRank))];
  rsFCRank <- rsFCRank[sort(names(rsFCRank))];
  plot2Groups(GroupA=maFCRank, GroupB=rsFCRank, xlab="MicroArray Rank", ylab="RNASeq Rank", main="Log2 Fold Change Rank comparison", file="Comp_FC_rank_across_tech.png");
  Rs_FCRank <- cor(x=maFCRank, y=rsFCRank, method="spearman");
  rm(maFCRank);
  rm(rsFCRank);
}

#subset by significance
if(args$diffExprsFlag)
{
  subBySig <- FALSE;
  if(subBySig)
  {
    print("Filtering genes by significance:");
    cutoff<- 0.0000001;
    filter <- efit.p.adj<cutoff & !is.na(efit.p.adj);
    maGenes<- as.matrix(cbind(efit.p.adj[filter], Data$ma[filter,]));#attach adjusted p value now
    filter <- res$padj<cutoff & !is.na(res$padj);
    rsGenes<- as.matrix(cbind(res$padj[filter], Data$rs_DESeq[filter,]));#attach adjusted p value now
  } else #subset by rank top X most significant
  {
    print("Filtering genes by fold-change:")
    cutoff <- 12000;
    filter <- rank(efit.p.adj)<=cutoff & !is.na(efit.p.adj);
    maGenes<- as.matrix(cbind(efit.p.adj[filter], Data$ma[filter,]));#attach adjusted p value now
    filter <- rank(res$padj)<=cutoff & !is.na(res$padj);
    rsGenes<- as.matrix(cbind(res$padj[filter], Data$rs_DESeq[filter,]));#attach adjusted p value now
  }
  
  #intersection of technology significant genes and output venn diagram gene lists
  filter <- intersect(row.names(maGenes), row.names(rsGenes));
  maIntGenes <- maGenes[filter,];
  rsIntGenes <- rsGenes[filter,];
  maUniGenes <- maGenes[setdiff(row.names(maGenes),row.names(maIntGenes)),];
  rsUniGenes <- rsGenes[setdiff(row.names(rsGenes),row.names(rsIntGenes)),];
  
  #top genes in each group
  head(sort(maIntGenes[,1]))
  head(sort(rsIntGenes[,1]))
  head(sort(maUniGenes[,1]))
  head(sort(rsUniGenes[,1]))
  #remove p.values from Gene lists
  maGenes <- maGenes[,-1];
  rsGenes <- rsGenes[,-1];
} else
{
  print("No gene filtering performed.")
  maGenes <- Data$ma;
  rsGenes <- Data$rs_DESeq;
  cutoff <- "allGenes";
}

library("ggplot2"); 

if(args$QCFlag)
{
  print("Outputing quality control figures:");
  
  #density_data <- data.frame();
  dataSets<-list(Data$MGH26,Data$MGH262,Data$MGH28,Data$MGH29,Data$MGH30,Data$MGH31, Data$Samples, Data$CSC6, Data$CSC8);
  dataSetNames <- c("MGH26", "MGH262", "MGH28", "MGH29", "MGH30", "MGH31","Samples","CSC6","CSC8");
  #breaks <- 100;
  for(i in 1:length(dataSets))
  {
    dataSet = dataSets[[i]];
    name <- dataSetNames[i];
    
    png(filename=paste0(name,"_normalization_check.png"));
    boxplot(x=dataSet,names=seq(1,dim(dataSet)[2]), outcex=0.5, outpch=20, main="Patient box plots", xlab="Patient", ylab="expression value");
    dev.off();
  }
}

#calculated correlation statistics

correlationHistogram <- function(data, method, breaks=100, file)
{
  corrMat <- cor(x=t(data), method=method, use="complete");
  corrMat[is.na(corrMat)]<-0;
  hist <- hist(x=corrMat,breaks=breaks,plot=FALSE);
  write.csv(x=corrMat,file=file);
  return(list(corrMat=corrMat, hist=hist));
}

library("igraph");

for(method in c("spearman"))
{
  print(paste0("Computing ", method, " correlation for:"));
  
  print("Constructing correlation matricies");
  
  density_data <- data.frame();
  dataSets<-list(Data$MGH26,Data$MGH28,Data$MGH29,Data$MGH30,Data$MGH31);#, Data$Samples, Data$Population, Data$Average);
  dataSetNames <- c("MGH26", "MGH28", "MGH29", "MGH30", "MGH31");#, "Samples", "Population", "Average");
  runningAvgCorMat <- matrix(0,nrow=dim(Data$MGH26)[1],ncol=dim(Data$MGH26)[1]);
  rownames(x=runningAvgCorMat)<- rownames(Data$MGH26);
  colnames(x=runningAvgCorMat)<- rownames(Data$MGH26);
  
  for(i in 1:length(dataSets))
  {
    dataSet = dataSets[[i]];
    name <- dataSetNames[i];
    print(paste0("Calculating correlation histogram for ", name));
    profile <- correlationHistogram(data=dataSet, method=method, file=paste0("../Data/coexpressionNetworks/", name,"_", method,"_int.txt"));
    density_data <- data.frame(cor=c(density_data$cor, profile$hist$mids),
                               density=c(density_data$density, profile$hist$counts/sum(profile$hist$counts)),
                               method=c(density_data$method, rep(x = name, times=length(profile$hist$counts))),stringsAsFactors=FALSE);
    
    #update running Average
    runningAvgCorMat <- runningAvgCorMat * (i-1)/i + profile$corrMat * (1/i);
    
    profile$corrMat<-NULL;
    
    print(paste0("\tTotal Hist counts: ", sum(profile$hist$counts)))
    print(paste0("\tNum genes: ", dim(Data[[i]])[1]))
  }
  
   #plot overlapping histogram of PCC
  print(paste0("Outputting comparative ", method, " histogram:"));
  
  # Density plots
  
  # png(filename=);
  ggplot(data=density_data, aes(x=cor, y=density, group=method, colour=method)) + 
    geom_line(size=1, aes(linetype=method)) +
    ggtitle(paste0(method, " density comparison"));
  ggsave(paste0("Comparative density of network ", method, ".png"));
  #dev.off();
  
  rm(density_data, dataSets, dataSetNames,i, method, profile);
  
  ##################################################################
  #coexpression networks direct comparison
  #maCorrMat <- cor(x=t(Data$ma), method=method, use="complete");
  #rsCorrMat <- cor(x=t(Data$rs_DESeq), method=method, use="complete");
  
  if(args$diffCoexFlag)
  {
    print("Calculating differential coexpression network.");
    
    diffCorrMat <- rsCorrMat - maCorrMat
    
    write.csv(x=diffCorrMat,file=paste0("Data/BRCA/Differential Network ",method,".txt"));
    #calc histogram
    hist <- hist(x=diffCorrMat,breaks=100,plot=FALSE);
    density_data = data.frame("cor"=hist$mids,"density"=hist$counts/sum(hist$counts),"method"=rep(x="diff",times=length(hist$counts)));
    
    #output density distribution
    #png(filename=paste0("Differential Network ", method, " density.png"));
    ggplot(data=density_data, aes(x=cor, y=density, group=method, colour=method)) + 
      geom_line(size=1, aes(linetype=method)) +
      ggtitle(paste0("Density of ", method,"(rs) - ", method,"(ma)"));
    ggsave(paste0("Differential Network ", method, " density.png"));
    #dev.off();
    
    rm(density_data, hist);
    #output venn diagram gene-edge lists
  }
  
  #create iGraph and plot
  
  print("Constructing microArray iGraph.");
  maGraph <-graph.adjacency(adjmatrix=maCorrMat*1000,mode="undirected", weighted=TRUE);
  rm(maCorrMat);
  #add vertex attributes to graph
  for(i in 1:length(V(maGraph)))
  {
    name <- V(maGraph)[i]$name;
    #V(maGraph)[i]$fc <- maFC[name];
    #V(maGraph)[i]$p <- efit.p.adj[name];
  }
  
  print("Outputting microArray iGraph.");
  write.graph(maGraph, file=paste0(method, "_maGraph_", cutoff, ".graphml"), format="graphml" );
  rm(maGraph);
  
  print("Constructing RNASeq iGraph.");
  rsGraph <-graph.adjacency(adjmatrix=rsCorrMat*1000,mode="undirected", weighted=TRUE);
  rm(rsCorrMat);
  for(i in 1:length(V(rsGraph)))
  {
    name <- V(rsGraph)[i]$name;
    #V(rsGraph)[i]$fc <- rsFC[name];
    #V(rsGraph)[i]$p <- res[name,]$padj;
  }
  
  print("Outputting RNASeq iGraph");
  write.graph(rsGraph, file=paste0(method, "_rsGraph_", cutoff, ".graphml"), format="graphml" );
  rm(rsGraph);
  
  if(args$diffCoexFlag)
  {
    
    print("Constructing differential coexpression iGraph");
    diffGraph <-graph.adjacency(adjmatrix=diffCorrMat*1000,mode="undirected", weighted=TRUE);
    rm(diffCorrMat);
    for(i in 1:length(V(diffGraph)))
    {
      name <- V(diffGraph)[i]$name;
      #V(diffGraph)[i]$RS_fc <- rsFC[name];
      #V(diffGraph)[i]$RS_p <- res[name,]$padj;
      #V(diffGraph)[i]$MA_fc <- maFC[name];
      #V(diffGraph)[i]$MA_p <- efit.p.adj[name];
    }
    print("Outputting differential coexpression iGraph");
    write.graph(diffGraph, file=paste0(method, "_diffGraph_", cutoff, ".graphml"), format="graphml" );
    rm(diffGraph);
  }
}
#calculate difference network and fold-change network
#difNet <- rsCorrMat - maCorrMat;
#fcNet <- rsCorrMat / maCorrMat;
#visualize network

##scale free networks via WGCNA
##Load WGCNA package
#library(WGCNA);
##Load additional necessary packages
#library(cluster);
#k=softConnectivity(datE=t(maGenes),power=6);
## Plot a histogram of k and a scale free topology plot
#sizeGrWindow(10,5);
#par(mfrow=c(1,2));
#png(filename="Data/BRCA/maPearson_WGCNA-power6_Hist.png");
#hist(k, main="Connectivity (MArray Pearson pow6)");
#dev.off();
#png(filename="Data/BRCA/maPearson_WGCNA-power6_ScaleFreePlot.png");
#scaleFreePlot(k, main="Check scale free topology (MArray Pearson pow6\n");
#dev.off();

#k=softConnectivity(datE=t(rsGenes),power=6);
## Plot a histogram of k and a scale free topology plot
#sizeGrWindow(10,5);
#par(mfrow=c(1,2));
#png(filename="Data/BRCA/rsPearson_WGCNA-power6_Hist.png");
#hist(k, main="Connectivity (RNASeq Pearson pow6)");
#dev.off();
#png(filename="Data/BRCA/rsPearson_WGCNA-power6_ScaleFreePlot.png");
#scaleFreePlot(k, main="Check scale free topology (RNASeq Pearson pow6)\n");
#dev.off();

quit();
