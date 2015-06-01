aggregateRSEM.TPM <- function(addresses)
{
  sample <- read.table(file=addresses[[1]],sep="\t",header=TRUE);
  numFeatures <- length(addresses);
  numRows<-dim(sample)[1];
  Mat <- matrix(data = rep(x=0,times=numFeatures*numRows),nrow = numRows,ncol=numFeatures)
  Mat[,1]<-sample$TPM
  
  if(length(addresses > 1))
  {
    for(i in 2:length(addresses))
    {
      sample <- read.table(file=addresses[[i]],sep="\t",header=TRUE);
      Mat[,i]<-sample$TPM
    }
  }
  
  row.names(Mat) <- as.character(sample$gene_id)
  return(Mat);
}

addresses <-  list.files(pattern = ".rsem.genes.results$")
addresses<- addresses[-60]
Data <- list();
Data$TPM <- aggregateRSEM.TPM(addresses = addresses)

Data$avgGeneTPM <- apply(X = log2(Data$TPM+1),MARGIN = 1, FUN = mean)
Data$normTPM <- log2(Data$TPM+1) - Data$avgGeneTPM
#boxplot.matrix(Data$TPM)
##boxplot.matrix(Data$TPM,use.cols = FALSE)
#boxplot.matrix(Data$normTPM)
##boxplot.matrix(Data$normTMP,use.cols = FALSE)

#filter out genes that are less than some threshold
genesFilter <- apply(X=Data$TPM,MARGIN = 1, FUN=function(x){return( mean(log2(x+1))>4.5)})
genesFilter <- genesFilter
h<-hist(x=apply(X=Data$TPM,MARGIN = 1, FUN=function(x){return( mean(log2(x+1)))}),breaks = 100,xlab = "Average log2 TPM",main = "Distribution Average log2")
cummulativeDist <- sum(h$counts)-cumsum(h$counts)

library("ggplot2");

p<-qplot(y=cummulativeDist,x=h$mids,geom = c("point","line"),
      main="Gene count after filtering", 
      ylab="Gene count",
      xlab="Average log2 TPM threshold")
ggsave(filename="GSE48865_TPMfilterGeneCount.png",path="../../../Figures/",plot = p,width = 4,height=4)

Data$filt <- Data$normTPM[genesFilter,]
write.table(x=row.names(Data$filt),file="../../../Data/coexpressionNetworks/GSE48865_geneOrger.txt",quote=FALSE)


correlationHistogram <- function(data, method, breaks=100, file)
{
  corrMat <- cor(x=t(data), method=method, use="complete");
  corrMat[is.na(corrMat)]<-0;
  hist <- hist(x=corrMat,breaks=breaks,plot=FALSE);
  write.csv(x=corrMat,file=file);
  return(list(corrMat=corrMat, hist=hist));
}



method <- "spearman";
density_data <- data.frame();
name <- "GSE48865_norm_TPM"
print(paste0("Calculating correlation histogram for ", name));
profile <- correlationHistogram(data=Data$filt, method=method, file=paste0("../../../Data/coexpressionNetworks/", name,"_", method,"_int.txt"));
density_data <- data.frame(cor=c(density_data$cor, profile$hist$mids),
                           density=c(density_data$density, profile$hist$counts/sum(profile$hist$counts)),
                           method=c(density_data$method, rep(x = name, times=length(profile$hist$counts))),stringsAsFactors=FALSE);

#plot overlapping histogram of PCC
print(paste0("Outputting comparative ", method, " histogram:"));

# Density plots

# png(filename=);
ggplot(data=density_data, aes(x=cor, y=density, group=method, colour=method)) + 
  geom_line(size=1, aes(linetype=method)) +
  ggtitle(paste0(method, " density comparison"));
ggsave(paste0("Comparative density of network ", method, ".png"));

rm(density_data, dataSets, dataSetNames,i, method, profile);