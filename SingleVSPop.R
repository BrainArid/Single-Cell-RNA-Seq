
Data <- readRDS("singleVSPop_Data.RDS");

intersectGenes <- intersect(row.names(Data$Average),row.names(Data$GSE48865_norm_filt_TPM))
unionGenes <- union(row.names(Data$Average),row.names(Data$GSE48865_norm_filt_TPM))

library("ggplot2");
require(grid);
plotPopVSc <- function(popData,scData,p,thresh,popName,scName,title,filename){

  Significant=rep(x="<1",times=length(popData))
  thresh<-sort(thresh,decreasing=TRUE);
  for(t in thresh)
  {
    Significant[p<t]<-paste0("<",t);
  }
  
  gs.pal <- colorRampPalette(c("#f8766d","#ffd682","#00bfc4"),space="Lab")
  
  myTheme<-   theme(panel.margin=unit(x=0.1,units="in"),
                    axis.ticks.margin=unit(x=0.1,units="in"),
                    axis.text=element_text(colour="BLACK"),
                    axis.title=element_text(face="bold"),
                    title=element_text(face="bold"));
  
  p<-ggplot(data=data.frame(
    Population=popData,
    SingleCell=scData, 
    Significance=factor(Significant,levels=rev(sort(unique(Significant))),ordered=TRUE)),
            aes(x=Population, y=SingleCell, color=Significance)) +
    scale_colour_manual(values=rev(gs.pal(length(thresh)+1))) +
    #guides(fill = guide_legend(reverse = TRUE)) +
    geom_point() +
    xlab(popName)+ 
    ylab(scName) +
    myTheme +
    ggtitle(title);
  #p;
  
  ggsave(filename=filename,plot=p);
}

#calculate the original square matrix coordinates from the linearized index of an upper triangular matrix
row_index <- function( i, m ){
  i<-i-1;
  row <- (-2*m - 1 + sqrt( (4*m*(m+1) - 8*i - 7) )) / -2;
  if( row == floor(row) ) row <- row -1;
  return( floor(row)+1);
}

column_index <- function( i, m ){
  i<-i-1;
  row = row_index( i+1, m)-1;
  return( (i - m * row + row*(row+1) / 2)+1);
}

#coexpression stuff
#Read in population coexpression network
filename="Data/coexpressionNetworks/GSE48865_norm_TPM_spearman_int.txt"
Data$GSE48865_coex  <- as.matrix(read.table(file=filename,sep=",",row.names=1,header=TRUE))
#Subset to intersection of genes
Data$GSE48865_coex <- Data$GSE48865_coex[match(intersectGenes,row.names(Data$GSE48865_coex)),match(intersectGenes,row.names(Data$GSE48865_coex))];
#consider only upperTri
Data$GSE48865_coex <- Data$GSE48865_coex[upper.tri(Data$GSE48865_coex)]
popName<-"GSE48865_Population";
popCount<-59;

#Analyze 1 single-cell tumor at a time
scFilenames=c("Data/coexpressionNetworks/MGH26_spearman_int.txt","Data/coexpressionNetworks/MGH28_spearman_int.txt","Data/coexpressionNetworks/MGH29_spearman_int.txt","Data/coexpressionNetworks/MGH30_spearman_int.txt","Data/coexpressionNetworks/MGH31_spearman_int.txt");
scNames=c("GSE57872_MGH26","GSE57872_MGH28","GSE57872_MGH29","GSE57872_MGH30","GSE57872_MGH31");
scCounts=c(dim(Data$MGH26)[2],dim(Data$MGH28)[2],dim(Data$MGH29)[2],dim(Data$MGH30)[2],dim(Data$MGH31)[2]);

#p-value upper tri
len <- length(Data$GSE48865_coex)
pValues <- vector(length=len,mode="numeric");

source("../Single Cell RNA Seq/Single-Cell-RNA-Seq/fisherZTransform.R");
for(i in 2:length(scFilenames))
{
  #Read in next single-cell coexpression network
  scFilename<-scFilenames[i];
  scName<-scNames[i];
  scCount<-scCounts[i];
    
  Data$GSE57872_coex  <- as.matrix(read.table(file=scFilename,sep=",",row.names=1,header=TRUE))
  #Subset to intersection of genes
  Data$GSE57872_coex <- Data$GSE57872_coex[match(intersectGenes,row.names(Data$GSE57872_coex)),match(intersectGenes,row.names(Data$GSE57872_coex))];
  #consider only upperTri
  Data$GSE57872_coex <- Data$GSE57872_coex[upper.tri(Data$GSE57872_coex)];

  #p-value per point
  pValues <- mapply(compareCorrelations, Data$GSE48865_coex, Data$GSE57872_coex, rep(popCount,len),rep(scCount,len));
    
  #plot points with p-value as color
  thresh <- c(0.05, 0.01)
  plotPopVSc(popData=Data$GSE48865_coex,
            scData=Data$GSE57872_coex,
            p=pValues,
            thresh=thresh,
            popName=popName,
            scName=scName,
            title=paste0("Per-gene-pair Co-Expression Comparison"),
            filename=paste0("Coexpression_PopVSc_",scName,".png"));
  
  sig_i<-which(pValues<thresh[1]);
  sig_r<-mapply(i=sig_i,m=rep(length(intersectGenes)-1,length(sig_i)),FUN=row_index);
  sig_c<-mapply(i=sig_i,m=rep(length(intersectGenes)-1,length(sig_i)),FUN=column_index)+1;
  sigEdges <- data.frame(gene1=intersectGenes[sig_r], 
                         gene2=intersectGenes[sig_c], 
                         pvalue=pValues[sig_i]);
  
  pValues <- p.adjust(pValues,method="bonferroni")
  
  sigEdges$bonferroni <- pValues[sig_i];
  #sort
  sigEdges<-sigEdges[order(sigEdges$bonferroni,decreasing=FALSE),];
  #output to file
  write.table(x=sigEdges,file=paste0("Coexpression_PopVSc_sigEdges_",scName,".txt"),quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE);
  
  #plot only points that meet signifigance bonferoni
  plotPopVSc(popData=Data$GSE48865_coex,
             scData=Data$GSE57872_coex,
             p=pValues,
             thresh=thresh,
             popName=popName,
             scName=scName,
             title=paste0("Per-gene-pair Co-Expression Comparison (Bonferroni)"),
             filename=paste0("Coexpression_PopVSc_",scName,"_Bonferroni.png"));
  
  png(filename=paste0("Coexpression_PopVSc_",scName,"_Smooth.png"));
  smoothScatter(x=as.vector(Data$GSE48865_coex),y=as.vector(Data$GSE57872_coex),
                xlab=popName,
                ylab=scName, 
                main="Per-gene pair coexpression Comparison");
  dev.off();
}
