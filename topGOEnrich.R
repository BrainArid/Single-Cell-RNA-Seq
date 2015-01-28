topGo_get_geneID2GO <- function()
{
  #source("http://bioconductor.org/biocLite.R")
  #biocLite("topGO")
  library("topGO")
  
  #read geneOntologyAnnotation
  # currently this list contains doubled-up genes like: UNQ5830/PRO19650/PRO19816 and MNB/DYRK
  # if you don't have this file get it here: http://geneontology.org/page/download-annotations
  temp <- read.table(file = "..//Data//Gene_GO_Annotation//gene_association.goa_human",blank.lines.skip = TRUE, comment.char = "!",header = FALSE,sep="\t",fill=TRUE)
  temp <- data.frame(Gene=temp[,3],GOTerm=temp[,5])
  allGenes <- unique(temp[,1])
  geneID2GO <- list()
  for(i in 1:length(allGenes))
    geneID2GO[[i]]<- as.character(temp[temp$Gene==allGenes[i],2]);
  names(geneID2GO)<-allGenes;
  rm(temp);
  rm(allGenes)
  return(geneID2GO);
}

topGOEnrich <- function(genesOfInterest, geneID2GO)
{ 
  #construct geneList contains 0 and 1 for not-interesting genes and interesting genes
  geneNames<-names(geneID2GO)
  geneList <- factor(as.integer(geneNames %in% genesOfInterest))
  names(geneList) <- geneNames
  str(geneList)
  rm(geneNames)
  
  #construct three part enrichment objects
  GO_BP <- new("topGOdata", description="Biological Process enrichment", 
               ontology = "BP", allGenes = geneList,
               annot = annFUN.gene2GO, gene2GO = geneID2GO)
  
  GO_MF <- new("topGOdata", description="Molecular Function enrichment", 
               ontology = "MF", allGenes = geneList,
               annot = annFUN.gene2GO, gene2GO = geneID2GO)
  
  GO_CC <- new("topGOdata", description="Cell Cycle enrichment", 
                ontology = "CC", allGenes = geneList,
                annot = annFUN.gene2GO, gene2GO = geneID2GO)
  
  GO_BP.results.cf <- runTest(GO_BP,algorithm="classic",statistic = "fisher");
  GO_MF.results.cf <- runTest(GO_MF,algorithm="classic",statistic = "fisher");
  GO_CC.results.cf <- runTest(GO_CC,algorithm="classic",statistic = "fisher");
  GO_BP.results.ct <- runTest(GO_BP,algorithm="classic",statistic = "t");
  GO_MF.results.ct <- runTest(GO_MF,algorithm="classic",statistic = "t");
  GO_CC.results.ct <- runTest(GO_CC,algorithm="classic",statistic = "t");
  GO_BP.results.w1f <- runTest(GO_BP,algorithm="weight01",statistic = "fisher");
  GO_MF.results.w1f <- runTest(GO_MF,algorithm="weight01",statistic = "fisher");
  GO_CC.results.w1f <- runTest(GO_CC,algorithm="weight01",statistic = "fisher");
  
  topGOsCount<-15;
  GO_BP.allRes <- GenTable(GO_BP, classicFisher=GO_BP.results.cf, classicT=GO_BP.results.ct, weight01Fisher=GO_BP.results.w1f, orderBy="weight01Fisher",topNodes=topGOsCount);
  GO_MF.allRes <- GenTable(GO_MF, classicFisher=GO_MF.results.cf, classicT=GO_MF.results.ct, weight01Fisher=GO_MF.results.w1f, orderBy="weight01Fisher",topNodes=topGOsCount);
  GO_CC.allRes <- GenTable(GO_CC, classicFisher=GO_CC.results.cf, classicT=GO_CC.results.ct, weight01Fisher=GO_CC.results.w1f, orderBy="weight01Fisher",topNodes=topGOsCount);

  return(list(BP=GO_BP.allRes,MF=GO_MF.allRes,CC=GO_CC.allRes));
}

#biomart to grab all map hgnc to GO terms
#  biocLite("biomaRt")
#  library("biomaRt")

#   ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
#   go=c("GO:0006306")
#   getBM(attributes="hgnc_symbol", "hgnc_id"
#         filters=c("go", "chromosome_name"),
#         values=list(go, chrom), mart=ensembl)
#   go=c("GO:0051330","GO:0000080","GO:0000114")
#   chrom=c(17,20,"Y")
#  getBM(attributes= c("hgnc_symbol","go_id"),
#        filters="hgnc_symbol", values= genesOfInterest,
#        mart=ensembl);