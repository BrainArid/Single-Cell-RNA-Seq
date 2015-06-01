
Data <- readRDS("singleVSPop_Data.RDS");
Data$meta
meta <- matrix(nrow=430,ncol=length(Data$meta[[1]]));
for(i in 1:430)
{
  meta[i,]<-as.character(Data$meta[[i]])
}

subtypes <- unique(meta[,11])
cancers <- unique(meta[,9])

for(subtype in subtypes)
{
  for(cancer in cancers)
  {
    print(paste(cancer, subtype, " ; count:",length(meta[meta[,11]==subtype & meta[,9]==cancer,11])))
  }
}