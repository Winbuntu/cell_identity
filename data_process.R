library(data.table)
library(randomForest)

#data = fread("GSE59739_DataTable.txt")

data =  fread("GSE81861_Cell_Line_FPKM.csv")
data = as.data.frame(data)

cell.type = sapply(strsplit(colnames(data)[-1], split='__',fixed=TRUE),function(x) x[2])


cell.label = sapply(strsplit(colnames(data)[-1], split='__',fixed=TRUE),function(x) x[3])

gene = sapply(strsplit(data$V1, split='_',fixed=TRUE),function(x) x[2])

gene  = data$V1

data = data[,-1]
#####################################
library(FactoMineR)

data.big = data[ (apply(data,1,function(x)  sum(  x>1   ) > dim(data)[2]/10 )) ,]


pca.res = PCA(t(  log10( data.big +1)  ),graph = F)


PC1 <- as.numeric(pca.res$ind$coord[,1])
PC2 <- as.numeric(pca.res$ind$coord[,2])

#sample.names = rownames(pca.res$ind$coord)

PCs <- data.frame(PC1,PC2,cell.type)

library(ggplot2)
#library(ggrepel)

P<-ggplot(PCs, aes(PC1,PC2,colour = factor(cell.type))) 
P +geom_point() + 
  xlab(paste("PC1",as.character(round(pca.res$eig[,2][1],2)),"%")) + 
  ylab(paste("PC2",as.character(round(pca.res$eig[,2][2],2)),"%")) +
  ggtitle(   "PCA on RNAseq"    ) 


########################################

library(NODES)

data.norm = pQ(data)

data.norm.big = data.norm[ (apply(data.norm,1,function(x)  sum(  x>1   ) > dim(data.norm)[2]/10 )) ,]

pca.res = PCA(t(  log10( data.norm.big +1)  ),graph = F)


PC1 <- as.numeric(pca.res$ind$coord[,1])
PC2 <- as.numeric(pca.res$ind$coord[,2])

#sample.names = rownames(pca.res$ind$coord)

PCs <- data.frame(PC1,PC2)

library(ggplot2)
#library(ggrepel)

P<-ggplot(PCs, aes(PC1,PC2,colour  = factor(sapply(strsplit(colnames(data.norm.big), split='__',fixed=TRUE),function(x) x[2])  )  ) ) 
P +geom_point() + 
  xlab(paste("PC1",as.character(round(pca.res$eig[,2][1],2)),"%")) + 
  ylab(paste("PC2",as.character(round(pca.res$eig[,2][2],2)),"%")) +
  ggtitle(   "PCA on RNAseq"    ) 




