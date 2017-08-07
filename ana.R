library(data.table)
library(randomForest)
library(plyr)


data =  fread("GSE81861_Cell_Line_FPKM.csv")
data = as.data.frame(data)

data$V1 = sapply(strsplit(data$V1, split='_',fixed=TRUE),function(x) x[2])
  
#data.dedup = aggregate(. ~ V1, data=data, FUN=sum)


##############################

data.dedup = data[ !duplicated(data$V1),]

rownames(data.dedup) = data.dedup$V1

data.dedup = data.dedup[,-1]




############################################

library(NODES)

#data.dedup = pQ(data.dedup)

data.dedup = pQ(data.dedup)

##############################################

cell.type = sapply(strsplit(colnames(data.dedup), split='__',fixed=TRUE),function(x) x[2])

bio.label = factor(substr(cell.type,1,4))

data.dedup.mean = t(by( t(data.dedup)   ,bio.label,FUN = colMeans))


#aggregate(data.dedup,by = bio.label,FUN = rowMeans,simplify = T)

data.dedup.mean = data.frame(rowMeans(data.dedup[,bio.label == "A549"]),
                             rowMeans(data.dedup[,bio.label == "GM12"]),
                             rowMeans(data.dedup[,bio.label == "H1_B"]),
                             rowMeans(data.dedup[,bio.label == "H143"]),
                             rowMeans(data.dedup[,bio.label == "HCT1"]),
                             rowMeans(data.dedup[,bio.label == "IMR9"]),
                             rowMeans(data.dedup[,bio.label == "K562"])
                             )

########################################



Pt.g.normalize <- function(x) {
  
  # normalize gene expression level. 
  # this fucntion convert gene expression level into probability,
  # so can be fitted into channon entropy formula
  
  #normalized.values = vector(mode = "numeric",length = length(x))
  normalized.values = x/sum(x)
  return(normalized.values)
  
}


Hg.compute <- function(Pt.g){
  
  #take Pt.g as input
  # compute shannon entropy for each gene
  
  sum(   -(Pt.g)*log2(Pt.g)   )
  
} 

Q.g.t.compute <- function(x){
  
  # compute Q statistics for each gene, 
  # based on shannon entropy and  normalized gene expression level.
  
  #Q.g.t.vector = vector(mode = "numeric",length = length(x))
  
  Pt.g = Pt.g.normalize(x)
  Hg = Hg.compute(Pt.g)
  #print(Hg)
  #print(Pt.g)
  #print(-log2(Pt.g))
  Q.g.t.vector = Hg * (-(log2(Pt.g)))
  return(Q.g.t.vector)
}

#Q.g.t.compute(c(5,5,5,5,50))


########################
# this function compute tissue specific shannon entropy 
# for each gene and each tissue

Q.gt.matrix.compute <- function(FPKM.table){
  
  Q.gt.matrix = matrix(0, nrow = nrow(FPKM.table), ncol = ncol(FPKM.table))
  
  for(i in c(1:nrow(FPKM.table))) {
    
    Q.gt.matrix[i,] = Q.g.t.compute(   as.numeric(FPKM.table[i,] )  )
  }
  #apply(FPKM.table,1,Q.g.t.compute)
  return(Q.gt.matrix)
}

##############################
## following are analysis part.

# we use these method to identify stage specific genes using data from:
# The landscape of accessible chromatin in mammalian preimplantation embryos
###  GSE66582, this is scRNA. This data file is from this record

Q.stat.matrix = Q.gt.matrix.compute(data.dedup.mean+1)

g1 = (Q.stat.matrix[,1]<1 & data.dedup.mean[,1]>2)
g2 = (Q.stat.matrix[,2]<1 & data.dedup.mean[,2]>2)
g3 = (Q.stat.matrix[,3]<1 & data.dedup.mean[,3]>2)

g4 = (Q.stat.matrix[,4]<1 & data.dedup.mean[,4]>2)
g5 = (Q.stat.matrix[,5]<1 & data.dedup.mean[,5]>2)
g6 = (Q.stat.matrix[,6]<1 & data.dedup.mean[,6]>2)

g7 = (Q.stat.matrix[,7]<1 & data.dedup.mean[,7]>2)

g =Reduce("|",list(g1,g2,g3,g4,g5,g6,g7))

data.dedup.g = data.dedup[g,]


data.dedup.g.rank = apply(data.dedup.g, 2, function(y) rank(y) / length(y))


#pca.res = PCA(t(  log10( data.dedup.g +0.1)  ),graph = F)
pca.res = PCA(t(  ( data.dedup.g.rank)  ),graph = F)


PC1 <- as.numeric(pca.res$ind$coord[,1])
PC2 <- as.numeric(pca.res$ind$coord[,2])
PC3 <- as.numeric(pca.res$ind$coord[,3])


#sample.names = rownames(pca.res$ind$coord)

PCs <- data.frame(PC1,PC2,PC3,cell.type)

library(ggplot2)
#library(ggrepel)

P<-ggplot(PCs, aes(PC1,PC2,colour = factor(bio.label))) 
P +geom_point() + 
  xlab(paste("PC1",as.character(round(pca.res$eig[,2][1],2)),"%")) + 
  ylab(paste("PC2",as.character(round(pca.res$eig[,2][2],2)),"%")) +
  ggtitle(   "PCA on RNAseq"    ) 

library(scatterplot3d)
#attach(mtcars)
scatterplot3d(PC1,PC2,PC3, main="3D Scatterplot",color = rainbow(7)[factor(bio.label)])


library(rgl)

plot3d(PC1,PC2,PC3, col=rainbow(7)[factor(bio.label)], size=3)

###########################################################







