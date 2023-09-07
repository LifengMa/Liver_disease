setwd("/media/ggj/ggjlab2/hezuo/gwh/myeloid/")

dyn.load("/home/ggj/Documents/glpk-5.0/src/.libs/libglpk.so.40") 
Sys.getenv("LD_LIBRARY_PATH")
Sys.setenv(LD_LIBRARY_PATH=paste0(Sys.getenv("LD_LIBRARY_PATH"), ":", "/home/ggj/Documents/glpk-5.0/src/.libs/"))

library(Seurat)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(gridExtra)
load("./myeloid_harmony_1129.RData")
DimPlot(pbmc.harmony,label=T)

# cell fraction (Fig 3C)-----------------------------------------------------------
##ggplot2
#pbmc.harmony$type <- factor(pbmc.harmony$type,levels=c("CA","Fibrosis","AD","Normal"))
unique(pbmc.harmony$Batch2)
Idents(pbmc.harmony) <- pbmc.harmony$seurat_clusters
final<-data.frame(Idents(pbmc.harmony),pbmc.harmony$Batch2)
colnames(final)<-c("cluster","state")
#final$state <- gsub("S1|S2","",final$state)
unique(final$state)
final$state <- as.factor(final$state)
#table(pbmc.harmony$type)

cluster1<-final[final$cluster==0,]
table(cluster1$state)
clustername<-data.frame(table(cluster1$cluster))
statename<-data.frame(table(cluster1$state))
clusternum<-nrow(data.frame(table(final$cluster)))
state<-data.frame(matrix(nrow = nrow(statename),ncol = 0),row.names = statename$Var1)
for (i in 1: clusternum ) {
  temp<-final[final$cluster==(i-1),]
  state0<-data.frame(table(temp$state))
  state<-cbind(state,state0[,2])
}
colnames(state) <- paste0("cluster",c(0:(length(colnames(state))-1)))

#state$num <- apply(state,1,sum)
state <- apply(state,1,function(x){x/sum(x)})
state <- as.data.frame(t(state))
state$patient <- rownames(state)

state.use <- melt(state)
colnames(state.use) <- c("donor","cluster","fraction")
state.use$type <- state.use$donor
state.use$type <- gsub("[0-9]","",state.use$type)
state.use$type <- gsub("ICC","CA",state.use$type)
state.use$type <- gsub("S","",state.use$type)
state.use$type <- gsub("Liver","Normal",state.use$type)
state.use$type <- gsub("Fib","Cirhosis",state.use$type)
state.use$type <- gsub("CAAD","ADJ",state.use$type)
unique(state.use$type)
state.use$fraction <- state.use$fraction*100
source("/media/ggj/ggjlab2/hezuo/gwh/errorbar.R")
colnames(state.use)
head(state.use)
#df2 <- data_summary(state.use,varname = "fraction",groupnames = c("type","cluster"))

df2 <- state.use
library(ggplot2)
df2$type <- factor(df2$type,levels=c("CA","Cirhosis","ADJ","Normal"))
library(RColorBrewer)
col_flg<-colorRampPalette(brewer.pal(8,"Paired"))(length(levels(as.factor(df2$cluster))))

ggplot(df2,aes(x=type,y=fraction,fill=cluster))+geom_bar(stat = "summary",fun="mean",position = position_dodge())+
  stat_summary(fun.data = 'mean_se',geom="errorbar",width=0.15,position = position_dodge(0.9))+
  #stat_compare_means(label="p.signif",method="wilcox.test",ref.group=".all.")+
  facet_grid(.~cluster)+ylab("Celluar fraction(%)")+scale_fill_manual(values=col_flg)+
  theme_bw()+theme(axis.text.x = element_text(angle=45,hjust=1))+NoLegend()


###Normal distribution test
library(dplyr)
shapiro.test(state.use$fraction)

#homogeneity test for variance
bartlett.test(state.use$fraction,state.use$type)

library(ggpubr)
library(gridExtra)
fig <- list()

for(i in 1:length(unique(df2$cluster))){
  temp <- df2[df2$cluster==unique(df2$cluster)[i],]
  #temp2 <- state.use[state.use$cluster==unique(df2$cluster)[i],]
  if(i==1){
    p <-  ggplot(temp,aes(x=type,y=fraction,fill=cluster))+geom_bar(stat = "summary",fun="mean",position = position_dodge())+
      stat_summary(fun.data = 'mean_se',geom="errorbar",width=0.15,position = position_dodge(0.9))+
      stat_compare_means(label="p.signif",method="wilcox.test",ref.group = 'Normal',hide.ns = TRUE)+
      #geom_hline(yintercept = mean(temp$fraction),linetype=2)+
      ylab("Celluar fraction(%)")+scale_fill_manual(values=col_flg[i])+
      ggtitle(unique(df2$cluster)[i])+xlab("")+
      theme_bw()+theme(axis.text.x = element_text(angle=45,hjust=1))+NoLegend()
  }
  else{
    p <-  ggplot(temp,aes(x=type,y=fraction,fill=cluster))+geom_bar(stat = "summary",fun="mean",position = position_dodge())+
      stat_summary(fun.data = 'mean_se',geom="errorbar",width=0.15,position = position_dodge(0.9))+
      stat_compare_means(label="p.signif",method="wilcox.test",ref.group="Normal",hide.ns = TRUE)+
      #geom_hline(yintercept = mean(temp$fraction),linetype=2)+
      ylab("")+scale_fill_manual(values=col_flg[i])+
      ggtitle(unique(df2$cluster)[i])+xlab("")+
      theme_bw()+theme(axis.text.x = element_text(angle=45,hjust=1))+NoLegend()
  }
  
  fig[[i]] <- p
}

fig[['nrow']] <- 1
fig[['ncol']] <- 15

pdf("./myeloid/cellfraction.pdf",width = 28,height = 4)
do.call('grid.arrange', fig)
dev.off()
#facet_grid(.~cluster)
#compare_means(data=temp2, fraction~type,label="p.signif",method="wilcox.test",ref.group = ".all.")



# URD (Fig 4A)---------------------------------------------------------------------

library(URD)
Idents(pbmc.harmony) <- pbmc.harmony$seurat_clusters
DimPlot(pbmc.harmony,label=T)
monocyte <- subset(pbmc.harmony,idents = c(1:9,11,12))
DimPlot(monocyte,label=T)
data <- createURD(monocyte@assays$RNA@counts,meta=monocyte@meta.data,min.cells = 3,min.counts = 3)
data <- findVariableGenes(data,set.object.var.genes = T,diffCV.cutoff = 0.3,do.plot=T)
data <- calcPCA(data,mp.factor = 2)
pcSDPlot(data)
set.seed(19)
data <- calcTsne(data)
plotDim(data,"celltype")
plotDim(data,"FABP5")
plotDim(data,"SPP1")
plotDim(data,"CD14")
plotDim(data,"C1QB")
data <- calcDM(data,knn=100,sigma=16)
plotDim(data,"seurat_clusters",transitions.plot=10000)
plotDim(data,"celltype",transitions.plot=1)
data@group.ids$cluster <- as.character(data@meta[rownames(data@group.ids),"seurat_clusters"])
root.cells <- cellsInCluster(data,clustering = "cluster",cluster=c('4'))
data.floods <- floodPseudotime(data , root.cells = root.cells,n=50,minimum.cells.flooded = 2,verbose=T)
data <- floodPseudotimeProcess(data,data.floods,floods.name = "pseudotime")
pseudotimePlotStabilityOverall(data)
plotDim(data,"pseudotime")
plotDists(data,"pseudotime","seurat_clusters")
plotDists(data,"pseudotime","celltype")
###find tips
data.tips <- urdSubset(data,cells.keep = cellsInCluster(data,"cluster",cluster = c("3","6")))
data.tips <- findVariableGenes(data.tips,set.object.var.genes = T,diffCV.cutoff = 0.3,do.plot=T)
data.tips <- calcPCA(data.tips,mp.factor = 1.5)
pcSDPlot(data.tips)
set.seed(20)
data.tips <- calcTsne(data.tips)
data.tips <- graphClustering(data.tips,num.nn=500,do.jaccard = T,method="Louvain")
plotDim(data.tips,"Louvain-500",point.size=3)
plotDim(data.tips,"CLEC9A")
plotDim(data.tips,"MARCO")
plotDim(data.tips,"FABP5")
plotDim(data.tips,"VEGFA")
plotDim(data.tips,"SPP1")
#Biased random walks
anno <- data@meta
head(anno)
anno <- anno[rownames(data.tips@group.ids),]
unique(anno$seurat_clusters)
anno$seurat_clusters <- as.character(anno$seurat_clusters)
anno$use <- anno$seurat_clusters
#anno$use <- gsub("7","1",anno$use)
data@group.ids[rownames(data.tips@group.ids),"tip.clusters"] <- data.tips@group.ids$'Louvain-500'
data@group.ids[rownames(data.tips@group.ids),"tip.clusters"] <- anno$use

data.ptlogistic <- pseudotimeDetermineLogistic(data,"pseudotime",optimal.cells.forward = 20,max.cells.back = 40,do.plot=T)


data.biased.tm <- as.matrix(pseudotimeWeightTransitionMatrix(data,"pseudotime",logistic.params = data.ptlogistic))
data.walks <- simulateRandomWalksFromTips(data,tip.group.id = "tip.clusters",root.cells = root.cells,transition.matrix = data.biased.tm,n.per.tip = 25000,
                                          root.visits = 1,max.steps = 5000,verbose = F)
data <- processRandomWalksFromTips(data,data.walks,verbose=F)

plotDim(data,"tip.clusters")

#Build tree
data.tree <- loadTipCells(data,"tip.clusters")
data.tree <- buildTree(data.tree,pseudotime="pseudotime",tips.use = c(2,4),divergence.method = "preference")

data.tree <- nameSegments(data.tree,segments = c("4","2"),segment.names = c("Macrophage_FABP5+","Dendritic cell"))
library(RColorBrewer)
col_flg<-colorRampPalette(brewer.pal(8,"Paired"))(14)
col_flg <- col_flg[c(0:9,11,12)]
pdf("./myeloid/URD.pdf",width = 8,height = 6)
plotTree(data.tree,"seurat_clusters",discrete.colors  = col_flg)
# plotTree(data.tree,"celltype")
dev.off()
plotTree(data.tree,"celltype")
gene <- c("FABP5","SPP1","TREM2","VEGFA","CD1C","FCER1A","WDFY4","CLEC9A")
plot2 <- list()
for(i in 1:length(gene)){
  if(i %in%c(9:10)){
    p <- plotTree(data.tree,gene[i],title=gene[i])
  }
  else{
    p <- plotTree(data.tree,gene[i],title=gene[i])+theme(axis.text.x =element_blank())
  }
  plot2[[i]] <- p
}
plot2[['nrow']] <- 2
plot2[['ncol']] <- 4

pdf("./myeloid/URD_gene.pdf",width = 10,height = 4)
do.call('grid.arrange', plot2)
dev.off()

plotTree(data.tree,"FABP5",title="FABP5")
plotTree(data.tree,"VEGFA",title="VEGFA")
plotTree(data.tree,"TREM2",title="TREM2")
plotTree(data.tree,"CD1C",title="CD1C")
plotTree(data.tree,"SPP1",title="SPP1")
plotTree(data.tree,"CLEC9A",title="CLEC9A")
plotTree(data.tree,"WDFY4",title="WDFY4")
plotTree(data.tree,"CD5L",title="CD5L")
plotTree(data.tree,"SLC40A1",title="SLC40A1")
plotTree(data.tree,"CD14",title="C1QA")

rm(pbmc.harmony)
save.image("./URD_1207.RData")


# URD_no kupffer cell ---------------------------------------------------------------------

library(URD)
data <- createURD(pbmc.harmony@assays$RNA@counts,meta=pbmc.harmony@meta.data,min.cells = 3,min.counts = 3)
data <- findVariableGenes(data,set.object.var.genes = T,diffCV.cutoff = 0.3,do.plot=T)
data <- calcPCA(data,mp.factor = 2)
pcSDPlot(data)
set.seed(19)
data <- calcTsne(data)
plotDim(data,"seurat_clusters")
plotDim(data,"FABP5")
plotDim(data,"SPP1")
plotDim(data,"CD14")

data <- calcDM(data,knn=100,sigma=16)
plotDim(data,"seurat_clusters",transitions.plot=10000)

data@group.ids$cluster <- as.character(data@meta[rownames(data@group.ids),"seurat_clusters"])
root.cells <- cellsInCluster(data,clustering = "cluster",cluster=c('3','4'))
data.floods <- floodPseudotime(data , root.cells = root.cells,n=50,minimum.cells.flooded = 2,verbose=T)
data <- floodPseudotimeProcess(data,data.floods,floods.name = "pseudotime")
pseudotimePlotStabilityOverall(data)
plotDim(data,"pseudotime")
plotDists(data,"pseudotime","seurat_clusters")

###find tips
data.tips <- urdSubset(data,cells.keep = cellsInCluster(data,"cluster",cluster = c("0","6","1")))
data.tips <- findVariableGenes(data.tips,set.object.var.genes = T,diffCV.cutoff = 0.3,do.plot=T)
data.tips <- calcPCA(data.tips,mp.factor = 1.5)
pcSDPlot(data.tips)
set.seed(20)
data.tips <- calcTsne(data.tips)
data.tips <- graphClustering(data.tips,num.nn=2000,do.jaccard = T,method="Louvain")
plotDim(data.tips,"Louvain-500",point.size=3)
plotDim(data.tips,"CLEC9A")
plotDim(data.tips,"MARCO")
plotDim(data.tips,"FABP5")
plotDim(data.tips,"VEGFA")
plotDim(data.tips,"SPP1")
#Biased random walks
anno <- data@meta
head(anno)
anno <- anno[rownames(data.tips@group.ids),]
unique(anno$seurat_clusters)
anno$seurat_clusters <- as.character(anno$seurat_clusters)
anno$use <- anno$seurat_clusters
#anno$use <- gsub("7","1",anno$use)
#data@group.ids[rownames(data.tips@group.ids),"tip.clusters"] <- data.tips@group.ids$'Louvain-2000'
data@group.ids[rownames(data.tips@group.ids),"tip.clusters"] <- anno$use
data.ptlogistic <- pseudotimeDetermineLogistic(data,"pseudotime",optimal.cells.forward = 20,max.cells.back = 40,do.plot=T)


data.biased.tm <- as.matrix(pseudotimeWeightTransitionMatrix(data,"pseudotime",logistic.params = data.ptlogistic))
data.walks <- simulateRandomWalksFromTips(data,tip.group.id = "tip.clusters",root.cells = root.cells,transition.matrix = data.biased.tm,n.per.tip = 25000,
                                          root.visits = 1,max.steps = 5000,verbose = F)
data <- processRandomWalksFromTips(data,data.walks,verbose=F)

plotDim(data,"tip.clusters")

#Build tree
data.tree <- loadTipCells(data,"tip.clusters")
data.tree <- buildTree(data.tree,pseudotime="pseudotime",tips.use = c(0,1,6),divergence.method = "preference")

data.tree <- nameSegments(data.tree,segments = c("0","1","6"),segment.names = c("Kupffer cell","Macrophage_FABP5+","Dendritic cell"))
pdf("./figure3/URD.pdf",width = 8,height = 6)
plotTree(data.tree,"seurat_clusters",discrete.colors  = col_flg)
dev.off()

gene <- c("FABP5","SPP1","TREM2","MARCO","CD5L","SLC40A1","CD1C","WDFY4","CLEC9A")
plot2 <- list()
for(i in 1:length(gene)){
  if(i %in%c(7:9)){
    p <- plotTree(data.tree,gene[i],title=gene[i])
  }
  else{
    p <- plotTree(data.tree,gene[i],title=gene[i])+theme(axis.text.x =element_blank())
  }
  plot2[[i]] <- p
}
plot2[['nrow']] <- 3
plot2[['ncol']] <- 3

pdf("./figure3/URD_gene.pdf",width = 15,height = 10)
do.call('grid.arrange', plot2)
dev.off()

plotTree(data.tree,"FABP5",title="FABP5")
plotTree(data.tree,"VEGFA",title="VEGFA")
plotTree(data.tree,"TREM2",title="TREM2")
plotTree(data.tree,"CD1C",title="CD1C")
plotTree(data.tree,"SPP1",title="SPP1")
plotTree(data.tree,"CLEC9A",title="CLEC9A")
plotTree(data.tree,"WDFY4",title="WDFY4")
plotTree(data.tree,"CD5L",title="CD5L")
plotTree(data.tree,"SLC40A1",title="SLC40A1")
plotTree(data.tree,"CD14",title="CD14")

rm(pbmc.harmony)
save.image("./URD.RData")


# split_type (Fig S3A)--------------------------------------------------------------
pbmc.harmony$type2 <- pbmc.harmony$Batch2
pbmc.harmony$type2 <- gsub("[0-9]","",pbmc.harmony$type2)
pbmc.harmony$type2 <- gsub("S","",pbmc.harmony$type2)
pbmc.harmony$type2 <- gsub("Liver","Normal",pbmc.harmony$type2)
pbmc.harmony$type2 <- gsub("Fib","Cirrhosis",pbmc.harmony$type2)
pbmc.harmony$type2 <- gsub("CA","HCC",pbmc.harmony$type2)
pbmc.harmony$type2 <- gsub("HCCAD","ADJ",pbmc.harmony$type2)
pbmc.harmony$type2 <- gsub("ICHCCD","ADJ",pbmc.harmony$type2)
unique(pbmc.harmony$type2)
pbmc.harmony$type2 <- factor(pbmc.harmony$type2,levels=c("HCC","ICC","Cirrhosis","ADJ","Normal"))

library(RColorBrewer)
col_flg<-colorRampPalette(brewer.pal(8,"Paired"))(length(levels(as.factor(pbmc.harmony$seurat_clusters))))

DimPlot(pbmc.harmony,split.by = "type2",cols = col_flg,pt.size = 0.1)
ggsave("./myeloid/split_type.pdf",width = 20,height = 6)
