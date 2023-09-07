setwd("/media/ggj/ggjlab2/hezuo/gwh/")
#dir.create("./cellchat")

library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(patchwork)

load("./SCT.rdata")
anno <- read.csv("./anno_0630.csv",row.names = 1)
anno$cellID <- rownames(anno)
cluster <- read.csv("./anno(3).csv",row.names = 1)
cluster$cellID <- rownames(cluster)
anno.use <- merge(anno,cluster,by='cellID')
colnames(anno.use)[7] <- "celltype"
table(anno.use$Type,anno.use$celltype)
#anno.use$Type <- gsub("ICC","CA",anno.use$Type)
table(anno.use$Type)
anno.use$celltype <- gsub("Epithelial/Cancer cell","Cancer cell",anno.use$celltype)
anno.use$celltype <- gsub("Macrophage|Dendritic cell|Kuppfer cell|Monocyte","Monocyte",anno.use$celltype)
anno.use$celltype <- gsub("NKT cell","T cell",anno.use$celltype)
anno.use2 <- anno.use[-which(anno.use$celltype=="Erythroid"),]
table(anno.use2$Type,anno.use2$celltype)
anno.use <- anno.use[anno.use$celltype%in%c("Cancer cell","Epithelial","Fibroblast","Hepatocyte","T cell","Plasma cell",
                                            "Endothelial"),]
####myeloid
anno.myeloid <- read.csv("./myeloid/myeloid_cell_anno.csv",row.names = 1) 
colnames(anno.myeloid) <- 'celltype'
anno.myeloid$cellID <- rownames(anno.myeloid)
###T cell
# anno.T <- read.csv("./Tcell/T_cell_anno.csv",row.names = 1)
# colnames(anno.T) <- 'celltype'
# anno.T$cellID <- rownames(anno.T)
# ###B cell
# anno.B <- read.csv("./Bcell/B_cell_anno.csv",row.names = 1)
# colnames(anno.B) <- 'celltype'
# anno.B$cellID <- rownames(anno.B)
# ###Endothelial
# anno.endo <- read.csv("./endothelial/endothelial_anno.csv",row.names = 1)
# colnames(anno.endo) <- 'celltype'
# anno.endo$cellID <- rownames(anno.endo)

anno.use <- anno.use[,c(7,1)]

anno <- rbind(anno.use,anno.myeloid)
#anno <- rbind(anno,anno.T)
#anno <- rbind(anno,anno.B)
#anno <- rbind(anno,anno.endo)
length(unique(anno$cellID))

table(anno$celltype)
#cell <- anno$cellID[duplicated(anno$cellID)]
#anno <- anno[-which(anno$cellID%in%cell),]
pbmc.use <- pbmc[,anno$cellID]
dim(pbmc.use)

colnames(pbmc.use)[1:3];anno$cellID[1:3]

dim(pbmc)
pbmc.use$celltype <- anno$celltype
table(pbmc.use$Type)
Idents(pbmc.use) <- pbmc.use$Type
unique(pbmc.use$celltype)
pbmc.use$celltype2 <- pbmc.use$celltype
pbmc.use$celltype2 <- gsub("Cancer cell|Hepatocyte","Epithelial",pbmc.use$celltype2)
unique(pbmc.use$celltype2)
table(pbmc.use$celltype2,pbmc.use$Type)
pbmc.use$type2 <- pbmc.use$Type 
pbmc.use$type2 <- gsub("ICC","CA",pbmc.use$type2)
Idents(pbmc.use) <- pbmc.use$type2

pbmc.use$type3 <- paste0(pbmc.use$type2,"_",pbmc.use$celltype2)
type3 <- as.data.frame(table(pbmc.use$type3))
head(type3)
type3 <- type3[type3$Freq>=10,]

Idents(pbmc.use) <- pbmc.use$type3
pbmc.use <- subset(pbmc.use,idents=c(type3$Var1))

####cellchat
Idents(pbmc.use) <- pbmc.use$type2
CA.seurat <- subset(pbmc.use,idents=c("CA"))
#ICC.seurat <- subset(pbmc.use,idents=c("ICC"))
Adj.seurat <- subset(pbmc.use,idents=c("AD"))
Cir.seurat <- subset(pbmc.use,idents=c("Fibrosis"))
Normal.seurat <- subset(pbmc.use,idents=c("Normal"))

table(pbmc.use$celltype2,pbmc.use$type2)

CA <- createCellChat(CA.seurat@assays$RNA@data,meta = CA.seurat@meta.data,group.by = "celltype2")
#ICC <- createCellChat(ICC.seurat@assays$RNA@data,meta = ICC.seurat@meta.data,group.by = "celltype2")
Adj <- createCellChat(Adj.seurat@assays$RNA@data,meta = Adj.seurat@meta.data,group.by = "celltype2")
Cir <- createCellChat(Cir.seurat@assays$RNA@data,meta = Cir.seurat@meta.data,group.by = "celltype2")
Normal <- createCellChat(Normal.seurat@assays$RNA@data,meta = Normal.seurat@meta.data,group.by = "celltype2")

dir.create("./cellchat/myeloid_new")
setwd("./cellchat/myeloid_new/")
type <- c("CA","Adj","Cir","Normal")
library(future)
plan('multisession',workers=4)
for(i in 1:length(type)){
  cellchat <- get(type[i])  
  cellchat@DB <- CellChatDB.human  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat,raw.use = T,population.size = TRUE) 
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat,slot.name = "netP")
  
  savename <- paste0(type[i],".RData")
  save(cellchat,file=savename)
}

#######merge
rm(list=ls())
gc()

library(reshape2)
file <- list.files(pattern="RData")
for(i in 1:length(file)){
  load(file[i])
  file1 <- colsplit(file[i],"\\.",names=c("n1","n2"))$n1
  assign(file1,cellchat)
}

coo.list <- list(CA=CA,Adj=Adj,Cirrhosis=Cir,Normal=Normal)
coo.list <- list(CA=CA,Normal=Normal)
cellchat <- mergeCellChat(coo.list,add.names = names(coo.list),cell.prefix = T)

p1 <- compareInteractions(cellchat,show.legend = F,group = c(1,2,3,4),measure = "count");p1
p2 <- compareInteractions(cellchat,show.legend = F,group = c(1,2,3,4),measure = "weight");p2
p <- p1+p2
p

par(mfrow=c(1,1))
netVisual_diffInteraction(cellchat,comparison = c(1,4),weight.scale = T,measure = "weight")
netVisual_heatmap(cellchat,comparison = c(4,1))

par(mfrow=c(2,2))
weight.max <- getMaxWeight(coo.list,attribute = c("idents","counts"))
for(i in 1:length(coo.list)){
  netVisual_circle(coo.list[[i]]@net$count,weight.scale=T,label.edge = F,
                   edge.weight.max = weight.max[2],edge.width.max = 12,
                   title.name = paste0("Number of interaction - ",names(coo.list)[i]))
}
edit(rankNet)
rankNet(cellchat,mode = "comparison",stacked=T,comparison = c(1,2),do.stat=TRUE)
netVisual_bubble(cellchat,comparison = c(3,4), sources.use = 7 ,targets.use = c(1:17))
netVisual_diffInteraction(cellchat)
###HCC
type <- c("CA","Adj","Cir","Normal")
for (i in 1:length(type)) {
  mat <- coo.list[[i]]@net$count
  groupSize <- as.numeric(table(coo.list[[i]]@idents))
  savename <- paste0(type[i],"_celltype.pdf")
  pdf(savename,width = 30,height = 30)
  par(mfrow=c(6,6))
  for(i in 1:nrow(mat)){
    mat2 <- matrix(0,nrow=nrow(mat),ncol=ncol(mat),dimnames = dimnames(mat))
    mat2[i,] <- mat[i,]
    netVisual_circle(mat2,vertex.weight = groupSize,weight.scale=T,
                     edge.weight.max = max(mat),
                     title.name = rownames(mat)[i])
  }
  dev.off()
}


mat <- coo.list[[1]]@net$count
groupSize <- as.numeric(table(coo.list[[1]]@idents))
pdf('./HCC_celltype.pdf',width = 30,height = 30)
par(mfrow=c(6,6))
for(i in 1:nrow(mat)){
  mat2 <- matrix(0,nrow=nrow(mat),ncol=ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale=T,
                   edge.weight.max = max(mat),
                   title.name = rownames(mat)[i])
}
dev.off()
