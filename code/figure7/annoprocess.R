setwd("/media/ggj/ggjlab2/hezuo/gwh/")

library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(patchwork)

load("./SCT.rdata")
all <- pbmc
anno <- read.csv("./anno_0630.csv",row.names = 1)
anno$cellID <- rownames(anno)
cluster <- read.csv("./anno(3).csv",row.names = 1)
cluster$cellID <- rownames(cluster)
anno.use <- merge(anno,cluster,by='cellID')
colnames(anno.use)[7] <- "celltype"
table(anno.use$Type,anno.use$celltype)
anno.use$Type <- gsub("ICC","CA",anno.use$Type)
anno.use$celltype <- gsub("Epithelial/Cancer cell","Cancer cell",anno.use$celltype)
anno.use$celltype <- gsub("Macrophage|Dendritic cell|Kuppfer cell|Monocyte","Monocyte",anno.use$celltype)
anno.use$celltype <- gsub("NKT cell","T cell",anno.use$celltype)
anno.use2 <- anno.use[-which(anno.use$celltype=="Erythroid"),]
table(anno.use2$Type,anno.use2$celltype)
anno.use <- anno.use[anno.use$celltype%in%c("Cancer cell","Epithelial","Fibroblast","Hepatocyte"),]
####myeloid
load("./myeloid/myeloid_harmony_1129.RData")
DimPlot(pbmc.harmony,label=T)
anno.myeloid <- read.csv("./myeloid/myeloid_cell_anno.csv",row.names = 1) 
colnames(anno.myeloid) <- 'celltype'
anno.myeloid$cellID <- rownames(anno.myeloid)
unique(anno.myeloid$celltype)
anno.myeloid$celltype <- gsub("Plasmacytoid Dendritic Cell_GZMB\\+","pDC",anno.myeloid$celltype)
anno.myeloid$celltype <- gsub("Plasmacytoid Dendritic Cell","pDC",anno.myeloid$celltype)
anno.myeloid$celltype <- gsub("Conventional Dendritic Cell","cDC",anno.myeloid$celltype)
anno.myeloid$celltype <- gsub("Neutrophil_DEFA3\\+|Neutrophil_CAMP\\+","Neutrophil",anno.myeloid$celltype)
###T cell
load("./Tcell/Tcell_harmony.RData")
DimPlot(pbmc.harmony,label=T)
library(openxlsx)
library(RColorBrewer)
anno <- read.xlsx("./Tcell/Tcell_marker_harmony.xlsx")
anno <- anno[,c(7,9)]
anno <- anno[!duplicated(anno),]
anno <- na.omit(anno)
anno
anno$celltype <- gsub("Helper T cell_RP\\+","Naive T cell",anno$celltype)
anno$celltype <- gsub("Helper T cell_TNF high","Activated T cell",anno$celltype)
anno$celltype <- gsub("T cell_RP\\+","Naive T cell",anno$celltype)
anno$celltype <- gsub("NK cell_XCL1\\+","NK cell",anno$celltype)
anno$celltype <- gsub("NK cell_GNLY\\+","NK cell",anno$celltype)
anno
new.cluster.ids<-as.character(anno$celltype)
names(new.cluster.ids) <- levels(pbmc.harmony)
pbmc.harmony <- RenameIdents(pbmc.harmony, new.cluster.ids)
pbmc.harmony@meta.data$celltype<-Idents(pbmc.harmony)
col_flg<-colorRampPalette(brewer.pal(8,"Set1"))(length(levels(as.factor(pbmc.harmony$celltype))))
DimPlot(object = pbmc.harmony, reduction = "umap",label  = T ,cols = col_flg  ,pt.size = 0.5)
anno <- as.data.frame(pbmc.harmony$celltype)
anno$cellID <- rownames(anno)
write.csv(anno,file="./Tcell/T_cell_anno.csv")
anno.T <- read.csv("./Tcell/T_cell_anno.csv",row.names = 1)
colnames(anno.T)[1] <- 'celltype'
#anno.T$cellID <- rownames(anno.T)
unique(anno.T$celltype)
###B cell
load("./Bcell/Bcell_harmony.RData")
DimPlot(pbmc.harmony,label=T)
library(openxlsx)
library(RColorBrewer)
anno <- read.xlsx("./Bcell/Bcell_marker_harmony.xlsx")
anno <- anno[,c(7,9)]
anno <- anno[!duplicated(anno),]
anno <- na.omit(anno)
anno

anno.B <- read.csv("./Bcell/B_cell_anno.csv",row.names = 1)
colnames(anno.B) <- 'celltype'
anno.B$cellID <- rownames(anno.B)
unique(anno.B$celltype)
anno.B$celltype <- colsplit(anno.B$celltype,"_",names = c("n1","n2"))$n1
###Endothelial

anno.endo <- read.csv("./endothelial/endothelial_anno.csv",row.names = 1)
colnames(anno.endo) <- 'celltype'
anno.endo$cellID <- rownames(anno.endo)
unique(anno.endo$celltype)
library(reshape2)
anno.endo$celltype <- gsub("vascular_J547endo","vascular_endo",anno.endo$celltype)
anno.endo$celltype <- colsplit(anno.endo$celltype,"_",names = c("n1","n2"))$n1
anno.endo$celltype <- gsub("vascular","VE",anno.endo$celltype)
anno.endo$celltype <- gsub("Endo","CAE",anno.endo$celltype)

anno.use <- anno.use[,c(7,1)]

anno <- rbind(anno.use,anno.myeloid)
anno <- rbind(anno,anno.T)
anno <- rbind(anno,anno.B)
anno <- rbind(anno,anno.endo)
length(unique(anno$cellID))

table(anno$celltype)
cell <- anno$cellID[duplicated(anno$cellID)]
anno <- anno[-which(anno$cellID%in%cell),]
anno <- anno[-which(anno$celltype=="Proliferating cell"),]
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
table(pbmc.use$celltype2,pbmc.use$type2)

saveRDS(pbmc.use,file="./nichenet/seuratobject.rds")
