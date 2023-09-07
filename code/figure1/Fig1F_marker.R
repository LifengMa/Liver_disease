#Figure 1F
setwd("/media/ggj/ggjlab2/hezuo/gwh/")

library(reshape2)
library(Seurat)
library(ggplot2)
library(dplyr)

load("./SCT.rdata")
anno <- read.csv("./anno_datu_1205.csv",row.names = 1)
marker <- read.csv("./markers_normalize.csv",row.names = 1)
marker <- marker[marker$p_val_adj<0.05,]
top10 <- marker %>% group_by(cluster) %>% top_n(15, avg_log2FC)
top10 <- top10[,c(2,6,7)]
top10 <- top10[-(grep("MT-",top10$gene)),]

anno.use <- anno[,c(1,6,8)]
anno.use <- anno.use[!duplicated(anno.use),]
colnames(anno.use)[1] <- "cluster"

anno.use <- merge(anno.use,top10,by="cluster")
anno.use <- anno.use[order(anno.use$celltype2),]
anno.use <- anno.use[,c(3:5)]
anno.use <- anno.use[,c(1,3,2)]
anno.use2 <- dcast(anno.use,celltype3~gene)
rownames(anno.use2) <- anno.use2$celltype3
anno.use2 <- anno.use2[,-1]

gene.use <- unique(anno.use$gene)
cluster.order <- unique(anno.use$celltype3)

pbmc$celltype <- anno$celltype3
Idents(pbmc) <- pbmc$celltype
DefaultAssay(pbmc) <- "RNA"
pbmc <- ScaleData(pbmc,features = gene.use)
avg <- AverageExpression(pbmc,slot="scale.data")
avg1 <- as.data.frame(avg$RNA)
#avg <- log(avg+1)

avg1 <- as.data.frame(t(avg1))
avg1 <- avg1[cluster.order,gene.use]
max(avg1);min(avg1)
avg1[avg1>(2)] <- (2)
avg1[avg1<(-1)] <- (-1)

library(RColorBrewer)
anno.use3 <- anno[,c(1,6,8)]
anno.use3 <- anno.use3[!duplicated(anno.use3),]
anno.use3 <- anno.use3[,c(2,3)] 
rownames(anno.use3) <- anno.use3$celltype3
anno.use3 <- anno.use3['celltype2']
colnames(anno.use3) <- "celltype"
anno.col <- list(
  celltype=c(
    `Dendritic cell`='#1f77b4',
    Endothelial='#ff7f0e',
    Epithelial='#279e68', 
    Erythroid='#d62728',
    Fibroblast='#aa40fc',
    Kupffer='#8c564b',
    Macrophage='#e377c2',
    `Mast cell`='#b5bd61',
    Monocyte='#17becf',
    Neutrophil='#aec7e8',
    `Plasma cell`='#ffbb78',
    `T cell`='#98df8a'
  ))

bk <- c(seq(-1,-0.1,by=0.01),seq(0,2,by=0.02))
pheatmap::pheatmap(t(avg1),cluster_rows = F,cluster_cols = F,border_color = "NA",show_rownames = F,
                   color = c(colorRampPalette(colors = c("#FFFACD","white"))(length(bk)/2),colorRampPalette(colors = c("white","darkgreen"))(length(bk)/2)),
                   breaks = bk)
pdf("./figures/marker.pdf",height = 10,width = 10)
pheatmap::pheatmap(t(avg1),cluster_rows = F,cluster_cols = F,border_color = "NA",show_rownames = F,
                   color = colorRampPalette(brewer.pal(n = 7, name =
                                                             "YlGnBu"))(100))
dev.off()
