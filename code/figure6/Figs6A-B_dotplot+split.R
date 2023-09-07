setwd("/media/ggj/ggjlab2/hezuo/gwh/endothelial/")

dyn.load("/home/ggj/Documents/glpk-5.0/src/.libs/libglpk.so.40") 
Sys.getenv("LD_LIBRARY_PATH")
Sys.setenv(LD_LIBRARY_PATH=paste0(Sys.getenv("LD_LIBRARY_PATH"), ":", "/home/ggj/Documents/glpk-5.0/src/.libs/"))

library(Seurat)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

load("./endothelial_harmony.RData")
DimPlot(pbmc.harmony,label=T)

#Dotplot
features <- c("CLEC1B","CLEC4G","FCGR2B","CLDN5",
              "PLVAP","RGCC","VWF","PECAM1","RAMP2",
              "MGP","CPE")
DotPlot(pbmc.harmony,features = features,cols = c("blue", "red"),dot.scale = 8)+theme(axis.text.x = element_text(angle=45,hjust=1))+
  ylab("Cluster")+xlab("Genes")
dir.create("./Endothelial")
ggsave("./Endothelial/Dotplot.pdf",width = 12,height = 6)




# split -------------------------------------------------------------------

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
ggsave("./Endothelial/split_type.pdf",width = 15,height = 6)

