setwd("/media/ggj/ggjlab2/hezuo/gwh/Tcell/")

dyn.load("/home/ggj/Documents/glpk-5.0/src/.libs/libglpk.so.40") 
Sys.getenv("LD_LIBRARY_PATH")
Sys.setenv(LD_LIBRARY_PATH=paste0(Sys.getenv("LD_LIBRARY_PATH"), ":", "/home/ggj/Documents/glpk-5.0/src/.libs/"))

library(Seurat)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

load("./Tcell_harmony.RData")
DimPlot(pbmc.harmony,label=T)

#Dotplot (Fig5B)
features <- c("CD3D","LTB",
              "TNF","JUN","JUNB",
              "CCL3","XCL1","XCL2","KLRF1",
              "CD8A","GZMK","CCL5",
              "FGFBP2","NKG7","GNLY",
              "CTLA4",
              "TNFRSF4","TIGIT")
DotPlot(pbmc.harmony,features = features,cols = c("blue", "red"),dot.scale = 8)+theme(axis.text.x = element_text(angle=45,hjust=1))+
  ylab("Cluster")+xlab("Genes")
dir.create("./Tcell")
ggsave("./Tcell/Dotplot.pdf",width = 12,height = 6)




# split (FigS5A) -------------------------------------------------------------------

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
ggsave("./Tcell/split_type.pdf",width = 20,height = 6)
