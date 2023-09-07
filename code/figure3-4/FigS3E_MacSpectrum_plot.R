setwd("/media/ggj/ggjlab2/hezuo/gwh/myeloid/MacSpectrum/")

library(ggplot2)
library(reshape2)
library(ggpubr)
result <- read.csv("./MPI_AMDI_table (1).csv",row.names = 1)

###MPI value means inflammatory features
# result$cluster <- colsplit(result$Feature,"_",names = c("n1","n2"))$n2
# result$cluster <- gsub("cluster0","Kupffer cell",result$cluster)
# result$cluster <- gsub("cluster3","Macrophage FABP5+",result$cluster)
# result$cluster <- gsub("cluster7","Macrophage CLU+",result$cluster)
# result$cluster <- gsub("cluster9","Macrophage",result$cluster)
result$Feature <- gsub("Macropahge","Macrophage",result$Feature)
comparisions <- list(c("Kupffer cell","Macrophage"),c("Macrophage","Macrophage_CLU+"),c("Macrophage_CLU+","Macrophage_FABP5+"))
ggplot(result,aes(x=Feature,y=MPI))+geom_boxplot()+theme_classic()+stat_compare_means(comparisons = comparisions,method="wilcox.test")

ggplot(result,aes(x=Feature,y=AMDI))+geom_boxplot()+theme_classic()+stat_compare_means(comparisons = comparisions,method="wilcox.test")

#####density plot
library(RColorBrewer)
col_flg<-colorRampPalette(brewer.pal(4,"Set1"))(4)
ggplot(result,aes(x=MPI,y=AMDI,col=Feature))+stat_density2d(contour_var="ndensity")+theme_bw()+scale_color_manual(values=col_flg)
ggsave("../myeloid/MacSpectrum.pdf",width = 6,height = 4)

#M2 features
M2 <- c("CCL4","CCL13","CCL18","CCL20","CCL22","CD276","CLEC7A","CTSA","CTSB","CTSC","CTSD",
        "FN1","IL4R","IRF4","LYVE1","MMP9","MMP14","MMP19","MSR1","TGFB1","TGFB2","TGFB3",
        "TNFSF8","TNFSF12","VEGFA","VEGFB","VEGFC")
M1 <- c("CCL5","CCR7","CD40","CD86","CXCL9","CXCL10","CXCL11","IDO1","IL1A","IL1B","IL6","IRF1",
        "IRF4","KYNU")
load("../myeloid_harmony_1129.RData")
pbmc.harmony <- AddModuleScore(pbmc.harmony,features = list(M2),name="M2")
score <- data.frame(score=pbmc.harmony$M21,celltype=pbmc.harmony$celltype)
ggplot(score,aes(x=celltype,y=score))+geom_boxplot()

pbmc.harmony <- AddModuleScore(pbmc.harmony,features = list(M1),name="M1")
score2 <- data.frame(score=pbmc.harmony$M11,celltype=pbmc.harmony$celltype)
ggplot(score2,aes(x=celltype,y=score))+geom_boxplot()
