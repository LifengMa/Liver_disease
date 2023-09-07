#Fig 4I
setwd("/media/ggj/ggjlab2/hezuo/gwh/myeloid/")

dyn.load("/home/ggj/Documents/glpk-5.0/src/.libs/libglpk.so.40") 
Sys.getenv("LD_LIBRARY_PATH")
Sys.setenv(LD_LIBRARY_PATH=paste0(Sys.getenv("LD_LIBRARY_PATH"), ":", "/home/ggj/Documents/glpk-5.0/src/.libs/"))

library(data.table)
library(Seurat)
library(reshape2)

auc <- fread("./pseudocell/results/Macrophage.auc_all.csv")
auc <- as.data.frame(auc)
auc[1:5,1:5]
rownames(auc) <- auc$Cell
auc <- auc[,-1]

anno <- read.csv("./pseudocell/macrophage_anno.csv",row.names = 1)
anno <- anno[rownames(auc),]

pbmc <- CreateSeuratObject(t(auc),meta.data = anno,min.cells = 0,min.features = 0)
dim(pbmc)

unique(pbmc$CellType)
pbmc$cluster <- colsplit(pbmc$CellType,"_",names = c("n1","n2"))$n2

Idents(pbmc) <- pbmc$cluster
avg <- AverageExpression(pbmc)
avg <- as.data.frame(avg)
colnames(avg) <- gsub("RNA.","",colnames(avg))
colnames(avg) <- gsub("cluster0","Kupffer cell",colnames(avg))
colnames(avg) <- gsub("cluster10","Macrophage CXCL9+",colnames(avg))
colnames(avg) <- gsub("cluster1","Macrophage FABP5+",colnames(avg))
colnames(avg) <- gsub("cluster7","Macrophage APOE+",colnames(avg))

pdf("./figure3/TF_pheatmap.pdf",width = 15,height = 10)
pheatmap::pheatmap(t(avg),scale = "column",clustering_method = "ward.D2",cluster_rows = T,border_color = NA,
                   fontsize = 8)
dev.off()



###rss
rss <- read.csv("./pseudocell/results/rss.type.csv",row.names = 1)
rss[1:5,1:5]
rss <- scale(rss)
pheatmap::pheatmap(rss,clustering_method = "ward.D2",cluster_cols = T,cluster_rows = T,border_color = NA,
                   fontsize = 8)


##plot
RSS<-read.csv("./pseudocell/results/rss.type.csv",row.names = 1)
# RSS$V1<-as.character(RSS$V1)
# RSS[1,1]<-"cluster"
# rownames(RSS)<-RSS$V1
# RSS<-RSS[,-1]

RSS<-as.data.frame(t(RSS))
# RSS$cluster<-as.character(RSS$cluster)
# rownames(RSS)<-RSS$cluster
# RSS<-RSS[,-1]

RSS<-na.omit(RSS)
rownames(RSS) <- gsub("\\...","",rownames(RSS))
#colnames(RSS)<-paste0("Cluster",colnames(RSS))
TF<-as.data.frame(rownames(RSS))
celltype<-unique(colnames(RSS))

dir.create("./pseudocell/results/figures/")
library(ggplot2)
library(ggrepel)
for (i in 1:length(colnames(RSS))) {
  temp<-(RSS[,i])
  temp<-as.data.frame(temp)
  colnames(temp)<-as.character(colnames(RSS)[i])
  temp<-cbind(TF,temp)
  colnames(temp)[1]<-"TF"
  colnames(temp)[2]<-"score"
  temp<-temp[order(temp$score,decreasing = T),]
  temp$rank<-rep(1:length(rownames(temp)))
  temp$text<-NA
  temp$text[1:10]<-as.character(temp$TF[1:10])
  temp$color<-"color2"
  temp$color[1:10]<-"color1"
  temp$score<-as.numeric(as.character(temp$score))
  filename<-paste0("./pseudocell/results/figures/",celltype[i],".pdf")
  p<-ggplot(temp,aes(x=rank,y=score))+
    geom_point(aes(color=color),size=1)+
    xlab("Regulons")+ylab("Specificity score")+
    geom_label_repel(aes(label=text,color=color),
                     nudge_x =0.05,force = 300)+
    scale_color_manual(values=c("#002060","#ADB9CA"))+
    scale_x_continuous(breaks=0:500*100)+
    ggtitle(celltype[i])+theme_bw()+
    theme(legend.position="none",plot.title = element_text(hjust = 0.5),
          axis.title.x =element_text(size=10), axis.title.y=element_text(size=10),
          axis.text=element_text(size=10))
  p
  ggsave(p,filename = filename,width = 4,height = 6)
}


#####
scenicfile <- paste0("./pseudocell/results/Macrophage.auc_all.csv")

scenic <- fread(scenicfile)
scenic <- as.data.frame(scenic)

#scenic <- dcast(scenic,TF~Tissue,mean)
rownames(scenic) <- scenic$Cell
scenic <- scenic[,-1]
scenic <- t(scenic)
scenic[1:5,1:5]
pbmc <- CreateSeuratObject(scenic,min.cells=0,min.features = 3)

anno <- read.csv("./pseudocell/macrophage_anno.csv",row.names = 1)
anno <- anno[colnames(scenic),]
library(reshape2)
colnames(anno) <- c("cellID","celltype")
# anno$type <- paste0(anno$stage,"-",anno$Final_lineage)
# pbmc$type <- anno$type
pbmc$celltype <- anno$celltype
Idents(pbmc) <- pbmc$celltype
table(pbmc$celltype)
avg <- AverageExpression(pbmc)
avg <- as.data.frame(avg$RNA)
TF <- unique(rownames(avg))
avg <- avg[TF,]
avg1 <- scale(t(avg))
avg1 <- as.data.frame(t(avg1))
avg1$TF <- rownames(avg1)

avg.use <- melt(avg1)

avg.TF<-avg.use[order(avg.use$variable,-avg.use$value  ),]
library(dplyr)
top20 <- avg.TF %>% group_by(variable) %>% top_n(15, value)

write.csv(top20,file="./pseudocell/results/auc_scale_top15.csv")

library(gplots)
scale_top20 <- read.csv("./pseudocell/results/auc_scale_top15.csv",row.names = 1)
unique(scale_top20$variable)
use <- scale_top20
avg1 <- avg[use$TF,unique(use$variable)]
avg1 <- scale(t(avg1))
max(avg1);min(avg1)
avg1 <- as.data.frame(avg1)
avg1 <- avg1[c(3,1,2,4),]
library(pheatmap)
pdf("./pseudocell/results/figures/heatmap.pdf",width = 6,height = 10)
pheatmap(t(avg1),clustering_method = 'ward.D2',
         cluster_cols= F,cellwidth = 30,cellheight = 8,fontsize = 8,
         col = colorpanel(100,low="blue",mid = "black", high="yellow"),
         border_color = 'NA')

dev.off()
