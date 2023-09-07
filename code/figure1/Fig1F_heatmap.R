#Figure 1F
setwd("/media/ggj/ggjlab2/hezuo/gwh/")

library(ggplot2)
library(reshape2)
library(tidyr)
library(reshape)

anno <- read.csv("./anno_datu_1205.csv",row.names = 1)
anno.use <- anno[,c(1,6,8)]
anno.use <- anno.use[!duplicated(anno.use),]
colnames(anno.use)[1] <- "cluster"

anno.use <- anno.use[order(anno.use$celltype2),]

cluster <- unique(anno.use$celltype3)
all <- as.data.frame(table(anno$type))

use <- NULL
anno$type <- factor(anno$type)
for(i in 1:length(cluster)){
  temp <- anno[anno$celltype3==cluster[i],]
  summary <- as.data.frame(table(temp$type))
  summary$clusternum <- length(rownames(temp))
  summary <- cbind(summary,all)
  
  summary <- summary[,c(1:3,5)]
  colnames(summary) <- c("type","typenum","clusternum","allnum")
  
  ratio <- (summary$typenum/summary$clusternum)/(summary$allnum/sum(summary$allnum))
  
  use.temp <- data.frame(ratio=ratio,type=summary$type,cluster=cluster[i])
  use <- rbind(use,use.temp)
}

use$cluster <- factor(use$cluster,levels=cluster)
use2 <- use
use2[use2$ratio>5,]$ratio <- 5
use2$ratio <- round(use2$ratio,1)
use2$type <- gsub("AD","Adj",use2$type)
use2$type <- gsub("Fibrosis","Cirrhosis",use2$type)
ggplot(use2,aes(x=cluster,y=type))+geom_tile(aes(fill=ratio))+theme_classic()+
  scale_fill_gradient(name="ratio",low='blue',high='red')+
  geom_text(aes(label=ratio),color="white",size=2)+
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())

###Pheatmap
colnames(use2)[1] <- "value"
use3 <- cast(use2,type~cluster)
rownames(use3) <- use3$type
use3 <- use3[,-1]

pdf("./figures/heatmap.pdf",width = 15,height = 4)
bk <- c(seq(0,0.9,by=0.01),seq(1,5,by=0.05))
pheatmap::pheatmap(use3,cluster_rows = F,cluster_cols = F,cellwidth = 18,cellheight = 18,
                   display_numbers = T,fontsize_number = 8,show_colnames = F,color = c(colorRampPalette(colors = c("#FFFACD","white"))(length(bk)/2),colorRampPalette(colors = c("white","darkgreen"))(length(bk)/2)),
                   breaks = bk)
dev.off()
