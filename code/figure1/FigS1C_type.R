#Figure S1C
setwd("/media/ggj/ggjlab2/hezuo/gwh/")

library(reshape2)
library(ggplot2)

emb <- read.csv("./datu_embeddings.csv",row.names = 1)
colnames(emb) <- c("tSNE_1","tSNE_2")
anno <- read.csv("./anno_datu_1205.csv",row.names = 1)
use <- cbind(anno,emb)
unique(use$type)
#ADJ
use$ADJ <- use$type
use$ADJ <- ifelse(use$ADJ=="AD","ADJ","NO")
ggplot(use,aes(x=tSNE_1,y=tSNE_2,col=ADJ))+geom_point(size=0.1)+
  scale_color_manual(values=c("#D62728","grey"))+ 
  theme(legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())+ggtitle("Adj")
ggsave("./figures/ADJ.pdf",width = 8,height = 8)
#Cirrhosis
use$Cirrhosis <- use$type
use$Cirrhosis <- ifelse(use$Cirrhosis=="Fibrosis","Cirrhosis","NO")
ggplot(use,aes(x=tSNE_1,y=tSNE_2,col=Cirrhosis))+geom_point(size=0.1)+
  scale_color_manual(values=c("#FF7F0E","grey"))+ 
  theme(legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())+ggtitle("Cirrhosis")
ggsave("./figures/Cirrhosis.pdf",width = 8,height = 8)
#HCC
use$HCC <- use$type
use$HCC <- ifelse(use$HCC=="HCC","HCC","NO")
ggplot(use,aes(x=tSNE_1,y=tSNE_2,col=HCC))+geom_point(size=0.1)+
  scale_color_manual(values=c("#17BECF","grey"))+ 
  theme(legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())+ggtitle("HCC")
ggsave("./figures/HCC.pdf",width = 8,height = 8)
#ICC
use$ICC <- use$type
use$ICC <- ifelse(use$ICC=="ICC","ICC","NO")
ggplot(use,aes(x=tSNE_1,y=tSNE_2,col=ICC))+geom_point(size=0.1)+
  scale_color_manual(values=c("#B5BD61","grey"))+ 
  theme(legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())+ggtitle("ICC")
ggsave("./figures/ICC.pdf",width = 8,height = 8)
#Normal
use$Normal <- use$type
use$Normal <- ifelse(use$Normal=="Normal","Normal","NO")
ggplot(use,aes(x=tSNE_1,y=tSNE_2,col=Normal))+geom_point(size=0.1)+
  scale_color_manual(values=c("grey","#E377C2"))+ 
  theme(legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())+ggtitle("Normal")
ggsave("./figures/Normal.pdf",width = 8,height = 8)
