setwd("/media/ggj/ggjlab2/hezuo/gwh/")
library(Seurat)
library(ggplot2)
library(openxlsx)
library(ggthemes)
library(ggsci)
library(scales)
library(reshape2)
dir.create("./scripts")

###cell proportion (Figure 1E)
anno <- read.csv("./anno_datu_1205.csv",row.names = 1)
col <- c('#1f77b4',
         '#ff7f0e',
         '#279e68',
         '#d62728',
         '#aa40fc',
         '#8c564b',
         '#e377c2',
         '#b5bd61',
         '#17becf',
         '#aec7e8',
         '#ffbb78',
         '#98df8a')
unique(anno$patient)
anno$patient <- colsplit(anno$batch,"_",names=c("n1","n2"))$n1
anno$patient <- gsub("Fib","Cirrhosis",anno$patient)
anno$patient <- gsub("Liver","Normal",anno$patient)
anno$patient <- gsub("CA","HCC",anno$patient)
anno$patient <- gsub("AD","Adj",anno$patient)
anno$patient <- gsub("Cirrhosis3S1","Cirrhosis3",anno$patient)
anno$patient <- gsub("Cirrhosis3S2","Cirrhosis4",anno$patient)
anno$patient <- factor(anno$patient,levels = c("HCC1","HCC2","HCC3","HCC4","HCC5","HCC6","HCC7S1",
                                               "HCC7S2","HCC8","HCC9","HCC10","ICC1","ICC2","ICC3",
                                               "Cirrhosis1S1","Cirrhosis1S2","Cirrhosis2","Cirrhosis3","Cirrhosis4",
                                               "HCC8Adj","HCC9Adj","HCC10Adj","ICC1Adj","ICC2Adj","ICC3Adj",
                                               "Normal1","Normal2","Normal4","Normal5S1","Normal5S2","Normal6S1","Normal6S2","Normal7S1","Normal7S2"))
write.csv(anno,file="anno_datu_1205.csv")
ggplot(anno,aes(x=patient,fill=celltype2))+geom_bar(position = "fill",lwd=0,width=0.6,color="white")+
  theme_classic()+ylab("Cell proportion")+xlab("Patient ID")+labs(fill="Cell type")+
  scale_y_continuous(expand = c(0.01,0))+scale_fill_manual(values=col)+
  theme(axis.text.x = element_text(angle = 45,hjust=1),
        axis.line = element_blank())
ggsave("./figures/barplot_0707.pdf",width = 8,height = 4)
