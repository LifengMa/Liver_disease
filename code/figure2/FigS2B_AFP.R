#Figure S2B
setwd("/media/ggj/ggjlab2/hezuo/gwh/cancer")

dyn.load("/home/ggj/Documents/glpk-5.0/src/.libs/libglpk.so.40") 
Sys.getenv("LD_LIBRARY_PATH")
Sys.setenv(LD_LIBRARY_PATH=paste0(Sys.getenv("LD_LIBRARY_PATH"), ":", "/home/ggj/Documents/glpk-5.0/src/.libs/"))

library(Seurat)
library(reshape2)
library(ggplot2)
library(openxlsx)
#library(BiocManager)
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)

load("./cancer.RData")
Idents(pbmc) <- pbmc$type
DimPlot(pbmc)
HCC <- subset(pbmc,idents = c("HCC"))
DimPlot(HCC)

DefaultAssay(HCC) <- "RNA"
HCC <- NormalizeData(HCC)
p1 <- FeaturePlot(HCC,features = c("AFP"))
cellname <- WhichCells(HCC,expression = AFP>4)

p2 <- DimPlot(HCC,cells.highlight = cellname)
CombinePlots(plots = list(p1,p2))
ggsave("AFP_usecell.pdf",width = 15,height = 8)

aa <- FindMarkers(HCC,ident.1 = cellname,logfc.threshold = 0,min.pct = 0,min.diff.pct = 0)
write.csv(aa,file="AFP_diffgene_all.csv")

gene <- rownames(aa[aa$avg_log2FC<=(-0.8)&aa$p_val_adj<0.05,])

eg = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
id = as.character(eg[,2])
head(id)
length(gene)
length(id)

##GO
ego <- enrichGO(gene = id,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

use<- ego@result[c(1:20),]
use <- use[use$p.adjust<=0.05,]

###geom_bar
ggplot(data=use)+
  geom_bar(aes(x=reorder(Description,c(Count)),y=Count, fill=-log10(p.adjust)), stat='identity') + 
  coord_flip() +
  scale_fill_gradient(expression(-log10(p.adjust)),low="#81c5f4",high = "firebrick3") +
  xlab("") +
  ylab("Gene count") +
  scale_y_continuous(expand=c(0, 0))+theme_classic()+
  theme(
    axis.text.x=element_text(color="black",size=10),
    axis.text.y=element_text(color="black", size=15),
    axis.title.x = element_text(color="black", size=rel(0.8)),
    legend.text=element_text(color="black",size=rel(1)),
    legend.title = element_text(color="black",size=rel(1))
  )
#geom_point
ggplot(data=use)+
  geom_point(aes(x=reorder(Description,c(Count)),y=Count, color=-log10(p.adjust),size=-log10(p.adjust))) + 
  coord_flip() +
  scale_color_gradient(expression(-log10(p.adjust)),low="#81c5f4",high = "firebrick3") +
  xlab("") +
  ylab("Gene count") +
  scale_y_continuous(expand=c(0, 0))+theme_classic()+
  theme(
    axis.text.x=element_text(color="black",size=10),
    axis.text.y=element_text(color="black", size=15),
    axis.title.x = element_text(color="black", size=rel(0.8)),
    legend.text=element_text(color="black",size=rel(1)),
    legend.title = element_text(color="black",size=rel(1))
  )
ggsave("./AFP-_GO.pdf",width = 10,height = 8)

#volcano
aa$gene <- rownames(aa)

logFC <-aa$avg_log2FC
logFC <- as.numeric(logFC)
adj <- aa$p_val_adj
gene<- aa$gene

data <- data.frame(logFC=logFC,padj=adj,gene=gene)
data$sig[(data$padj > 0.05|data$padj=="NA")|(data$logFC <0.25)& data$logFC > -0.25] <- "no"
data$sig[data$padj <= 0.05 & data$logFC >= 0.25] <- "up"
data$sig[data$padj <= 0.05 & data$logFC <= -0.25] <- "down"

# 
x_lim <- max(logFC,-logFC)
# 
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
#pdf(file = "./pdf/volcano_MDDvsNormal.pdf",width=8,height=8)
theme_set(theme_classic())
data$gene<-as.character(data$gene)
data$type <- data$sig
data$label=ifelse((data$logFC>1|data$logFC<(-1))&data$padj<0.05,data$gene,"")
#data$label=ifelse(data$gene%in%c("GAPDH","Ku70"),data$gene,"")
data$type <- as.character(data$type)
data$logFC <- as.numeric(data$logFC)

p <- ggplot(data,aes(round(logFC,4),-1*log10(padj),
                     color = type))+geom_point()+
  xlim(-5,5) +  labs(x="log2(FoldChange)",y="-log10(pvalue_adj)")
p <- p + scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+
  geom_hline(yintercept=-log10(0.05),linetype=4)+
  geom_vline(xintercept=c(-0.25,0.25),linetype=4)
p <- p +theme(panel.grid =element_blank())+
  theme(axis.line = element_line(size=0))+ylim(0,300)
p <- p  +guides(colour = FALSE)
#p <- p +theme(axis.text=element_text(size=20),axis.title=element_text(size=20))
p<-p+geom_text_repel(data = data, aes(x = data$logFC, 
                                      y = -log10(data$padj), 
                                      label = label),
                     size = 3,box.padding = unit(0.5, "lines"),
                     point.padding = unit(0.8, "lines"), 
                     segment.color = "black", 
                     show.legend = FALSE)
p


ggsave("./volcano_AFP.pdf",width = 12,height = 8)
rm(pbmc,pbmc.markers,seob_epi)
save.image("./AFP_result.RData")
