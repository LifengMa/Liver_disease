setwd("/media/ggj/ggjlab2/hezuo/gwh/myeloid/")

dyn.load("/home/ggj/Documents/glpk-5.0/src/.libs/libglpk.so.40") 
Sys.getenv("LD_LIBRARY_PATH")
Sys.setenv(LD_LIBRARY_PATH=paste0(Sys.getenv("LD_LIBRARY_PATH"), ":", "/home/ggj/Documents/glpk-5.0/src/.libs/"))

library(Seurat)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

load("./myeloid_harmony_1129.RData")
DimPlot(pbmc.harmony,label=T)

#Dotplot (Fig 3B)
features <- c("CD5L","MARCO","C1QB",
              "IL1B","CXCL2","CCL4",
              "FCER1A","CD1C",
              "SPP1","GPNMB","FABP5",
             
              
              "CD14","VCAN",
              "FCN1",
              "CLEC9A","IDO1","IRF8",
              "APOE","CTSD",
              "FCGR3B","GZMA","NKG7",
              "PCNA","HMGB2","STMN1",
              "CXCL10","CXCL9",
              "PTGDS","LILRA4","IRF7",
              "CAMP","LTF","MMP8",
              "DEFA3","DEFA4","MPO")
DotPlot(pbmc.harmony,features = features,cols = c("blue", "red"),dot.scale = 8)+theme(axis.text.x = element_text(angle=45,hjust=1))+
  ylab("Cluster")+xlab("Genes")
ggsave("./myeloid/Dotplot.pdf",width = 12,height = 6)



#vlion plot (Fig S3BC)
features2 <- c("FABP5","FABP4","MARCO","CD9","PLIN2","APOE","APOC1","CXCL10","CXCL9","SPP1","C1QC","TREM2")
Idents(pbmc.harmony) <- pbmc.harmony$type
unique(pbmc.harmony$type)
sub1 <- subset(pbmc.harmony,idents=c("Normal","CA"))
Idents(sub1) <- sub1$celltype
#cellname <- WhichCells(sub1,idents = 14)
#sub1 <- subset(sub1,cells = setdiff(colnames(sub1),cellname))
VlnPlot(sub1,features = features2,split.by = "type",pt.size = 0,split.plot = T,adjust=1.5,ncol=4)
ggsave("./myeloid/Vlnplot_CA.pdf",width = 25,height = 20)


sub2 <- subset(pbmc.harmony,idents=c("Normal","Fibrosis"))
Idents(sub2) <- sub2$celltype
#cellname <- WhichCells(sub2,idents = 14)[c(1,2)]
#sub2 <- subset(sub2,cells = setdiff(colnames(sub2),cellname))
VlnPlot(sub2,features = features2,split.by = "type",pt.size = 0,split.plot = T,adjust=1.5,ncol=4)
ggsave("./myeloid/Vlnplot_Fib.pdf",width = 25,height = 20)

sub3 <- subset(pbmc.harmony,idents=c("AD","CA"))
Idents(sub3) <- sub3$celltype
sub3$type <- factor(sub3$type,levels = c("CA","AD"))
#cellname <- WhichCells(sub1,idents = 14)
#sub1 <- subset(sub1,cells = setdiff(colnames(sub1),cellname))
VlnPlot(sub3,features = features2,split.by = "type",pt.size = 0,split.plot = T,adjust=1.5,ncol=4)
ggsave("./myeloid/Vlnplot_CA_Adj.pdf",width = 25,height = 20)
# volcano (Fig S4B)-----------------------------------------------------------------
Idents(pbmc.harmony) <- pbmc.harmony$seurat_clusters
DimPlot(pbmc.harmony,label=T)
aa <- FindMarkers(pbmc.harmony,ident.1 = 3,ident.2 = 0,logfc.threshold = 0,min.pct = 0,min.diff.pct = 0)
aa.use <- aa[aa$p_val_adj<0.05&abs(aa$avg_log2FC)>0.25,]
write.csv(aa.use,"./diffgene_TAMvsKupffer.csv")
#volcano
aa$avg_log2FC <- (-aa$avg_log2FC)
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
data$label=ifelse(data$gene%in%c("FABP5","VEGFA"),data$gene,data$label)

data$type <- as.character(data$type)
data$logFC <- as.numeric(data$logFC)

p <- ggplot(data,aes(round(logFC,4),-1*log10(padj),
                     color = type))+geom_point()+
  xlim(-3.5,3.5) +  labs(x="log2(FoldChange)",y="-log10(pvalue_adj)")
p <- p + scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+
  geom_hline(yintercept=-log10(0.05),linetype=4)+
  geom_vline(xintercept=c(-0.25,0.25),linetype=4)
p <- p +theme(panel.grid =element_blank())+
  theme(axis.line = element_line(size=0))+ylim(0,350)
p <- p  +guides(colour = FALSE)
#p <- p +theme(axis.text=element_text(size=20),axis.title=element_text(size=20))
p1<-p+geom_text_repel(data = data, aes(x = data$logFC, 
                                      y = -log10(data$padj), 
                                      label = label),
                     size = 3,box.padding = unit(0.5, "lines"),
                     point.padding = unit(0.8, "lines"), 
                     segment.color = "black", 
                     show.legend = FALSE,
                     max.overlaps = 30)
p1


ggsave("./myeloid/volcano_TAMvsKupffer.pdf",width = 15,height = 10)

###GO (Fig S4C)
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
aa <- read.csv("./diffgene_TAMvsKupffer.csv",row.names = 1)
cluster.use <- aa[aa$p_val_adj<0.05&aa$avg_log2FC<(-0.8),]
cluster.use <- aa[aa$p_val_adj<0.05&aa$avg_log2FC>(0.8),]
symbol=as.character(rownames(cluster.use))
eg = bitr(symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
id = as.character(eg[,2])
head(id)
#length(gene)
length(id)

##GO
ego <- enrichGO(gene = id,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
View(ego@result)
#Kupffer
use<- ego@result[c("GO:0006956","GO:0006959","GO:0002253","GO:0006958"),]
use<- ego@result[c("GO:0001666","GO:0032496","GO:0050900","GO:2001233","GO:0045765","GO:1901342","GO:0006096"),]

head(use)
use$p.adjust_log <- -log10(use$p.adjust)


cluster1 <- use
cluster1$cluster <- "Marcophage FABP5+"
cluster2 <- use
cluster2$cluster <- "Kupffer cell"
cluster2$p.adjust_log <- -cluster2$p.adjust_log
use <- rbind(cluster1,cluster2)
use$p.adjust_log <- abs(use$p.adjust_log)
# cluster7 <- use
# cluster7$cluster <- "APOE+ macrophage"
#use<-use[order(aaa$p.value ),]
#dotplot(ego)
#library(topGO)
#plotGOgraph(ego)
library(ggpubr)
ggdotchart(use,x="Description",y="p.adjust_log",color = "cluster",
           palette = c("#FFA500","#08519C"),group = "cluster",
           add="segments",add.params = list(color="lightgray",size=1.5),
           dot.size = 5,sorting = "descending",rotate = T,
           ggtheme = theme_pubclean())+ylab("-log10(p.adjust)")+xlab("GO term")+
  theme(axis.text.y = element_text(size=12))
ggsave("./myeloid/volcano_GO_dotchart.pdf",width = 10,height = 6)




# use$hjust <- ifelse(use$cluster=="Kupffer cell",0,1)
# use$nudge_y <- ifelse(use$cluster=="Kupffer cell",-0.01,0.01)
# use <- use[order(use$p.adjust_log),]
# use$Description <- factor(use$Description,levels=use$Description)
# 
# ggplot(use, aes(Description, p.adjust_log,fill=cluster)) + 
#   geom_bar(stat = 'identity',alpha = 0.7) + 
#   scale_fill_manual(values = c("#FFA500","#08519C"))+
#   geom_text(data = use, aes(label = Description, y = nudge_y),
#             nudge_x =0,nudge_y =0,hjust =use$hjust,
#             size = 5)+
#   labs(x=("GO term"),
#     y=("-log10(Pvalue_adjust)"))+
#   scale_y_continuous(limits=c(-20,20))+
#   coord_flip() + 
#   theme_classic() + 
#   # theme(panel.grid =element_blank())+
#   # theme(panel.border = element_rect(size = 0.6)
#   #       #panel.border = element_blank()
#   # )+
#   theme(plot.title = element_text(hjust = 0.5,size = 18),
#         axis.text.y = element_blank(),
#         axis.title = element_text(hjust = 0.5,size = 18),
#         axis.line = element_blank(),
#         axis.ticks.y = element_blank()
#         #legend.position = limt
#   )

# GSVA (Fig S4D)--------------------------------------------------------------------


library(GSVA)
library(msigdbr)
Idents(pbmc.harmony) <- pbmc.harmony$seurat_clusters
sub3 <- subset(pbmc.harmony,idents = c(0,3))
expr <- as.matrix(sub3@assays$RNA@data)

#reactome
msgdC2 <- msigdbr(species = "Homo sapiens",category = "C2",subcategory = "REACTOME")
Reactomeset <- msgdC2%>% split(x=.$gene_symbol,f=.$gs_description)
Reactome <- gsva(expr,gset.idx.list =Reactomeset,kcdf="Gaussian",method="zscore",parallel.sz=8,min.sz>4 )

#hallmark
msgdH <- msigdbr(species = "Homo sapiens",category = "H")
Hset <- msgdH%>% split(x=.$gene_symbol,f=.$gs_description)
HALLMARK <- gsva(expr,gset.idx.list =Hset,kcdf="Gaussian",method="zscore",parallel.sz=8,min.sz>4 )


library(limma)
source("./limma.R")
meta <- sub3@meta.data$seurat_clusters
meta <- gsub("3","Macrophage",meta)
meta <- gsub("0","Kupffer",meta)
Diff <- de_gsva(HALLMARK,meta,compare = "Macrophage-Kupffer")

idiff <- Diff[["Macrophage-Kupffer"]]
df <- data.frame(ID=rownames(idiff),score=idiff$t)

df$group <- sapply(1:nrow(idiff),function(x){
  if(idiff[x,"logFC"]>0& idiff[x,"adj.P.Val"]<0.05){return("up")}
  else if(idiff[x,"logFC"]<0 &idiff[x,"adj.P.Val"]<0.05){return("down")}
  else{return("noSig")}
})

df$hjust <- ifelse(df$score>0,1,0)
df$nudge_y <- ifelse(df$score>0,-0.1,0.1)

sortdf <- df[order(df$score),]
sortdf$ID <- factor(sortdf$ID,levels=sortdf$ID)

sortdf.use <- sortdf[c(1,2,5,8,38,39,43,44,45,46,48,49,50),]
#sortdf.use$ID <- gsub("\\.","",sortdf.use$ID)
limt = max(abs(df$score))
ggplot(sortdf.use, aes(ID, score,fill=group)) + 
  geom_bar(stat = 'identity',alpha = 0.7) + 
  scale_fill_manual(breaks=c("down","noSig","up"),
                    values = c("#FFA500","grey","#08519C"))+
  geom_text(data = sortdf.use, aes(label = ID, y = nudge_y),
            nudge_x =0,nudge_y =0,hjust =sortdf.use$hjust,
            size = 5)+
  labs(
    y=("t value of GSVA score"))+
  scale_y_continuous(limits=c(-limt,limt))+
  coord_flip() + 
  theme_classic() + 
  # theme(panel.grid =element_blank())+
  # theme(panel.border = element_rect(size = 0.6)
  #       #panel.border = element_blank()
  # )+
  theme(plot.title = element_text(hjust = 0.5,size = 18),
        axis.text.y = element_blank(),
        axis.title = element_text(hjust = 0.5,size = 18),
        axis.line = element_blank(),
        axis.ticks.y = element_blank()
        #legend.position = limt
  )

ggsave("./myeloid/GSVA.pdf",width=25,height=10)
write.csv(df,"HALLMARK_GSVA_result.csv")

####addmodule (Fig S4E)

##Genes encoding proteins involved in metabolism of fatty acids
genelist1 <- list(intersect(Hset[["Genes defining epithelial-mesenchymal transition, as in wound healing, fibrosis and metastasis."]],rownames(pbmc.harmony)))
pbmc.harmony <- AddModuleScore(pbmc.harmony,features = genelist1,name="EMT")
#pbmc.harmony$fatty1 <- scale(pbmc.harmony$fatty1)
mydata <- FetchData(pbmc.harmony,vars=c("UMAP_1","UMAP_2","hypoxia1"))
ggplot(mydata,aes(x=UMAP_1,y=UMAP_2,color=hypoxia1))+geom_point(size=0.5)+scale_color_gradientn(values=seq(0,1,0.2),colours=c('blue','grey','yellow','red'))

FeaturePlot(pbmc.harmony,features = "hypoxia1",cols =   c("#00ff00","lightgrey","#ff0000"),pt.size = 0.5)
use <- data.frame(score=pbmc.harmony$EMT1,celltype=pbmc.harmony$celltype,cluster=pbmc.harmony$seurat_clusters)
use <- use[which(use$cluster%in%c(0,3,7,9)),]
library(RColorBrewer)
col_flg<-colorRampPalette(brewer.pal(4,"Set1"))(4)
ggplot(use,aes(x=celltype,y=score,fill=celltype))+geom_boxplot()+theme_classic()+scale_fill_manual(values = col_flg)+
  theme(axis.text.x = element_text(angle = 45,hjust=1),legend.position = 'none')+ggtitle("EMT")
ggsave("./myeloid/GSVA_score_EMT.pdf",width = 4,height = 6)


##
genelist1 <- list(intersect(Hset[["Genes up-regulated during formation of blood vessels (angiogenesis)."]],rownames(pbmc.harmony)))
pbmc.harmony <- AddModuleScore(pbmc.harmony,features = genelist1,name="angiogenesis")
#pbmc.harmony$fatty1 <- scale(pbmc.harmony$fatty1)
mydata <- FetchData(pbmc.harmony,vars=c("UMAP_1","UMAP_2","angiogenesis1"))
ggplot(mydata,aes(x=UMAP_1,y=UMAP_2,color=angiogenesis1))+geom_point(size=0.5)+scale_color_gradientn(values=seq(0,1,0.2),colours=c('blue','grey','yellow','red'))
pbmc.harmony@meta.data[,"angiogenesis1"] <- ifelse(pbmc.harmony@meta.data[,"angiogenesis1"] > quantile(pbmc.harmony@meta.data[,"angiogenesis1"],0.75),"High","Low")
DimPlot(pbmc.harmony,group.by  = "angiogenesis1")

use <- data.frame(score=pbmc.harmony$angiogenesis1,celltype=pbmc.harmony$celltype,cluster=pbmc.harmony$seurat_clusters)
use <- use[which(use$cluster%in%c(0,3,7,9)),]
library(RColorBrewer)
col_flg<-colorRampPalette(brewer.pal(4,"Set1"))(4)
ggplot(use,aes(x=celltype,y=score,fill=celltype))+geom_boxplot()+theme_classic()+scale_fill_manual(values = col_flg)+
  theme(axis.text.x = element_text(angle = 45,hjust=1),legend.position = 'none')+ggtitle("Angiogenesis")
ggsave("./myeloid/GSVA_score_angiogenesis.pdf",width = 4,height = 6)

##
genelist1 <- list(intersect(Hset[["Genes up-regulated in response to low oxygen levels (hypoxia)."]],rownames(pbmc.harmony)))
pbmc.harmony <- AddModuleScore(pbmc.harmony,features = genelist1,name="hypoxia")
#pbmc.harmony$fatty1 <- scale(pbmc.harmony$fatty1)
mydata <- FetchData(pbmc.harmony,vars=c("UMAP_1","UMAP_2","hypoxia1"))
ggplot(mydata,aes(x=UMAP_1,y=UMAP_2,color=hypoxia1))+geom_point(size=0.5)+scale_color_gradientn(values=seq(0,1,0.2),colours=c('blue','grey','yellow','red'))

pbmc.harmony@meta.data[,"hypoxia1"] <- ifelse(pbmc.harmony@meta.data[,"hypoxia1"] > quantile(pbmc.harmony@meta.data[,"hypoxia1"],0.75),"High","Low")
DimPlot(pbmc.harmony,group.by  = "hypoxia1")
use <- data.frame(score=pbmc.harmony$hypoxia1,celltype=pbmc.harmony$celltype,cluster=pbmc.harmony$seurat_clusters)
use <- use[which(use$cluster%in%c(0,3,7,9)),]
library(RColorBrewer)
col_flg<-colorRampPalette(brewer.pal(4,"Set1"))(4)
ggplot(use,aes(x=celltype,y=score,fill=celltype))+geom_boxplot()+theme_classic()+scale_fill_manual(values = col_flg)+
  theme(axis.text.x = element_text(angle = 45,hjust=1),legend.position = 'none')+ggtitle("hypoxia")
ggsave("./myeloid/GSVA_score_hypoxia.pdf",width = 4,height = 6)
