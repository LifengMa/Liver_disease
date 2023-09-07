#Figure S2A
setwd("/media/ggj/ggjlab2/hezuo/gwh/cancer/")

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

DimPlot(pbmc,label=T)

Idents(pbmc) <- pbmc$seurat_clusters
DimPlot(pbmc,label=T)

#marker <- read.csv("./CA_ICC_SCTmarkers.csv",row.names = 1)

DefaultAssay(pbmc) <- "RNA"
pbmc <- NormalizeData(pbmc)
pbmc.markers <- FindAllMarkers(object =pbmc, only.pos = TRUE, min.pct = 0.25, 
                               logfc.threshold  = 0.25)
pbmc.markers<-pbmc.markers[order(pbmc.markers$cluster,-pbmc.markers$avg_log2FC,pbmc.markers$p_val  ),]
WriteXLS::WriteXLS(pbmc.markers, "./cancer_RNAmarkers.xlsx",row.names=T)

###
anno <- pbmc@meta.data
table(anno$seurat_clusters)
anno.cancer <- anno
#anno.cancer <- anno[grep("cancer",anno$anno2),]
final <- NULL
for(i in 1:length(unique(anno.cancer$seurat_clusters))){
  temp <- anno.cancer[anno.cancer$seurat_clusters==unique(anno.cancer$seurat_clusters)[i],]
  marker <- pbmc.markers[pbmc.markers$cluster==unique(temp$seurat_clusters),]
  gene <- marker[marker$avg_log2FC>=1&marker$p_val_adj<0.05,]$gene
  
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
  
  use<- ego@result[c(1:5),]
  if(unique(temp$seurat_clusters)%in%c(4,5,7,12)){
    use$cluster <- paste0("ICC cancer_",unique(temp$seurat_clusters))
  }
  else{
    use$cluster <- paste0("HCC cancer_",unique(temp$seurat_clusters))
  }
  
  final <- rbind(final,use)
}

library(RColorBrewer)
col_flg<-colorRampPalette(brewer.pal(8,"Paired"))(length(unique(final$cluster)))
final1 <- final
unique(final1$cluster)
final1$cluster <- factor(final1$cluster,levels = c(paste0("HCC cancer_",c(0,1,2,3,6,8,9,10,11,13)),paste0("ICC cancer_",c(4,5,7,12))))
ggplot(data=final1)+
  geom_bar(aes(x=reorder(Description,-log10(pvalue)),y=-log10(pvalue), fill=cluster), stat='identity') + 
  coord_flip() +facet_wrap(~cluster,scales = "free",nrow=7)+
  scale_fill_manual(values = col_flg)+
  #scale_fill_gradient(expression(-log10(p.adjust)),low="#81c5f4",high = "firebrick3") +
  xlab("") +
  ylab("-log10(pvalue)") +
  scale_y_continuous(expand=c(0, 0))+theme_classic()+
  theme(
    axis.text.x=element_text(color="black",size=10),
    axis.text.y=element_text(color="black", size=15),
    axis.title.x = element_text(color="black", size=rel(0.8)),
    legend.text=element_text(color="black",size=rel(1)),
    legend.title = element_text(color="black",size=rel(1))
  )
ggsave("./GO_cancercell2.pdf",width = 25,height=20)
