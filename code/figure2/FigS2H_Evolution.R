#Fig S2H
setwd("/media/ggj/ggjlab2/hezuo/gwh/new_data/")
library(Seurat)
library(infercnv)
library(Seurat)
library(openxlsx)

#Normal
load("./Normal/Normal.RData")
anno <- read.csv("./Normal/Normalanno.csv")

which(anno$V1=="Hepatocyte")
which(anno$V1=="Epithelial")
new.cluster.ids<-as.character(anno$V1)
#new.cluster.ids<-c("Cancer cell1", "T cell2","Plasmocyte3","DC4", "Cycling5", "Plasmocyte6", "Macrophage7", "Cancer cell8", "Monocyte9",
#"T cell10","Unknown11","Macrophage12","Cancer ceell13","Cancer cell14","RPL/S high15","Endothelial16","Cancer cell17","Cancer cell18","Cancer cell19","HSP high20","Plasmocyte21","Unknown22")

names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc@meta.data$celltype<-Idents(pbmc)

pbmc.use <- subset(pbmc,idents=c("Hepatocyte"))
dim(pbmc.use)

##inferCNV

setwd("/media/ggj/ggjlab/RData/202009/liver/inferCNV/")

load("../CA_RData/CA3.RData")
dim(pbmc)
DimPlot(pbmc,label = T)
#CA<-GetAssayData(pbmc@assays$RNA,slot="counts")
FeaturePlot(pbmc,features = c("IGKC"))
DefaultAssay(pbmc) <- "RNA"
pbmc <- NormalizeData(pbmc)
pbmc.markers <- FindAllMarkers(object =pbmc, only.pos = TRUE, min.pct = 0.25, 
                               logfc.threshold  = 0.25)
pbmc.markers<-pbmc.markers[order(pbmc.markers$cluster,-pbmc.markers$avg_log2FC,pbmc.markers$p_val  ),]

#annofile<-read.table("./gencode_v19_gene_pos.txt",sep="\t")
marker<-read.xlsx("../marker/CA3-marker.xlsx",1)
marker<-marker[,c(6,8)]
celltype<-unique(marker)
celltype<-na.omit(celltype)
colnames(celltype)<-c("cluster","anno")
new.cluster.ids<-as.character(celltype$anno)
#new.cluster.ids<-c("Cancer cell1", "T cell2","Plasmocyte3","DC4", "Cycling5", "Plasmocyte6", "Macrophage7", "Cancer cell8", "Monocyte9",
#"T cell10","Unknown11","Macrophage12","Cancer ceell13","Cancer cell14","RPL/S high15","Endothelial16","Cancer cell17","Cancer cell18","Cancer cell19","HSP high20","Plasmocyte21","Unknown22")

names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc@meta.data$celltype<-Idents(pbmc)
unique(Idents(pbmc))
pbmc.CA1 <- subset(pbmc,idents=c("cancer cell"))

pbmc <- merge(pbmc.CA1,pbmc.use)
dim(pbmc)
CA <- as.matrix(GetAssayData(pbmc@assays$RNA,slot = "counts"))

annofile<-as.data.frame(pbmc$celltype)
colnames(annofile)<-"anno"
table(annofile$anno)


#annofile<-as.data.frame(Idents(pbmc))
#annofile$`Idents(pbmc)`<-paste0("Cluster",annofile$`Idents(pbmc)`)
setwd("/media/ggj/ggjlab2/hezuo/gwh/")
dir.create("infercnv")
setwd("./infercnv/")
write.table(annofile,file="./annofile_CA3_normal_evolution.txt",sep="\t",quote = F,col.names = F)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=CA, # 可以直接提供矩阵对象
                                    annotations_file="./annofile_CA3_normal_evolution.txt",
                                    delim="\t",
                                    gene_order_file="/media/ggj/ggjlab/RData/202009/liver/inferCNV//gencode_v19_gene_pos.txt",
                                    ref_group_names=c("Hepatocyte"))

# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.001,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="output_CA3+normal_evolution",  # 输出文件???                             cluster_by_groups=T,   # 聚类
                             denoise=T, #去噪
                             HMM=T,# 是否基于HMM预测CNV
                             cluster_by_groups = F,# 是否基于HMM预测CNV
                             analysis_mode = "subclusters",
                             tumor_subcluster_partition_method = c("random_trees"),
                             num_threads = 16) # 是否基于HMM预测CNV
