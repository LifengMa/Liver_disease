#Fig S2G

setwd("/media/ggj/ggjlab/RData/202009/liver/CA_RData/")

library(Seurat)
library(infercnv)
library(reshape2)
library(openxlsx)

load("./CA7S1.RData")
marker<-read.xlsx("../marker/CA7S1-marker.xlsx",1)
marker<-marker[,c(6,8)]
celltype<-unique(marker)
celltype<-na.omit(celltype)
colnames(celltype)<-c("cluster","anno")
new.cluster.ids<-as.character(celltype$anno)
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc@meta.data$celltype<-Idents(pbmc)
unique(Idents(pbmc))

data <- as.data.frame(GetAssayData(pbmc@assays$RNA,slot="counts"))
anno <- as.data.frame(Idents(pbmc))
colnames(anno) <- "CellType"

dir.create("./Pseudocell_infercnv/")
fwrite(anno,file="./Pseudocell_infercnv/CA7S1_anno.csv",row.names = T)
fwrite(data,file="./Pseudocell_infercnv/CA7S1_dge.csv",row.names = T)

##pseudocell
source("/media/ggj/ggjlab2/index_microscopy/zhongqi/pseudocell_function.R")

Pseudocell_analysis_pipeline(
  "./Pseudocell_infercnv/CA7S1_dge.csv",
  "./Pseudocell_infercnv/CA7S1_anno.csv", 
  
  sep=",",
  pseudocell_size = 10,
  colname_str = "CA7S1_",
  save_out_file = TRUE, 
  out_file_name_str = "./Pseudocell_infercnv/CA7S1") 


#####infercnv

CA <- readRDS("./Pseudocell_infercnv/CA7S1_pseudocell_10_.Rds")
anno <- read.csv("./Pseudocell_infercnv/CA7S1_pseudocell_10_.pheno.csv",row.names = 1)
rownames(anno) <- colnames(CA)
annofile<-as.data.frame(anno$Celltype)
colnames(annofile)<-"anno"
rownames(annofile) <- colnames(CA)
table(annofile$anno)

write.table(annofile,file="/media/ggj/ggjlab2/hezuo/gwh/infercnv/annofile_CA1_evolution_10_1020.txt",sep="\t",quote = F,col.names = F)


#colnames(CA)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=CA, # 可以直接提供矩阵对象
                                    annotations_file="/media/ggj/ggjlab2/hezuo/gwh/infercnv/annofile_CA1_evolution_10_1020.txt",
                                    delim="\t",
                                    gene_order_file="/media/ggj/ggjlab/RData/202009/liver/inferCNV/gencode_v19_gene_pos.txt",
                                    ref_group_names=setdiff(unique(annofile$anno),"CA7S1_cancer cell"))
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             #out_dir="output_CA1+normal_evolution",  # 输出文件???                             cluster_by_groups=T,   # 聚类
                             denoise=T, #去噪
                             HMM=T,# 是否基于HMM预测CNV
                             cluster_by_groups = F,# 是否基于HMM预测CNV
                             analysis_mode = "subclusters",
                             tumor_subcluster_partition_method = c("random_trees"),
                             num_threads = 24,
                             out_dir='/media/ggj/ggjlab2/hezuo/gwh/infercnv/CA7S1_10') # 是否基于HMM预测CNV


