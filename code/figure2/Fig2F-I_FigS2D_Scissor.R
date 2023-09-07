setwd("/media/ggj/ggjlab2/hezuo/gwh/cancer/Scissor/")

library(Scissor)
library(reshape2)
library(ggplot2)

TCGA <- read.table("/media/ggj/ggjlab2/hezuo/gwh/myeloid/CIBERSORT/TCGA.txt.gz",header = T,row.names = 1)
TCGA[1:5,1:5]
TCGA_anno <- read.csv("/media/ggj/ggjlab/RData/202009/liver/TCGA/LIHCanno.csv",row.names = 1) 
TCGA_clinical <- read.delim("/media/ggj/ggjlab/RData/202009/liver/TCGA/clinical.tsv")
###grep tumor sample
TCGA_anno <- TCGA_anno[TCGA_anno$sample%in%c("01A"),]
#filter data
TCGA.use <- TCGA[,TCGA_anno$colnames.aa.]
#process clinical information
TCGA_clinical.use <- TCGA_clinical[,c("case_submitter_id","days_to_death","days_to_last_follow_up","vital_status","ajcc_pathologic_stage")] 
TCGA_clinical.use <- TCGA_clinical.use[!duplicated(TCGA_clinical.use),]
TCGA_clinical.use$days_to_death <- gsub("'--","0",TCGA_clinical.use$days_to_death)
TCGA_clinical.use$days_to_last_follow_up <- gsub("'--","0",TCGA_clinical.use$days_to_last_follow_up)
TCGA_clinical.use$OS <- apply(cbind(TCGA_clinical.use$days_to_death,TCGA_clinical.use$days_to_last_follow_up),1,max)
TCGA_clinical.use$status <- ifelse(TCGA_clinical.use$vital_status=="Dead",1,0)
TCGA_clinical.use2 <- TCGA_clinical.use[,c(1,6,7)]
TCGA_clinical.use2$case_submitter_id <- gsub("-","\\.",TCGA_clinical.use2$case_submitter_id)

###order data
length(unique(TCGA_anno$anno))#369
length(intersect(TCGA_clinical.use2$case_submitter_id,TCGA_anno$anno))#369
setdiff(TCGA_clinical.use2$case_submitter_id,TCGA_anno.use$anno)
rownames(TCGA_clinical.use2) <- TCGA_clinical.use2$case_submitter_id
TCGA_clinical.use3 <- TCGA_clinical.use2[TCGA_anno$anno,]
setdiff(TCGA_clinical.use2$case_submitter_id,TCGA_clinical.use3$case_submitter_id)

#check
all(colnames(TCGA.use)==TCGA_anno$colnames.aa.)
colnames(TCGA.use)[1:5]
rownames(TCGA_clinical.use3)[1:5]

#####
phenotype <- TCGA_clinical.use3[,2:3]
colnames(phenotype) <- c("time","status")
phenotype$time <- as.integer(phenotype$time)
phenotype$status <- as.double(phenotype$status)
head(phenotype)

load("../cancer.RData")
pbmc@graphs$RNA_snn <- pbmc@graphs$SCT_snn
Idents(pbmc) <- pbmc$type
DimPlot(pbmc,label=T)
HCC <- subset(pbmc,idents = "CA")
DimPlot(HCC,reduction = "tsne")
infosl <- Scissor(TCGA.use,GetAssayData(HCC@assays$RNA,slot="counts"),phenotype,alpha = NULL,cutoff = 0.25,family = "cox",Save_file = "Scissor_HCC_survival_0.25.RData")

length(infosl$Scissor_pos)

DimPlot(HCC,label = T)
Scissor_select <- rep(0,ncol(HCC))
names(Scissor_select) <- colnames(HCC)
Scissor_select[infosl$Scissor_pos] <- 1
Scissor_select[infosl$Scissor_neg] <- 2
HCC <- AddMetaData(HCC, metadata = Scissor_select, col.name = "scissor")
DimPlot(HCC, reduction = 'tsne', group.by = 'scissor', cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))
ggsave("./Scissor_result.pdf",width=8,height = 6)
save(infosl,HCC,file="./Scissor_HCC_survival_0.25_result.RData")
###reliability significance test
load("./Scissor_HCC_survival_0.25.RData")
numbers <- length(infosl$Scissor_pos)+length(infosl$Scissor_neg)
result1 <- reliability.test(X,Y,network,alpha = 0.01,family = "cox",cell_num = numbers,n=10,nfold=10)
