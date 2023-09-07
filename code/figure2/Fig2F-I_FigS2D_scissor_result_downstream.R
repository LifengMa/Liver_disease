setwd("/media/ggj/ggjlab2/hezuo/gwh/cancer/Scissor/")
library(Seurat)
library(ggplot2)
library(reshape2)
library(data.table)
library(gprofiler2)
library(openxlsx)
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)

load('./Scissor_HCC_survival_0.25_result.RData')
DimPlot(HCC, reduction = 'tsne', group.by = 'scissor', cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))

DefaultAssay(HCC) <- "RNA"
HCC <- NormalizeData(HCC)
aa <- FindMarkers(HCC,group.by = "scissor",ident.1 = '1',ident.2 = '2',min.pct = 0,min.diff.pct = 0,logfc.threshold = 0)
write.csv(aa,file="./diffgene_all.csv")
# volcano (Figure 2H)-----------------------------------------------------------------

#volcano
#aa$avg_log2FC <- (-aa$avg_log2FC)
aa$gene <- rownames(aa)

logFC <-aa$avg_log2FC
logFC <- as.numeric(logFC)
adj <- aa$p_val_adj
gene<- aa$gene

data <- data.frame(logFC=logFC,padj=adj,gene=gene)
data$sig[(data$padj > 0.05|data$padj=="NA")|(data$logFC <0.25)& data$logFC > -0.25] <- "no"
data$sig[data$padj <= 0.05 & data$logFC >= 0.25] <- "up"
data$sig[data$padj <= 0.05 & data$logFC <= -0.24] <- "down"

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
data$label=ifelse((data$logFC>1.5|data$logFC<(-1.5))&data$padj<0.05,data$gene,"")
#data$label=ifelse(data$gene%in%c("FABP5","VEGFA"),data$gene,data$label)

data$type <- as.character(data$type)
data$logFC <- as.numeric(data$logFC)

p <- ggplot(data,aes(round(logFC,4),-1*log10(padj),
                     color = type))+geom_point()+
  xlim(-3,3) +  labs(x="log2(FoldChange)",y="-log10(pvalue_adj)")
p <- p + scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+
  geom_hline(yintercept=-log10(0.05),linetype=4)+
  geom_vline(xintercept=c(-0.25,0.25),linetype=4)
p <- p +theme(panel.grid =element_blank())+
  theme(axis.line = element_line(size=0))+ylim(0,250)
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


ggsave("./volcano.pdf",width = 8,height = 6)


# GO (Figure 2I)----------------------------------------------------------------------

bb <- FindMarkers(HCC,group.by = "scissor",ident.1 = '1',ident.2 = '2')
cluster.use <- bb[bb$p_val_adj<0.05&bb$avg_log2FC>0.8,]
symbol=as.character(rownames(cluster.use))
eg = bitr(symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
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
View(ego@result)

use<- ego@result[c("GO:0001666","GO:0006096","GO:0045765","GO:2001233"),]

head(use)
use$p.adjust_log <- -log10(use$p.adjust)


#use<-use[order(aaa$p.value ),]
#dotplot(ego)
#library(topGO)
#plotGOgraph(ego)
library(ggpubr)
ggdotchart(use,x="Description",y="p.adjust_log",add="segments",color = "indianred1",add.params = list(color="lightgray",size=1),
           dot.size = 5,sorting = "descending",rotate = T,
           ggtheme = theme_pubclean())+
  ylab("-log10(p.adjust)")+xlab("GO term")
ggsave("./GO.pdf",width = 6,height = 4)

# vlnplot (Figure 2G)-----------------------------------------------------------------
library(ggpubr)
library(gridExtra)
features <- c("AKR1B10","LDHA","ANGPTL8","ENO1")
fig <- list()
for (i in 1:length(features)) {
  temp <- subset(HCC,features=features[i])
  temp <- as.data.frame(GetAssayData(temp@assays$RNA,slot="data"))
  temp <- as.data.frame(t(temp))
  colnames(temp) <- 'Expression level'
  temp$type <- HCC$scissor
  temp <- temp[temp$type%in%c(1,2),]
  temp$type <- gsub("1","Scissor+ cell",temp$type)
  temp$type <- gsub("2","Scissor- cell",temp$type)
  
  comparision <- list(c("Scissor+ cell","Scissor- cell"))
  p <- ggplot(temp,aes(x=type,y=`Expression level`,fill=type))+geom_violin()+
    scale_fill_manual(values=c('royalblue','indianred1'))+
    stat_compare_means(comparisons = comparision,method="wilcox.test")+ggtitle(features[i])+
    theme(legend.position = 'none')
  
  fig[[i]] <- p
  }
fig[['nrow']] <- 2
fig[['ncol']] <- 2

pdf("./vlnplot.pdf",width = 8,height = 8)
do.call('grid.arrange', fig)
dev.off()

save.image("./downstream_result.RData")

#################fraction###############
library(ggplot2)
use <- as.data.frame(table(HCC$scissor,HCC$patient))
use <- use[use$Var1==1,]
use <- use[use$Freq>0,]
use$Var2 <- gsub("CA","HCC",use$Var2)
data <- use

dt = data[order(data$Freq, decreasing = TRUE),]
myLabel = as.vector(dt$Var2)   
myLabel = paste(myLabel, "(", round(dt$Freq / sum(dt$Freq) * 100, 2), "%)", sep = "")   
dt$Var1 <- myLabel

library(RColorBrewer)
HCC$patient <- gsub("CA","HCC",HCC$patient)
col <- colorRampPalette(brewer.pal(6,"Set1"))(14)
DimPlot(HCC,group.by = "patient",reduction = "tsne",cols = col)

names(col) <- levels(as.factor(HCC$patient))
col.use <- col[use$Var2]
p = ggplot(dt, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 1) +    
  coord_polar(theta = "y") + 
  theme_bw()+
  labs(x = "", y = "", title = "") + 
  theme(axis.ticks = element_blank()) + 
  theme(legend.title = element_blank(), legend.position = "top") + 
  #scale_fill_discrete(breaks = dt$Freq, labels = myLabel) + 
  theme(axis.text.x = element_blank())+scale_fill_manual(values=as.character(col.use))
#geom_text(aes(y = Freq/2 + c(0, cumsum(Freq)[-length(Freq)]), x = sum(Freq)/20, label = myLabel), size = 5)   ## 在图中加上百分比：x 调节标签到圆心的距离, y 调节标签的左右位置
p
ggsave("./patient_piechart.pdf",width = 8,height = 8)


################signature gene logfc>1 & pval <0.05 (Figure S2D)#######
load("./downstream_result.RData")
rm(list=setdiff(ls(),c("bb")))
diffgene <- bb[bb$p_val_adj<=0.05&bb$pct.1>0.15,]
diffgene <- diffgene[order(diffgene$avg_log2FC,decreasing = T),]
write.csv(diffgene,"./diffgene_Scissor.csv")
HCC_signature <-  bb[bb$avg_log2FC>=0.25&bb$p_val_adj<0.05&bb$pct.1>0.25,]
HCC_signature <- list(rownames(HCC_signature))
HCC_signature 

###load TCGA LIHC data
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
TCGA_clinical.use = TCGA_clinical.use[grepl("Stage",TCGA_clinical.use$ajcc_pathologic_stage),]
TCGA_clinical.use$ajcc_pathologic_stage = gsub("[ABCD]$","",TCGA_clinical.use$ajcc_pathologic_stage)
TCGA_clinical.use$case_submitter_id <- gsub("-","\\.",TCGA_clinical.use$case_submitter_id)
TCGA_clinical.use <- merge(TCGA_clinical.use,TCGA_anno,by.x='case_submitter_id',by.y='anno')

TCGA_clinical.use2 <- TCGA_clinical.use[,c(1,6,7)]
TCGA_clinical.use2$case_submitter_id <- gsub("-","\\.",TCGA_clinical.use2$case_submitter_id)
table(TCGA_clinical.use$ajcc_pathologic_stage)

###ssgsea
library(GSVA)
TCGA <- TCGA[,TCGA_clinical.use$colnames.aa.]
ssgsea <- gsva(as.matrix(TCGA), HCC_signature, method="ssgsea", kcdf="Gaussian", abs.ranking=TRUE,parallel.sz=20)
ssgsea <- gsva(as.matrix(TCGA), HCC_signature, method="gsva", kcdf="Gaussian",parallel.sz=20)
ssgsea <- as.data.frame(t(ssgsea))
ssgsea$sample <- rownames(ssgsea)
colnames(ssgsea)[1] <- 'HCC_WS_signature'

ssgsea <- merge(ssgsea,TCGA_clinical.use,by.x='sample',by.y='colnames.aa.')

ggplot(ssgsea, aes(x = ajcc_pathologic_stage, y = HCC_WS_signature,fill = ajcc_pathologic_stage)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2,size=0)+
  theme_classic()+ylab("HCC_WS_signature")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position = "none")

#MCP counter
library(MCPcounter)
DEGs <- data.frame(DEGs = HCC_signature[[1]],cellType = 'HCC')

MCPEstimates=MCPcounter.estimate(TCGA,featuresType="affy133P2_probesets",probesets=DEGs)
MCPEstimates <- as.data.frame(t(MCPEstimates))
MCPEstimates$sample <- rownames(MCPEstimates)
MCPEstimates$HCC <- scale(MCPEstimates$HCC)
colnames(MCPEstimates)[1] <- 'HCC_WS_signature'

MCPEstimates <- merge(MCPEstimates,TCGA_clinical.use,by.x='sample',by.y='colnames.aa.')
ggplot(MCPEstimates, aes(x = ajcc_pathologic_stage, y = HCC_WS_signature,fill = ajcc_pathologic_stage)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2,size=0)+
  theme_classic()+ylab("HCC_WS_signature")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position = "none")

# smoothHR analysis:
library(phenoTest)
ssgsea$OS <- as.numeric(ssgsea$OS)
pdf("HCC_WS_signature_cox.pdf",width = 10,height = 8)
smoothCoxph(ssgsea$OS,ssgsea$status,ssgsea$HCC_WS_signature,xlab="HCC_WS_signature",main="HCC_WS_signature")
dev.off()
gene.cox <- coxph(Surv(OS, status==1) ~ HCC_WS_signature, ssgsea)
ggforest(gene.cox,data=ssgsea)
summary(gene.cox)

ssgsea$stage <- ssgsea$ajcc_pathologic_stage
ssgsea$stage <- gsub("Stage IV","4",ssgsea$stage)
ssgsea$stage <- gsub("Stage III","3",ssgsea$stage)
ssgsea$stage <- gsub("Stage II","2",ssgsea$stage)
ssgsea$stage <- gsub("Stage I","1",ssgsea$stage)

ssgsea$stage <- as.numeric(ssgsea$stage)
gene.cox2 <- coxph(Surv(OS, status==1) ~ HCC_WS_signature+stage, ssgsea)
ggforest(gene.cox2,data=ssgsea)
summary(gene.cox2)
#Survival
library(survminer)
library(survival)
surv <- ssgsea
surv$OS <- as.numeric(surv$OS)
res.cut <- surv_cutpoint(surv,
                         time="OS",
                         event="status",
                         variables = c( "HCC_WS_signature"))
rownames(surv) <- surv$sample
surv1 <- surv_categorize(res.cut)
surv1$sample <- surv$sample
surv <- surv1

anno <- surv

attach(anno)
Surv(OS,status)
fit <- survfit(Surv(OS,status)~HCC_WS_signature,
               data=anno)
fit


#######validation in ICGC 
ICGC_specimen <- read.table("/media/ggj/ggjlab2/hezuo/gwh/ICGC/specimen.tsv.gz",sep="\t",header=T)

ICGC <- fread("/media/ggj/ggjlab2/hezuo/gwh/ICGC/exp_seq.tsv.gz")
ICGC[1:5,1:22]
unique(ICGC$submitted_sample_id)
ICGC <- ICGC[grep("Cancer",ICGC$submitted_sample_id),]

ICGC_clinical <- read.table("/media/ggj/ggjlab2/hezuo/gwh/ICGC/donor.tsv.gz",sep="\t",header=T)
rownames(ICGC_clinical) <- ICGC_clinical$icgc_donor_id
ICGC_clinical <- ICGC_clinical[which(ICGC_clinical$donor_survival_time>=0&ICGC_clinical$donor_tumour_stage_at_diagnosis%in%c(1,2,3,4)),]
ICGC_clinical$event <- ifelse(ICGC_clinical$donor_vital_status=="alive",0,1) 

ICGC <- ICGC[ICGC$icgc_donor_id%in%ICGC_clinical$icgc_donor_id,]
head(ICGC)
ICGC <- as.data.frame(ICGC)
ICGC.use <- ICGC[,c(1,8,9)]
class(ICGC.use)
ICGC.use <- as.data.table(ICGC.use)
ICGC.use$type <- paste0(ICGC.use$icgc_donor_id,"_",ICGC.use$gene_id)
ICGC.use[1:5,1:4]
length(unique(ICGC.use$type))
ICGC.use <- ICGC.use[!duplicated(ICGC.use$type),]
ICGC.use <- dcast(ICGC.use,gene_id~icgc_donor_id,value.var = 'normalized_read_count',drop=FALSE, fill=0)
ICGC.use[1:5,1:5]
ICGC.use <- as.data.frame(ICGC.use)
rownames(ICGC.use) <- ICGC.use$gene_id
ICGC.use <- ICGC.use[,-1]
ICGC.use <- ICGC.use[rowSums(ICGC.use)>0,]
###ssgsea
library(GSVA)
ssgsea2 <- gsva(as.matrix(ICGC.use), HCC_signature, method="gsva", kcdf="Gaussian",parallel.sz=20)
ssgsea2 <- gsva(as.matrix(ICGC.use), HCC_signature,method="ssgsea", kcdf="Gaussian", abs.ranking=TRUE,parallel.sz=20)
ssgsea2 <- as.data.frame(t(ssgsea2))
ssgsea2$sample <- rownames(ssgsea2)
colnames(ssgsea2)[1] <- 'HCC_WS_signature'

ssgsea2 <- merge(ssgsea2,ICGC_clinical,by.x='sample',by.y='icgc_donor_id')
#ssgsea2 <- ssgsea2[ssgsea2$donor_tumour_stage_at_diagnosis%in%c(1,2,3,4),]
colnames(ssgsea2)
ggplot(ssgsea2, aes(x = donor_tumour_stage_at_diagnosis, y = HCC_WS_signature,fill = donor_tumour_stage_at_diagnosis)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2,size=0)+
  theme_classic()+ylab("HCC_WS_signature")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position = "none")

#MCP counter
library(MCPcounter)
DEGs <- data.frame(DEGs = HCC_signature[[1]],cellType = 'HCC')

MCPEstimates2=MCPcounter.estimate(ICGC.use,featuresType="affy133P2_probesets",probesets=DEGs)
MCPEstimates2 <- as.data.frame(t(MCPEstimates2))
MCPEstimates2$sample <- rownames(MCPEstimates2)
MCPEstimates2$HCC <- scale(MCPEstimates2$HCC)
colnames(MCPEstimates2)[1] <- 'HCC_WS_signature'

MCPEstimates2 <- merge(MCPEstimates2,ICGC_clinical,by.x='sample',by.y='icgc_donor_id')
ggplot(MCPEstimates2, aes(x = donor_tumour_stage_at_diagnosis, y = HCC_WS_signature,fill = donor_tumour_stage_at_diagnosis)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2,size=0)+
  theme_classic()+ylab("HCC_WS_signature")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position = "none")
MCPEstimates2$OS <- MCPEstimates2$donor_survival_time
# smoothHR analysis:
library(phenoTest)
ssgsea2$OS <- ssgsea2$donor_survival_time
ssgsea2$OS <- as.numeric(ssgsea2$OS)
pdf("../../ICGC/ICGC_WS_signature_cox.pdf",width = 10,height = 8)
smoothCoxph(ssgsea2$OS,ssgsea2$event,ssgsea2$HCC_WS_signature,xlab="HCC_WS_signature",main="HCC_WS_signature")
#smoothCoxph(MCPEstimates2$OS,MCPEstimates2$event,MCPEstimates2$HCC_WS_signature,xlab="HCC_WS_signature",main="HCC_WS_signature")
dev.off()

#Survival
library(survminer)
library(survival)

surv <- ssgsea2
#surv <- MCPEstimates2
surv$OS <- as.numeric(surv$OS)
colnames(surv)
res.cut <- surv_cutpoint(surv,
                         time="OS",
                         event="event",
                         variables = c("HCC_WS_signature"))
plot(res.cut)
rownames(surv) <- surv$sample
surv1 <- surv_categorize(res.cut)
#surv1 <- surv
# upper <- quantile(surv1$HCC_WS_signature,0.75)
# lower <- quantile(surv1$HCC_WS_signature,0.25)
# surv1$HCC_WS_signature1 <- ifelse(surv1$HCC_WS_signature>=upper,"high","middle")
# surv1$HCC_WS_signature1[surv1$HCC_WS_signature<=lower]<-"low"
# surv1 <- surv1[surv1$HCC_WS_signature1%in%c("low","high"),]
surv1$sample <- surv$sample
surv1$stage <- surv$donor_tumour_stage_at_diagnosis
surv1$gender <- surv$donor_sex
surv1$age <- surv$donor_age_at_diagnosis
surv <- surv1

anno <- surv

attach(anno)
Surv(OS,event)
fit <- survfit(Surv(OS,event==1)~HCC_WS_signature,
               data=anno)
fit
pdf("../../ICGC/ICGC_survival.pdf",
    width = 10,height = 10)
p1 <- ggsurvplot(fit,data=anno,pval=TRUE,palette="aaas",conf.int = F,surv.median.line="hv",
                 risk.table = T)
p1
print(p1, newpage = FALSE)
dev.off()

##cox
surv_cox <- ssgsea2
surv_cox$sample <- surv_cox$sample
surv_cox$stage <- surv_cox$donor_tumour_stage_at_diagnosis
surv_cox$gender <- surv_cox$donor_sex
surv_cox$age <- surv_cox$donor_age_at_diagnosis
surv_cox$gender <- ifelse(surv_cox$gender=="male",1,2)
surv_cox$stage <- as.numeric(surv_cox$stage)
cox = coxph(Surv(OS,event==1) ~ HCC_WS_signature + stage+age+gender, data=surv_cox)
cox.summary = summary(cox)
cox.summary
ggforest(cox,data=surv_cox,main = "ICGC")
ggsave("../../ICGC/ICGC_WS_signature_forest.pdf",width = 10,height = 6)

rm(ICGC,TCGA)
gc()
save.image("../../ICGC/result.RData")
install.packages('forestplot')
library(forestplot)