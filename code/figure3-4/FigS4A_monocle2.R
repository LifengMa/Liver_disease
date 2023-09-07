#Fig S4A
setwd("/media/ggj/ggjlab2/hezuo/gwh/myeloid/")

library(Seurat)
library(monocle)

load("./myeloid_harmony_1129.RData")
dim(pbmc.harmony)
pbmc <- pbmc.harmony
sample_sheet <- pbmc@meta.data
celltype <- unique(sample_sheet$celltype)
celltype.use <- celltype[c(2:6,8:12,14)]
sample_sheet <- sample_sheet[sample_sheet$celltype%in%celltype.use,]
cellname <- sample(rownames(sample_sheet),10000)
sample_sheet <- sample_sheet[cellname,]

expr_matrix <- GetAssayData(pbmc@assays$RNA,slot="counts")
expr_matrix <- expr_matrix[,rownames(sample_sheet)]
gene_annotation<-data.frame(rownames(expr_matrix),row.names = rownames(expr_matrix))
colnames(gene_annotation)<-"gene_short_name"

pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)
HSMM <- newCellDataSet(as.matrix(expr_matrix),
                       phenoData = pd, featureData = fd,
                       #expressionFamily = VGAM::gaussianff()
)
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
HSMM <- detectGenes(HSMM, min_expr = 0.1)
print(head(fData(HSMM)))
print(head(pData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 5))
print(head(pData(HSMM)))
valid_cells <- row.names(subset(pData(HSMM)))
pData(HSMM)$Total_mRNAs <- Matrix::colSums(exprs(HSMM))
HSMM <- HSMM[,pData(HSMM)$Total_mRNAs < 1e6]
upper_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) +
                     2*sd(log10(pData(HSMM)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) -
                     2*sd(log10(pData(HSMM)$Total_mRNAs)))
qplot(Total_mRNAs, data = pData(HSMM), color = type, geom =
        "density") +
  geom_vline(xintercept = lower_bound) +
  geom_vline(xintercept = upper_bound)
HSMM <- HSMM[,pData(HSMM)$Total_mRNAs > lower_bound &
               pData(HSMM)$Total_mRNAs < upper_bound]

#######################1

#method1
HSMM_myo<-HSMM
diff_test_res <- differentialGeneTest(HSMM_myo[expressed_genes,],
                                      fullModelFormulaStr = "~type")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.0001))
#method2
HSMM_myo <- detectGenes(HSMM_myo, min_expr = 0.1)
fData(HSMM_myo)$use_for_ordering <-
  fData(HSMM_myo)$num_cells_expressed > 0.05 * ncol(HSMM_myo)
#method3
HSMM_myo<-HSMM
disp_table <- dispersionTable(HSMM_myo[expressed_genes,])
ordering_genes<- subset(disp_table,
                        mean_expression >= 0.1&
                          dispersion_empirical >= 1* dispersion_fit)$gene_id

#method4
library(dplyr)
pbmc.markers <- readxl::read_xlsx("./myeloid_marker_harmony.xlsx")
pbmc.markers <- pbmc.markers[which(pbmc.markers$cluster%in%c(2:6,8:12,14)),]
top50 <- pbmc.markers %>% group_by(cluster) %>% top_n(20, avg_log2FC)
top50gene<-top50$gene[!duplicated(top50$gene)]
ordering_genes<-top50gene
#Next
mitogene <- rownames(pbmc)[grep("^MT-|^RP",rownames(pbmc))]
mitogene2 <- mitogene[grep("^MT-",mitogene)]
ordering_genes <- setdiff(ordering_genes,mitogene)
length(ordering_genes)
HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes)
plot_ordering_genes(HSMM_myo)
HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2)
# GM_state <- function(cds){
#   if(length(unique(pData(cds)$type))>1){
#     T0_counts <- table(pData(cds)$type)["St54"]
#     return((names(T0_counts)[which
#                   (T0_counts==max(T0_counts))]))
#   }
#   else{
#     return(1)
#   }
# }
HSMM_myo <- orderCells(HSMM_myo)
par(mfrow =c(1,2))

load("./color.B.Rdata")
color <- color.B$celltype[celltype.use]

dir.create("./monocle2")
library(RColorBrewer)
col_flg<-colorRampPalette(brewer.pal(8,"Paired"))(15)
names(col_flg) <- levels(pbmc)
col_flg <- col_flg[unique(HSMM_myo$celltype)]

p=plot_cell_trajectory(HSMM_myo, color_by="celltype",show_tree = T,show_branch_points=T,cell_size=1)
p+theme(legend.text=element_text(size = 10),legend.position="right")+scale_color_manual(values=col_flg)


plot_cell_trajectory(HSMM_myo, color_by="Pseudotime",show_tree = T,show_branch_points=F,cell_size=1)+
  theme(legend.text=element_text(size = 10),legend.position="right")
p2 <- plot_complex_cell_trajectory(HSMM_myo,color_by="celltype")+scale_color_manual(values=col_flg)
p2
plot_cell_trajectory(HSMM_myo,color_by="State")

######reorder
HSMM_myo2 <- orderCells(HSMM_myo, root_state = 1)
HSMM_myo2 <- HSMM_myo


p1 <- plot_cell_trajectory(HSMM_myo2, color_by="celltype",show_tree = T,show_branch_points=T,cell_size=1)+
  theme(legend.text=element_text(size = 10),legend.position="none",
        panel.border = element_blank())+scale_color_manual(values=col_flg)
p1

p2 <- plot_cell_trajectory(HSMM_myo2, color_by="Pseudotime",show_tree = T,show_branch_points=F,cell_size=1)+
  theme(legend.text=element_text(size = 10),legend.position="right")
p2

HSMM_myo2$type <- HSMM_myo2$celltype
unique(HSMM_myo2$type)
HSMM_myo2$type <- gsub("Plasmacytoid Dendritic Cell_GZMB\\+|Plasmacytoid Dendritic Cell","pDC",HSMM_myo2$type)
HSMM_myo2$type <- gsub("Conventional Dendritic Cell","cDC",HSMM_myo2$type)
HSMM_myo2$type <- gsub("Monocyte_S100A9\\+|Monocyte_FCGR3B\\+|Monocyte_DNAJB1\\+|Monocyte_IL1B\\+|Monocyte_CXCL10\\+","Monocyte",HSMM_myo2$type)
HSMM_myo2$type <- gsub("Macropahge_CLU\\+|Macropahge_FABP5\\+","Macrophage",HSMM_myo2$type)

col_flg2 <- col_flg<-colorRampPalette(brewer.pal(4,"Set1"))(4)
p3 <- plot_cell_trajectory(HSMM_myo2, color_by="type",show_tree = T,show_branch_points=T,cell_size=1)+
  theme(legend.text=element_text(size = 10),legend.position="right")+scale_color_manual(values=col_flg2)
p3

p1|p2|p3
ggsave("./monocle2/pseudotime.pdf",width = 18,height = 6)


save(HSMM_myo2,file="./monocle2/monocle2.RData")
###################2

diff_test_res <- differentialGeneTest(HSMM_myo2[ordering_genes],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
diff_test_res[,c("gene_short_name", "pval", "qval")]
sig_gene_names <- row.names(subset(diff_test_res,  qval < 1e-6 & pval<1e-6 & num_cells_expressed >200))
length(sig_gene_names )
sig_gene_names
dir.create("./monocle2")
pdf("./monocle2/heatmap.pdf",  width = 10,height = 25)
plot_pseudotime_heatmap(HSMM_myo2[sig_gene_names,],
                        num_clusters = 3,
                        cores = 1,
                        show_rownames = T,
                        return_heatmap = T)
dev.off()

####
BEAM_res <- BEAM(HSMM_myo2[expressed_genes,], branch_point = 1, cores = 8) 
#这里用的是ordergene，也就是第六步dpFeature找出来的基因。如果前面用的是seurat的marker基因，记得改成express_genes
#BEAM_res <- BEAM(cds, branch_point = 1, cores = 2) #对2829个基因进行排序，运行慢
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
#BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
head(BEAM_res)
sig_gene_names <- row.names(subset(BEAM_res,  qval < 1e-6 & pval<1e-6 &num_cells_expressed>200))
sig_gene_names <- setdiff(sig_gene_names,mitogene2)
pdf("./monocle2/heatmap.pdf",width = 8,height = 20)
plot_genes_branched_heatmap(HSMM_myo2[sig_gene_names,],
                            branch_point = 1, #绘制的是哪个分支
                            num_clusters = 4, #分成几个cluster，根据需要调整
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)#有632个gene，太多了
dev.off()
save(HSMM_myo,diff_test_res,file = "./Monocle.RData")
save.image(file="./monocle2/monocle2.RData")

p <- plot_genes_branched_heatmap(HSMM_myo2[sig_gene_names,],
                                 branch_point = 2, #绘制的是哪个分支
                                 num_clusters = 4, #分成几个cluster，根据需要调整
                                 cores = 1,
                                 use_gene_short_name = T,
                                 show_rownames = F)#有632个gene，太多了
module <- factor(paste0('cluster',cutree(p$tree_row,3)))
annotation_col <- data.frame(module)
rownames(annotation_col) <- rownames(HSMM_myo[sig_gene_names,])
annotation_col$gene <- rownames(annotation_col)
write.csv(annotation_col,"./monocle2/gene_module.csv")
