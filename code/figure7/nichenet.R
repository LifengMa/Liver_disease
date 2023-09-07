#library(devtools)

#install_local("/home/ggj/Documents/nichenetr-master.zip")

setwd("/media/ggj/ggjlab2/hezuo/gwh/")
dir.create("./nichenet")
setwd("./nichenet/")

library(nichenetr)
library(Seurat)
library(tidyverse)
library(future)
plan("multisession",workers=1)
options(future.globals.maxSize=100000*1024^2)
rm(list=ls())
gc()

ligand_target_matrix <- readRDS("./dataset/ligand_target_matrix.rds")
lr_network = readRDS("./dataset/lr_network.rds")
weighted_networks = readRDS("./dataset/weighted_networks.rds")

scRNA <- readRDS("./seuratobject.rds")
head(scRNA)

Idents(scRNA) <- "celltype2"
unique(scRNA$celltype2)
table(scRNA$celltype2,scRNA$type2)
unique(scRNA$type2)
nichenet_output = nichenet_seuratobj_aggregate(scRNA, 
                                               top_n_ligands = 30,
                                               receiver = c("Epithelial"), 
                                               sender = c("Kupffer cell", "pDC", "cDC",
                                                          "Monocyte_DNAJB1+", "Monocyte_IL1B+","Monocyte_S100A9+",
                                                          "Monocyte_CXCL10+","Monocyte_FCGR3B+","Macropahge_FABP5+","Macrophage",
                                                          "Neutrophil","Macropahge_CLU+"),
                                               condition_colname = "type2", 
                                               condition_oi = c("CA"), 
                                               condition_reference = c("AD","Normal"), 
                                               ligand_target_matrix = ligand_target_matrix, 
                                               lr_network = lr_network, 
                                               weighted_networks = weighted_networks, 
                                               organism = "human",
                                               assay_oi="RNA",
                                               lfc_cutoff = 0.1,
                                               geneset = 'up')

p = nichenet_output$ligand_target_heatmap
p
pdf("./CAvsADN_TAE.pdf",width = 30,height = 10)
nichenet_output$ligand_activity_target_heatmap
dev.off()

p = DotPlot(scRNA , 
            features = nichenet_output$top_receptors, 
            split.by = "type2") + RotatedAxis()
p



############
seurat_obj <- scRNA
seurat_obj@meta.data$celltype_aggregate = paste(seurat_obj@meta.data$celltype2, seurat_obj@meta.data$type2,sep = "_") # user adaptation required on own dataset
DimPlot(seurat_obj, group.by = "celltype_aggregate")
seurat_obj@meta.data$celltype_aggregate %>% table() %>% sort(decreasing = TRUE)
celltype_id = "celltype_aggregate" # metadata column name of the cell type of interest
seurat_obj = SetIdent(seurat_obj, value = seurat_obj[[celltype_id]])

lr_network = lr_network %>% mutate(bonafide = ! database %in% c("ppi_prediction","ppi_prediction_go"))
lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor, bonafide)

organism = "human"
table(seurat_obj$celltype2,seurat_obj$type2)
###Tcell
niches = list(
  "CA_niche" = list(
    "sender" = c("Kupffer cell_CA", "pDC_CA", "cDC_CA",
                 "Monocyte_DNAJB1+_CA", "Monocyte_IL1B+_CA","Monocyte_S100A9+_CA",
                 "Monocyte_CXCL10+_CA","Monocyte_FCGR3B+_CA","Macropahge_FABP5+_CA","Macrophage_CA",
                 "Neutrophil_CA","Macropahge_CLU+_CA"),
                 #"Activated T cell_CA",
                 #"NK cell_CA","NKT cell_CA","Naive T cell_CA","Naive T cell_C5_CA","CD8 T cell_CA","Treg cell_CA","B cell_C0_CA",
                 #"B cell_C1_CA","B cell_C2_CA","B cell_C3_CA","B cell_C4_CA","B cell_C5_CA"),
    "receiver" = c("Naive T cell_CA","CD8 T cell_CA","Activated T cell_CA","Treg_CA","NK cell_CA",
                   "Epithelial_CA","Fibroblast_CA","TAE_CA","VE_CA"
                   )),
  "AD_niche" = list(
    "sender" = c("Kupffer cell_AD", "pDC_AD", "cDC_AD",
                 "Monocyte_DNAJB1+_AD", "Monocyte_IL1B+_AD","Monocyte_S100A9+_AD",
                 "Monocyte_CXCL10+_AD","Monocyte_FCGR3B+_AD","Macropahge_FABP5+_AD","Macrophage_AD",
                 "Neutrophil_AD","Macropahge_CLU+_AD"),
                 #"Activated T cell_AD",
                 #"NK cell_AD","NKT cell_AD","Naive T cell_AD","Naive T cell_C5_AD","CD8 T cell_AD","Treg cell_AD","B cell_C0_AD",
                 #"B cell_C1_AD","B cell_C2_AD","B cell_C3_AD","B cell_C4_AD","B cell_C5_AD"),
    "receiver" =c(
      "Naive T cell_AD","CD8 T cell_AD","Activated T cell_AD","Treg_AD","NK cell_AD",
      "Epithelial_AD","Fibroblast_AD","TAE_AD","VE_AD","LSEC_AD"
                  )),
  "Cir_niche" = list(
    "sender" = c("Kupffer cell_Cir", "pDC_Cir", "cDC_Cir",
                 "Monocyte_DNAJB1+_Cir", "Monocyte_IL1B+_Cir","Monocyte_S100A9+_Cir",
                 "Monocyte_CXCL10+_Cir","Monocyte_FCGR3B+_Cir","Macropahge_FABP5+_Cir","Macrophage_Cir",
                 "Neutrophil_Cir","Macropahge_CLU+_Cir"),
    #"Activated T cell_CA",
    #"NK cell_CA","NKT cell_CA","Naive T cell_CA","Naive T cell_C5_CA","CD8 T cell_CA","Treg cell_CA","B cell_C0_CA",
    #"B cell_C1_CA","B cell_C2_CA","B cell_C3_CA","B cell_C4_CA","B cell_C5_CA"),
    "receiver" = c("Naive T cell_Cir","CD8 T cell_Cir","Activated T cell_Cir","Treg_Cir","NK cell_Cir",
                   "Epithelial_Cir","Fibroblast_Cir","TAE_Cir","VE_Cir"
    )),
  "Normal_niche" = list(
    "sender" = c("Kupffer cell_Normal", "pDC_Normal", "cDC_Normal",
                 "Monocyte_DNAJB1+_Normal", "Monocyte_IL1B+_Normal","Monocyte_S100A9+_Normal",
                 "Monocyte_CXCL10+_Normal","Monocyte_FCGR3B+_Normal","Macropahge_FABP5+_Normal","Macrophage_Normal",
                 "Neutrophil_Normal","Macropahge_CLU+_Normal"),
    #"Activated T cell_CA",
    #"NK cell_CA","NKT cell_CA","Naive T cell_CA","Naive T cell_C5_CA","CD8 T cell_CA","Treg cell_CA","B cell_C0_CA",
    #"B cell_C1_CA","B cell_C2_CA","B cell_C3_CA","B cell_C4_CA","B cell_C5_CA"),
    "receiver" = c("Naive T cell_Normal","CD8 T cell_Normal","Activated T cell_Normal","Treg_Normal"
    ))
) # user adaptation required on own dataset

###Fibroblast
niches = list(
  "CA_niche" = list(
    "sender" = c("Kupffer cell_CA", "pDC_CA", "cDC_CA",
                 "Monocyte_DNAJB1+_CA", "Monocyte_IL1B+_CA","Monocyte_S100A9+_CA",
                 "Monocyte_CXCL10+_CA","Monocyte_FCGR3B+_CA","Macropahge_FABP5+_CA","Macrophage_CA",
                 "Neutrophil_CA","Macropahge_CLU+_CA"),
    #"Activated T cell_CA",
    #"NK cell_CA","NKT cell_CA","Naive T cell_CA","Naive T cell_C5_CA","CD8 T cell_CA","Treg cell_CA","B cell_C0_CA",
    #"B cell_C1_CA","B cell_C2_CA","B cell_C3_CA","B cell_C4_CA","B cell_C5_CA"),
    "receiver" = c("Fibroblast_CA"
    )),
  "AD_niche" = list(
    "sender" = c("Kupffer cell_AD", "pDC_AD", "cDC_AD",
                 "Monocyte_DNAJB1+_AD", "Monocyte_IL1B+_AD","Monocyte_S100A9+_AD",
                 "Monocyte_CXCL10+_AD","Monocyte_FCGR3B+_AD","Macropahge_FABP5+_AD","Macrophage_AD",
                 "Neutrophil_AD","Macropahge_CLU+_AD"),
    #"Activated T cell_AD",
    #"NK cell_AD","NKT cell_AD","Naive T cell_AD","Naive T cell_C5_AD","CD8 T cell_AD","Treg cell_AD","B cell_C0_AD",
    #"B cell_C1_AD","B cell_C2_AD","B cell_C3_AD","B cell_C4_AD","B cell_C5_AD"),
    "receiver" =c(
      "Fibroblast_AD"
    )),
  "Cir_niche" = list(
    "sender" = c("Kupffer cell_Fibrosis", "pDC_Fibrosis", "cDC_Fibrosis",
                 "Monocyte_DNAJB1+_Fibrosis", "Monocyte_IL1B+_Fibrosis","Monocyte_S100A9+_Fibrosis",
                 "Monocyte_CXCL10+_Fibrosis","Monocyte_FCGR3B+_Fibrosis","Macropahge_FABP5+_Fibrosis","Macrophage_Fibrosis",
                 "Neutrophil_Fibrosis","Macropahge_CLU+_Fibrosis"),
    #"Activated T cell_CA",
    #"NK cell_CA","NKT cell_CA","Naive T cell_CA","Naive T cell_C5_CA","CD8 T cell_CA","Treg cell_CA","B cell_C0_CA",
    #"B cell_C1_CA","B cell_C2_CA","B cell_C3_CA","B cell_C4_CA","B cell_C5_CA"),
    "receiver" = c("Fibroblast_Fibrosis"
    )),
  "Normal_niche" = list(
    "sender" = c("Kupffer cell_Normal", "pDC_Normal", "cDC_Normal",
                 "Monocyte_DNAJB1+_Normal", "Monocyte_IL1B+_Normal","Monocyte_S100A9+_Normal",
                 "Monocyte_CXCL10+_Normal","Monocyte_FCGR3B+_Normal","Macropahge_FABP5+_Normal","Macrophage_Normal",
                 "Neutrophil_Normal","Macropahge_CLU+_Normal"),
    #"Activated T cell_CA",
    #"NK cell_CA","NKT cell_CA","Naive T cell_CA","Naive T cell_C5_CA","CD8 T cell_CA","Treg cell_CA","B cell_C0_CA",
    #"B cell_C1_CA","B cell_C2_CA","B cell_C3_CA","B cell_C4_CA","B cell_C5_CA"),
    "receiver" = c("Fibroblast_Normal"
    ))
) # user adaptation required on own dataset

###all immune
niches = list(
  "CA_niche" = list(
    "sender" = c("Kupffer cell_CA", "pDC_CA", "cDC_CA",
                 "Monocyte_DNAJB1+_CA", "Monocyte_IL1B+_CA","Monocyte_S100A9+_CA",
                 "Monocyte_CXCL10+_CA","Monocyte_FCGR3B+_CA","Macropahge_FABP5+_CA","Macrophage_CA",
                 "Neutrophil_CA","Macropahge_CLU+_CA",
                 "Activated T cell_CA","CD8 T cell_CA","Naive T cell_CA","Treg_CA",
                 "B cell_CA",
                 "NK cell_CA"),
    #"Activated T cell_CA",
    #"NK cell_CA","NKT cell_CA","Naive T cell_CA","Naive T cell_C5_CA","CD8 T cell_CA","Treg cell_CA","B cell_C0_CA",
    #"B cell_C1_CA","B cell_C2_CA","B cell_C3_CA","B cell_C4_CA","B cell_C5_CA"),
    "receiver" = c("TAE_CA","Epithelial_CA","Fibroblast_CA",
                   "VE_CA"
    )),
  "AD_niche" = list(
    "sender" = c("Kupffer cell_AD", "pDC_AD", "cDC_AD",
                 "Monocyte_DNAJB1+_AD", "Monocyte_IL1B+_AD","Monocyte_S100A9+_AD",
                 "Monocyte_CXCL10+_AD","Monocyte_FCGR3B+_AD","Macropahge_FABP5+_AD","Macrophage_AD",
                 "Neutrophil_AD","Macropahge_CLU+_AD",
                 "Activated T cell_AD","CD8 T cell_AD","Naive T cell_AD","Treg_AD",
                 "B cell_AD",
                 "NK cell_AD"),
    #"Activated T cell_AD",
    #"NK cell_AD","NKT cell_AD","Naive T cell_AD","Naive T cell_C5_AD","CD8 T cell_AD","Treg cell_AD","B cell_C0_AD",
    #"B cell_C1_AD","B cell_C2_AD","B cell_C3_AD","B cell_C4_AD","B cell_C5_AD"),
    "receiver" =c(
      "TAE_AD","Epithelial_AD","Fibroblast_AD",
      "VE_AD","LSEC_AD"
    )),
  "Cir_niche" = list(
    "sender" = c("Kupffer cell_Fibrosis", "pDC_Fibrosis", "cDC_Fibrosis",
                 "Monocyte_DNAJB1+_Fibrosis", "Monocyte_IL1B+_Fibrosis","Monocyte_S100A9+_Fibrosis",
                 "Monocyte_CXCL10+_Fibrosis","Monocyte_FCGR3B+_Fibrosis","Macropahge_FABP5+_Fibrosis","Macrophage_Fibrosis",
                 "Neutrophil_Fibrosis","Macropahge_CLU+_Fibrosis",
                 "Activated T cell_Fibrosis","CD8 T cell_Fibrosis","Naive T cell_Fibrosis","Treg_Fibrosis",
                 "B cell_Fibrosis",
                 "NK cell_Fibrosis"),
    #"Activated T cell_CA",
    #"NK cell_CA","NKT cell_CA","Naive T cell_CA","Naive T cell_C5_CA","CD8 T cell_CA","Treg cell_CA","B cell_C0_CA",
    #"B cell_C1_CA","B cell_C2_CA","B cell_C3_CA","B cell_C4_CA","B cell_C5_CA"),
    "receiver" = c("CAE_Fibrosis","Epithelial_Fibrosis","Fibroblast_Fibrosis",
                   "VE_Fibrosis","TAE_Fibrosis","LSEC_Fibrosis"
    )),
  "Normal_niche" = list(
    "sender" = c("Kupffer cell_Normal", "pDC_Normal", "cDC_Normal",
                 "Monocyte_DNAJB1+_Normal", "Monocyte_IL1B+_Normal","Monocyte_S100A9+_Normal",
                 "Monocyte_CXCL10+_Normal","Monocyte_FCGR3B+_Normal","Macropahge_FABP5+_Normal","Macrophage_Normal",
                 "Neutrophil_Normal","Macropahge_CLU+_Normal",
                 "Activated T cell_Normal","CD8 T cell_Normal","Naive T cell_Normal","Treg_Normal",
                 "B cell_Normal",
                 "NK cell_Normal"),
    #"Activated T cell_CA",
    #"NK cell_CA","NKT cell_CA","Naive T cell_CA","Naive T cell_C5_CA","CD8 T cell_CA","Treg cell_CA","B cell_C0_CA",
    #"B cell_C1_CA","B cell_C2_CA","B cell_C3_CA","B cell_C4_CA","B cell_C5_CA"),
    "receiver" = c("CAE_Normal","Epithelial_Normal","Fibroblast_Normal",
                   "VE_Normal","TAE_Normal","LSEC_Normal"
    ))
) # user adaptation required on own dataset

###all immune
niches = list(
  "CA_niche" = list(
    "sender" = c("Kupffer cell_CA", "pDC_CA", "cDC_CA",
                 "Monocyte_DNAJB1+_CA", "Monocyte_IL1B+_CA","Monocyte_S100A9+_CA",
                 "Monocyte_CXCL10+_CA","Monocyte_FCGR3B+_CA","Macropahge_FABP5+_CA","Macrophage_CA",
                 "Neutrophil_CA","Macropahge_CLU+_CA",
                 "Activated T cell_CA","CD8 T cell_CA","Naive T cell_CA","Treg_CA",
                 "B cell_CA",
                 "NK cell_CA"),
    #"Activated T cell_CA",
    #"NK cell_CA","NKT cell_CA","Naive T cell_CA","Naive T cell_C5_CA","CD8 T cell_CA","Treg cell_CA","B cell_C0_CA",
    #"B cell_C1_CA","B cell_C2_CA","B cell_C3_CA","B cell_C4_CA","B cell_C5_CA"),
    "receiver" = c("Epithelial_CA"
    )),
  "AD_niche" = list(
    "sender" = c("Kupffer cell_AD", "pDC_AD", "cDC_AD",
                 "Monocyte_DNAJB1+_AD", "Monocyte_IL1B+_AD","Monocyte_S100A9+_AD",
                 "Monocyte_CXCL10+_AD","Monocyte_FCGR3B+_AD","Macropahge_FABP5+_AD","Macrophage_AD",
                 "Neutrophil_AD","Macropahge_CLU+_AD",
                 "Activated T cell_AD","CD8 T cell_AD","Naive T cell_AD","Treg_AD",
                 "B cell_AD",
                 "NK cell_AD"),
    #"Activated T cell_AD",
    #"NK cell_AD","NKT cell_AD","Naive T cell_AD","Naive T cell_C5_AD","CD8 T cell_AD","Treg cell_AD","B cell_C0_AD",
    #"B cell_C1_AD","B cell_C2_AD","B cell_C3_AD","B cell_C4_AD","B cell_C5_AD"),
    "receiver" =c(
      "Epithelial_AD"
    ))
) # user adaptation required on own dataset

assay_oi = "SCT" # other possibilities: RNA,...
DE_sender = calculate_niche_de(seurat_obj = seurat_obj %>% subset(features = lr_network$ligand %>% unique()), niches = niches, type = "sender", assay_oi = assay_oi) # only ligands important for sender cell types
DE_receiver = calculate_niche_de(seurat_obj = seurat_obj %>% subset(features = lr_network$receptor %>% unique()), niches = niches, type = "receiver", assay_oi = assay_oi) # only receptors now, later on: DE analysis to find targets
DE_sender = DE_sender %>% mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))
DE_receiver = DE_receiver %>% mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))
expression_pct = 0.10
DE_sender_processed = process_niche_de(DE_table = DE_sender, niches = niches, expression_pct = expression_pct, type = "sender")
DE_receiver_processed = process_niche_de(DE_table = DE_receiver, niches = niches, expression_pct = expression_pct, type = "receiver")

specificity_score_LR_pairs = "min_lfc"
DE_sender_receiver = combine_sender_receiver_de(DE_sender_processed, DE_receiver_processed, lr_network, specificity_score = specificity_score_LR_pairs)

###spatial
include_spatial_info_sender = FALSE# if not spatial info to include: put this to false # user adaptation required on own dataset
include_spatial_info_receiver = FALSE # if spatial info to include: put this to true # user adaptation required on own dataset

spatial_info = tibble(celltype_region_oi = "CAF_High", celltype_other_region = "myofibroblast_High", niche =  "pEMT_High_niche", celltype_type = "sender") # user adaptation required on own dataset
specificity_score_spatial = "lfc"

if(include_spatial_info_sender == FALSE & include_spatial_info_receiver == FALSE){
  spatial_info = tibble(celltype_region_oi = NA, celltype_other_region = NA) %>% mutate(niche =  niches %>% names() %>% head(1), celltype_type = "sender")
}
if(include_spatial_info_sender == TRUE){
  sender_spatial_DE = calculate_spatial_DE(seurat_obj = seurat_obj %>% subset(features = lr_network$ligand %>% unique()), spatial_info = spatial_info %>% filter(celltype_type == "sender"))
  sender_spatial_DE_processed = process_spatial_de(DE_table = sender_spatial_DE, type = "sender", lr_network = lr_network, expression_pct = expression_pct, specificity_score = specificity_score_spatial)
  
  # add a neutral spatial score for sender celltypes in which the spatial is not known / not of importance
  sender_spatial_DE_others = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "sender", lr_network = lr_network)
  sender_spatial_DE_processed = sender_spatial_DE_processed %>% bind_rows(sender_spatial_DE_others)
  
  sender_spatial_DE_processed = sender_spatial_DE_processed %>% mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))
  
} else {
  # # add a neutral spatial score for all sender celltypes (for none of them, spatial is relevant in this case)
  sender_spatial_DE_processed = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "sender", lr_network = lr_network)
  sender_spatial_DE_processed = sender_spatial_DE_processed %>% mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))  
  
}

if(include_spatial_info_receiver == TRUE){
  receiver_spatial_DE = calculate_spatial_DE(seurat_obj = seurat_obj %>% subset(features = lr_network$receptor %>% unique()), spatial_info = spatial_info %>% filter(celltype_type == "receiver"))
  receiver_spatial_DE_processed = process_spatial_de(DE_table = receiver_spatial_DE, type = "receiver", lr_network = lr_network, expression_pct = expression_pct, specificity_score = specificity_score_spatial)
  
  # add a neutral spatial score for receiver celltypes in which the spatial is not known / not of importance
  receiver_spatial_DE_others = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "receiver", lr_network = lr_network)
  receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% bind_rows(receiver_spatial_DE_others)
  
  receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
  
} else {
  # # add a neutral spatial score for all receiver celltypes (for none of them, spatial is relevant in this case)
  receiver_spatial_DE_processed = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "receiver", lr_network = lr_network)
  receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
}
####
lfc_cutoff = 0.1 # recommended for 10x as min_lfc cutoff. 
specificity_score_targets = "min_lfc"

DE_receiver_targets = calculate_niche_de_targets(seurat_obj = seurat_obj, niches = niches, lfc_cutoff = lfc_cutoff, expression_pct = expression_pct, assay_oi = assay_oi) 
## [1] "Calculate receiver DE between: Malignant_High and Malignant_Low"
## [1] "Calculate receiver DE between: Malignant_Low and Malignant_High"
DE_receiver_processed_targets = process_receiver_target_de(DE_receiver_targets = DE_receiver_targets, niches = niches, expression_pct = expression_pct, specificity_score = specificity_score_targets)

background = DE_receiver_processed_targets  %>% pull(target) %>% unique()
geneset_niche1 = DE_receiver_processed_targets %>% filter(receiver == niches[[1]]$receiver[1] & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
geneset_niche2 = DE_receiver_processed_targets %>% filter(receiver == niches[[2]]$receiver[1] & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
geneset_niche1 %>% setdiff(rownames(ligand_target_matrix))
geneset_niche2 %>% setdiff(rownames(ligand_target_matrix))
length(geneset_niche1)
length(geneset_niche2)

top_n_target = 250

niche_geneset_list = list(
  "CA_niche" = list(
    "receiver" = niches[[1]]$receiver[2],
    "geneset" = geneset_niche1,
    "background" = background),
  "AD_niche" = list(
    "receiver" = niches[[2]]$receiver[2],
    "geneset" = geneset_niche2 ,
    "background" = background)
)

ligand_activities_targets = get_ligand_activities_targets(niche_geneset_list = niche_geneset_list, ligand_target_matrix = ligand_target_matrix, top_n_target = top_n_target)

####
features_oi = union(lr_network$ligand, lr_network$receptor) %>% union(ligand_activities_targets$target) %>% setdiff(NA)

dotplot = suppressWarnings(Seurat::DotPlot(seurat_obj %>% subset(idents = niches %>% unlist() %>% unique()), features = features_oi, assay = assay_oi))
exprs_tbl = dotplot$data %>% as_tibble()
exprs_tbl = exprs_tbl %>% rename(celltype = id, gene = features.plot, expression = avg.exp, expression_scaled = avg.exp.scaled, fraction = pct.exp) %>%
  mutate(fraction = fraction/100) %>% as_tibble() %>% select(celltype, gene, expression, expression_scaled, fraction) %>% distinct() %>% arrange(gene) %>% mutate(gene = as.character(gene))

exprs_tbl_ligand = exprs_tbl %>% filter(gene %in% lr_network$ligand) %>% rename(sender = celltype, ligand = gene, ligand_expression = expression, ligand_expression_scaled = expression_scaled, ligand_fraction = fraction) 
exprs_tbl_receptor = exprs_tbl %>% filter(gene %in% lr_network$receptor) %>% rename(receiver = celltype, receptor = gene, receptor_expression = expression, receptor_expression_scaled = expression_scaled, receptor_fraction = fraction)
exprs_tbl_target = exprs_tbl %>% filter(gene %in% ligand_activities_targets$target) %>% rename(receiver = celltype, target = gene, target_expression = expression, target_expression_scaled = expression_scaled, target_fraction = fraction)

exprs_tbl_ligand = exprs_tbl_ligand %>%  
  mutate(scaled_ligand_expression_scaled = scale_quantile_adapted(ligand_expression_scaled)) %>% mutate(ligand_fraction_adapted = ligand_fraction) %>% mutate_cond(ligand_fraction >= expression_pct, ligand_fraction_adapted = expression_pct)  %>% mutate(scaled_ligand_fraction_adapted = scale_quantile_adapted(ligand_fraction_adapted))

exprs_tbl_receptor = exprs_tbl_receptor %>% mutate(scaled_receptor_expression_scaled = scale_quantile_adapted(receptor_expression_scaled))  %>% mutate(receptor_fraction_adapted = receptor_fraction) %>% mutate_cond(receptor_fraction >= expression_pct, receptor_fraction_adapted = expression_pct)  %>% mutate(scaled_receptor_fraction_adapted = scale_quantile_adapted(receptor_fraction_adapted))

exprs_sender_receiver = lr_network %>% 
  inner_join(exprs_tbl_ligand, by = c("ligand")) %>% 
  inner_join(exprs_tbl_receptor, by = c("receptor")) %>% inner_join(DE_sender_receiver %>% distinct(niche, sender, receiver))

ligand_scaled_receptor_expression_fraction_df = exprs_sender_receiver %>% group_by(ligand, receiver) %>% mutate(rank_receptor_expression = dense_rank(receptor_expression), rank_receptor_fraction  = dense_rank(receptor_fraction)) %>% mutate(ligand_scaled_receptor_expression_fraction = 0.5*( (rank_receptor_fraction / max(rank_receptor_fraction)) + ((rank_receptor_expression / max(rank_receptor_expression))) ) )  %>% distinct(ligand, receptor, receiver, ligand_scaled_receptor_expression_fraction, bonafide) %>% distinct() %>% ungroup() 

prioritizing_weights = c("scaled_ligand_score" = 5,
                         "scaled_ligand_expression_scaled" = 1,
                         "ligand_fraction" = 1,
                         "scaled_ligand_score_spatial" = 2, 
                         "scaled_receptor_score" = 0.5,
                         "scaled_receptor_expression_scaled" = 0.5,
                         "receptor_fraction" = 1, 
                         "ligand_scaled_receptor_expression_fraction" = 1,
                         "scaled_receptor_score_spatial" = 0,
                         "scaled_activity" = 0,
                         "scaled_activity_normalized" = 1,
                         "bona_fide" = 1)

output = list(DE_sender_receiver = DE_sender_receiver, 
              ligand_scaled_receptor_expression_fraction_df = ligand_scaled_receptor_expression_fraction_df, 
              sender_spatial_DE_processed = sender_spatial_DE_processed, receiver_spatial_DE_processed = receiver_spatial_DE_processed,
              ligand_activities_targets = ligand_activities_targets, DE_receiver_processed_targets = DE_receiver_processed_targets, exprs_tbl_ligand = exprs_tbl_ligand,  exprs_tbl_receptor = exprs_tbl_receptor, exprs_tbl_target = exprs_tbl_target)
prioritization_tables = get_prioritization_tables(output, prioritizing_weights)

#Visualization of the Differential NicheNet output
top_ligand_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% group_by(ligand) %>% top_n(1, prioritization_score) %>% ungroup() %>% select(ligand, receptor, niche) %>% rename(top_niche = niche)
top_ligand_receptor_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% group_by(ligand, receptor) %>% top_n(1, prioritization_score) %>% ungroup() %>% select(ligand, receptor, niche) %>% rename(top_niche = niche)

ligand_prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, prioritization_score) %>% group_by(ligand, niche) %>% top_n(1, prioritization_score) %>% ungroup() %>% distinct() %>% inner_join(top_ligand_niche_df) %>% filter(niche == top_niche) %>% group_by(niche) %>% top_n(50, prioritization_score) %>% ungroup() # get the top50 ligands per niche

receiver_oi = "CD8 T cell_CA" 

filtered_ligands = ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% pull(ligand) %>% unique()

prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% 
  top_n(3, prioritization_score) %>% ungroup() 
lfc_plot = make_ligand_receptor_lfc_plot(receiver_oi, prioritized_tbl_oi, prioritization_tables$prioritization_tbl_ligand_receptor, plot_legend = FALSE, heights = NULL, widths = NULL)
lfc_plot

exprs_activity_target_plot = make_ligand_activity_target_exprs_plot(receiver_oi, prioritized_tbl_oi, 
                                                                    prioritization_tables$prioritization_tbl_ligand_receptor, 
                                                                    prioritization_tables$prioritization_tbl_ligand_target, 
                                                                    output$exprs_tbl_ligand,  output$exprs_tbl_target, lfc_cutoff, 
                                                                    ligand_target_matrix, plot_legend = FALSE, heights = NULL, widths = NULL)
exprs_activity_target_plot$combined_plot

##top20
filtered_ligands = ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% 
  top_n(20, prioritization_score) %>% pull(ligand) %>% unique()

prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
  filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% 
  distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>%
  filter(receiver == receiver_oi) %>% top_n(5, prioritization_score) %>% ungroup() 

exprs_activity_target_plot = make_ligand_activity_target_exprs_plot(receiver_oi, prioritized_tbl_oi,  
                                                                    prioritization_tables$prioritization_tbl_ligand_receptor,
                                                                    prioritization_tables$prioritization_tbl_ligand_target, 
                                                                    output$exprs_tbl_ligand,  output$exprs_tbl_target,
                                                                    lfc_cutoff, ligand_target_matrix, plot_legend = FALSE, heights = NULL, widths = NULL)

pdf("./test.pdf",width = 30,height = 12)
exprs_activity_target_plot$combined_plot
dev.off()
####Circos plot of prioritized ligand-receptor pairs
library(RColorBrewer)
filtered_ligands = ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% top_n(20, prioritization_score) %>% pull(ligand) %>% unique()

prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup() 

colors_sender = brewer.pal(n = prioritized_tbl_oi$sender %>% unique() %>%
                             sort() %>% length(), name = 'Spectral') %>% 
  magrittr::set_names(prioritized_tbl_oi$sender %>% unique() %>% sort())
colors_receiver = c("lavender")  %>% magrittr::set_names(prioritized_tbl_oi$receiver %>% unique() %>% sort())

circos_output = make_circos_lr(prioritized_tbl_oi, colors_sender, colors_receiver)

rm(seurat_obj)
gc()
save.image("./result_TAE.RData")
