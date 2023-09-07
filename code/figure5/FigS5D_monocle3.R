setwd("/media/ggj/ggjlab2/hezuo/gwh/Tcell/")
library(Seurat)
library(monocle3)
load("./Tcell_harmony.RData")
DimPlot(pbmc.harmony,label=T)
pbmc.harmony <- subset(pbmc.harmony,idents=c(0,1,3,6))
data <- GetAssayData(pbmc.harmony,assay = "RNA",slot = "counts")
cell_metadata <- pbmc.harmony@meta.data
gene_annotation <- data.frame(gene_short_name=rownames(data))
rownames(gene_annotation) <- rownames(data)

cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds,num_dim = 50)

cds <- reduce_dimension(cds,preprocess_method = "PCA")
p1 <- plot_cells(cds,reduction_method = "UMAP",color_cells_by = "type")
p1

cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(pbmc.harmony,reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds,reduction_method = "UMAP",color_cells_by = "type")
p2

cds <- cluster_cells(cds,resolution = 0.0004,random_seed=18,verbose=T)
p1 <- plot_cells(cds)
p1

cds <- learn_graph(cds, verbose =T,use_partition=F,close_loop=F)
p <- plot_cells(cds, color_cells_by = "seurat_clusters", label_groups_by_cluster=T,
                label_leaves=F, label_branch_points=F,cell_size = 0.5,group_label_size=4)

p
get_earliest_principal_node <- function(cds, time_bin="2"){
  cell_ids <- which(cds@clusters@listData[["UMAP"]][["clusters"]] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

cds <- order_cells(cds)
plot_cells(cds,color_cells_by = "pseudotime",label_branch_points = FALSE,label_leaves = F)
ggsave("./Tcell_monocle3.pdf",width = 10,height = 8)
