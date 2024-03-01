#### Load packages
library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)

#### Construct CDS object
expression_matrix <- as(as.matrix(scRNA@assays$RNA@counts), 'sparseMatrix')
cell_metadata <- scRNA@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix))
rownames(gene_annotation) <- rownames(expression_matrix)
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

#### Preprocessing,dimensionality reduction and clustering
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds,preprocess_method = "PCA") 
cds <- reduce_dimension(cds, reduction_method="tSNE")
cds <- cluster_cells(cds) 

#### Learn the trajectory graph
cds <- learn_graph(cds, use_partition = F)  
plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster=FALSE,
               label_leaves=FALSE, label_branch_points=FALSE,group_label_size=5)
mycds1 <- order_cells(cds)

#### Visualization
plot_cells(mycds1, label_cell_groups = F, 
           color_cells_by = "pseudotime")
