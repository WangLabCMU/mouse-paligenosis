#### Load packages
library(Seurat)
library(GSVA)
library(msigdbr)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(patchwork)
library(pheatmap)

#### Load data
scRNA <- readRDS("scRNA.rds")
expr <- as.matrix(scRNA@assays$RNA@counts)
meta <- scRNA@meta.data[,c("seurat_clusters")]
genesets <- read.gmt("h.all.v7.4.symbols test.gmt")
s.sets <- split(genesets$gene_symbol, genesets$gs_name)

#### Run GSVA
es.matrix = gsva(expr, s.sets, method="gsva", kcdf="Poisson") 

#### Plot the heatmap
pheatmap(es.matrix, 
         show_rownames=1, 
         show_colnames=0, 
         annotation_col=meta,
         fontsize_row=5, 
         width=15, 
         height=12)
