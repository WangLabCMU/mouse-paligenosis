#### Load packages
suppressMessages(library(ROGUE))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))
library(Seurat)
library(dplyr)

#### Prepare input data
scRNA <- readRDS("./scRNA.rds")
expr <- as.matrix(scRNA@assays$RNA@counts)
meta <- as.data.frame(scRNA@meta.data)
expr <-matr.filter(expr,min.cells=10,min.genes=10)

#### Perform calculations for Rogue
ent.res <-SE_fun(expr)
SEplot(ent.res)
rogue.value <- CalculateRogue(ent.res,platform="UMI")
rogue.res <-rogue(expr,
                  labels=meta$annotation,
                  samples = meta$sampleinf,
                  platform="UMI",
                  span = 0.6)

#### Plot boxplot
rogue.boxplot(rogue.res)+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
