#### Load packages
library(Seurat)
library(ReactomePA)
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(biomaRt)
library(enrichplot)
library(SeuratData)

#### Load data and identify DEG
gogmt<-read.gmt("./c5.go.v2023.2.Hs.symbols.gmt")
DEG=FindMarkers(scRNA,ident.1 = "0",ident.2 = "1",logfc.threshold = 0.25)
geneList = DEG$avg_log2FC
names(geneList) = rownames(DEG)
geneList<-sort(geneList,decreasing = T)

#### Run GSEA
GO<-GSEA(geneList,TERM2GENE = gogmt, minGSSize = 3,maxGSSize = 1000,pvalueCutoff = 1)
GO <- GO[order(GO@result$NES,decreasing = T),]

#### Data processing and visualization
sortdf <- GO
f <- sortdf[sortdf$NES >= 2,]
z <- sortdf[sortdf$NES <= -2,]
sortf <-  f[order(f$NES,decreasing = T),]
sortz <-  z[order(z$NES,decreasing = T),]
sortdf1 <- rbind(sortf,sortz)
sortdf1$Description <- factor(sortdf1$Description, levels = sortdf1$Description,ordered = T)
levels(sortdf1$Description)
sortdf1$group=ifelse(sortdf1$NES>0,"0","1")
table(sortdf1$group)
ggplot(sortdf1, aes(Description, NES, fill=group)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  theme_bw() + 
  theme(panel.grid =element_blank())+
  theme(panel.border = element_rect(size = 0.6)) +
  scale_fill_manual("Group", values = c("0" = "#F0552C", "1" = "#00A0E9"))
