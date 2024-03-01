#### Load packages
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
library(harmony)
library(foreach)
library(pheatmap)
library(cowplot)

#### Set up the analysis environment and initialize SCENIC
scRNA=readRDS("scRNA.rds")
exprMat <- as.matrix(scRNA@assays$RNA@counts)
mydbDIR <- "G:\\SCENIC_database/"
mydbs <- c( "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather",
            "hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
names(mydbs) <- c("10kb","500bp")
scenicOptions <- initializeScenic(org="hgnc", 
                                  nCores=14,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "os")
saveRDS(scenicOptions, "int/scenicOptions.rds")
genesKept <- geneFiltering(exprMat, scenicOptions, 
                           minCountsPerGene = 3 * 0.01 * ncol(exprMat), 
                           minSamples = ncol(exprMat) * 0.01)
exprMat_filtered <- exprMat[genesKept, ]

#### Computed correlation matrix
runCorrelation(exprMat_filtered, scenicOptions)

#### TF-Targets correlation regression analysis
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions, nParts = 3)
scenicOptions@settings$nCores <- 1
cellInfo <- data.frame(scRNA@meta.data)
colnames(cellInfo)[which(colnames(cellInfo)=="sampleinf")] <- "sample"
colnames(cellInfo)[which(colnames(cellInfo)=="stim")] <- "stim"
colnames(cellInfo)[which(colnames(cellInfo)=="celltype")] <- "celltype"
cellInfo <- cellInfo[,c("sample","stim","celltype")]
saveRDS(cellInfo, file="int/cellInfo.Rds")
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
exprMat_all <- as.matrix(scRNA@assays$RNA@counts)
exprMat_all <- log2(exprMat_all+1)
runSCENIC_3_scoreCells(scenicOptions, exprMat=exprMat_all)
runSCENIC_4_aucell_binarize(scenicOptions, exprMat=exprMat_all)


#### Import the original regulonAUC matrix
AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
rmove1=colnames(AUCmatrix)[grep("_extended",colnames(AUCmatrix))]
AUCmatrix=AUCmatrix[,!colnames(AUCmatrix)%in%rmove1]
colnames(AUCmatrix)=gsub(" \\(.*\\)$","",colnames(AUCmatrix))
colnames(AUCmatrix)

#### Screening gene expression
TF_dat=cbind(celltype=scRNA@active.ident,AUCmatrix)
TF_dat = TF_dat %>% group_by(celltype) %>% 
  summarise(across(.cols = everything(), mean)) %>% 
  column_to_rownames(.,var = "celltype") %>% t() %>% as.data.frame()

#### Visualize the results using heatmap
pheatmap(TF_dat)
