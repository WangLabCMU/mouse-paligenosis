#### Load packages
library(SLICE)
library(Seurat)

#### Load data and construct object
load(paste("hs_km.Rda", sep=""))
Idents(scRNA) <- scRNA@meta.data$RNA_snn_res.0.6
expr <- GetAssayData(scRNA,slot="data",assay="RNA")
sc<- construct(exprmatrix=as.data.frame(as.matrix(expr)), 
                    cellidentity=scRNA@active.ident)

#### Calculate entropy
sc<- getEntropy(sc, km=km,    
                calculation="bootstrap",   
                B.num=100,         
                exp.cutoff=0.01,       
                B.size=1000,   
                clustering.k=floor(sqrt(1000/2)),  
                random.seed=201602)      

#### Visualization
plotEntropies(sc)
