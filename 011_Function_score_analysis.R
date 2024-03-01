#### Load gene signatures
selected_gene = read.csv("selected gene.csv")
ciliated_genes <- selected_gene$Lysosome

#### calculation function score
scRNA <- scRNA[scRNA@assays$RNA@data@Dimnames[[1]] %in% ciliated_genes,]
colMeans = colMeans(scRNA)
colMeans <- colMeans[rownames(scRNA@meta.data)]
scRNA$colMeans <- colMeans

#### Plot Violin plot
VlnPlot(scRNA, features = "colMeans", pt.size = 0)+geom_boxplot(width=0.07,position=position_dodge(0.9))
ggsave('VlnPlot.pdf', width = 7, height = 5)
