#### Identify differential expression gene
diff.wilcox = FindAllMarkers(scRNA,logfc.threshold = 0.25)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)

#### Plot the top10 gene heatmap
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(object = scRNA,features = top10$gene)

#### Plot the violin plot of specific gene
VlnPlot(scRNA,features = 'PTPRC')

#### Plot the feature plot of specific gene
FeaturePlot(object = scRNA,features = 'PTPRC')
