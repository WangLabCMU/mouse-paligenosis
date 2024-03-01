#### Load packages
library(clustree)
library(cluster)
library(tidyverse)

#### Calculate silhouette coefficient
sce = scRNA
cell_dists <- dist(sce@reductions$harmony@cell.embeddings, method = "euclidean")
cluster_info <- sce@meta.data[, grepl(paste0(DefaultAssay(sce),"_snn_res"), colnames(sce@meta.data))] %>%
  mutate_all(as.character) %>% mutate_all(as.numeric)
silhouette_res <- apply(cluster_info, 2, function(x){
  si <- silhouette(x, cell_dists)
  if(!any(is.na(si))) { mean(si[, 'sil_width']) } else { NA } 
})

#### Plot line chart
si_dat <- data.frame(resolution = names(silhouette_res) %>% str_remove("RNA_snn_res.") %>% as.numeric(),
                     silhouette = silhouette_res)

ggplot(si_dat, aes(x = resolution, y = silhouette)) +
  geom_line() +
  geom_vline(xintercept = si_dat$resolution[which.max(si_dat$silhouette)], linetype = "dashed", color = "red") +
  geom_point(size = 3, shape = 21, fill = "red", color = "black") +
  labs(x = "Resolution", y = "Silhouette Coefficient") +
  ggtitle("Silhouette Coefficient vs Resolution") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_continuous(breaks = seq(0, max(si_dat$resolution), 0.2))+theme_bw()
