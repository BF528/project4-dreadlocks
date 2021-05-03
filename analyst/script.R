setwd("/projectnb/bf528/users/dreadlocks/project_4")

library(magrittr)
library(rtracklayer)
library(Seurat)
library(tximport)
library(patchwork)
library(tidyverse)

sample.markers <- FindAllMarkers(sample, only.pos = T, min.pct = 0.25, 
                                 logfc.threshold = 0.25)
top_m <- sample.markers[sample.markers$p_val_adj <= 0.05,] %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

sum(sample.markers$p_val_adj <= 0.05)

new.cluster.ids <- c("alpha_1", "alpha_2", "beta", 
                     "ductal", "acinar", "delta",
                     "endothelial_1", "endothelial_2")
names(new.cluster.ids) <- levels(sample)
sample <- RenameIdents(sample, new.cluster.ids)
DimPlot(sample, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# list the top 1
top_1 <- sample.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
g_n <- top_1$gene
VlnPlot(sample, features = g_n)

DoHeatmap(sample, features = top_m$gene)

# save file
write.csv(sample.markers[sample.markers$p_val_adj <= 0.05,],
          file = "analyst/sample_markers.csv",
          row.names = T)