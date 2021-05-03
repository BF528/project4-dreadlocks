setwd("/projectnb/bf528/users/dreadlocks/project_4")

library(magrittr)
library(rtracklayer)
library(Seurat)
library(tximport)
library(patchwork)
library(tidyverse)

# load gtf file
ref_gtf <- import("curator/data/ref/gencode.v37.chr_patch_hapl_scaff.annotation.gtf.gz",
                  format = "GTF")

# load the alevin file path
#/projectnb/bf528/project_4_scrnaseq/GSM2230760__salmon_quant
file <- file.path("curator/result/salmon_count/GSM2230758/alevin/quants_mat.gz")
#file <- file.path("/projectnb/bf528/project_4_scrnaseq/GSM2230760__salmon_quant/alevin/quants_mat.gz")
#file <- file.path("/projectnb/bf528/users/frazzled/project_4/Data_curator/salmon/salmon_output/alevin/quants_mat.gz")
file.exists(file)

# load alevin file
txi <- tximport(file, type = "alevin")
cat("before filtering: ", dim(txi$counts))

# recreate a count matrix 
count_matrix <- txi$counts
gene_symbol_index <- match(rownames(count_matrix),
                           ref_gtf[ref_gtf$type == "gene",]$gene_id)
rownames(count_matrix) <- ref_gtf[ref_gtf$type == "gene",]$gene_name[gene_symbol_index]

sample <- CreateSeuratObject(counts = count_matrix, project = "GSM2230758", 
                             min.cells = 3, min.features = 200)
sample

sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^MT-")
head(sample@meta.data, 20)
VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

CvsM <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "percent.mt")
CvsF <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CvsM + CvsF

sample <- subset(sample, 
                 subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
#sample <- NormalizeData(sample, normalization.method = "LogNormalize", 
#                        scale.factor = 10000)
sample <- NormalizeData(sample)

sample <- FindVariableFeatures(sample, 
                               selection.method = "vst", 
                               nfeatures = 2000)

top10 <- head(VariableFeatures(sample), 10)
VP <- VariableFeaturePlot(sample)
VPt <- LabelPoints(plot = VP, points = top10, repel = TRUE)

all.genes <- rownames(sample)
sample <- ScaleData(sample, features = all.genes)
#sample <- ScaleData(sample, 
#                    vars.to.regress = "percent.mt")

sample <- RunPCA(sample, features = VariableFeatures(object = sample))
VizDimLoadings(sample, dims = 1:2, reduction = "pca")
DimPlot(sample, reduction = "pca")
ElbowPlot(sample)

sample <- FindNeighbors(sample, dims = 1:10)
sample <- FindClusters(sample, resolution = 0.4)
sample <- RunUMAP(sample, dims = 1:10)
DimPlot(sample, reduction = "umap")
barplot(table(sample$seurat_clusters)/sum(table(sample$seurat_clusters)),
        ylim = c(0, 0.5), xlab = "Cluster", ylab = "Proportion")