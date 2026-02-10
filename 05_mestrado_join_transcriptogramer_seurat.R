# Instalar e carregar pacotes necessários
library(Seurat)
library(transcriptogramer)
library(rsparse)
library(Matrix)
library(biomaRt)
library(Ropj)
library(vroom)
library(tidyverse)
library(ggplot2)
library(parallel)
library(ComplexHeatmap)
library(magrittr)
library(dplyr)
library(SingleCellExperiment)
library(scDblFinder)
library(patchwork)
library(SingleR)

# Detectar o número total de núcleos disponíveis no sistema
num_cores <- detectCores()

# Exibir o número de núcleos disponíveis
print(num_cores)

# Definir o número de núcleos globalmente
options(mc.cores = num_cores)

setwd("~/mestrado-single-cell/mestrado-single-cell")

# 1. Pré-processamento dos Dados de Single-Cell com Seurat

# Carregar dados de scRNA-seq (exemplo usando dados integrados no Seurat)

expression_data_notreated_batch1 <- Read10X("mtx_conversions/notreated_batch1/outs/filtered_feature_bc_matrix")
expression_data_notreated_batch2 <- Read10X("mtx_conversions/notreated_batch2/outs/filtered_feature_bc_matrix")
expression_data_TGFbeta1_1day_batch2 <- Read10X("mtx_conversions/TGFbeta1_1day_batch2/outs/filtered_feature_bc_matrix")
expression_data_TGFbeta1_2day_batch2 <- Read10X("mtx_conversions/TGFbeta1_2day_batch2/outs/filtered_feature_bc_matrix")
expression_data_TGFbeta1_3day_batch2 <- Read10X("mtx_conversions/TGFbeta1_3day_batch2/outs/filtered_feature_bc_matrix")
expression_data_TGFbeta1_4day_batch1 <- Read10X("mtx_conversions/TGFbeta1_4day_batch1/outs/filtered_feature_bc_matrix")
expression_data_TGFbeta1_8day_batch1 <- Read10X("mtx_conversions/TGFbeta1_8day_batch1/outs/filtered_feature_bc_matrix")
# data_teste <- Read10X(data.dir = "mtx_conversions/notreated_batch1/outs/filtered_feature_bc_matrix/")

# matrix_notreated_batch1 <- readMM("mtx_conversions/notreated_batch1/outs/filtered_feature_bc_matrix/matrix.mtx")
# matrix_notreated_batch2 <- readMM("mtx_conversions/notreated_batch2/outs/filtered_feature_bc_matrix/matrix.mtx")
# matrix_TGFbeta1_1day_batch2 <- readMM("mtx_conversions/TGFbeta1_1day_batch2/outs/filtered_feature_bc_matrix/matrix.mtx")
# matrix_TGFbeta1_2day_batch2 <- readMM("mtx_conversions/TGFbeta1_2day_batch2/outs/filtered_feature_bc_matrix/matrix.mtx")
# matrix_TGFbeta1_3day_batch2 <- readMM("mtx_conversions/TGFbeta1_3day_batch2/outs/filtered_feature_bc_matrix/matrix.mtx")
# matrix_TGFbeta1_4day_batch1 <- readMM("mtx_conversions/TGFbeta1_4day_batch1/outs/filtered_feature_bc_matrix/matrix.mtx")
# matrix_TGFbeta1_8day_batch1 <- readMM("mtx_conversions/TGFbeta1_8day_batch1/outs/filtered_feature_bc_matrix/matrix.mtx")

# rownames(matrix_notreated_batch1) <- genes_notreated_batch1$V1
# colnames(matrix_notreated_batch1) <- barcodes_notreated_batch1$V1
# 
# rownames(matrix_notreated_batch2) <- genes_notreated_batch2$V1
# colnames(matrix_notreated_batch2) <- barcodes_notreated_batch2$V1
# 
# rownames(matrix_TGFbeta1_1day_batch2) <- genes_TGFbeta1_1day_batch2$V1
# colnames(matrix_TGFbeta1_1day_batch2) <- barcodes_TGFbeta1_1day_batch2$V1
# 
# rownames(matrix_TGFbeta1_2day_batch2) <- genes_TGFbeta1_2day_batch2$V1
# colnames(matrix_TGFbeta1_2day_batch2) <- barcodes_TGFbeta1_2day_batch2$V1
# 
# rownames(matrix_TGFbeta1_3day_batch2) <- genes_TGFbeta1_3day_batch2$V1
# colnames(matrix_TGFbeta1_3day_batch2) <- barcodes_TGFbeta1_3day_batch2$V1
# 
# rownames(matrix_TGFbeta1_4day_batch1) <- genes_TGFbeta1_4day_batch1$V1
# colnames(matrix_TGFbeta1_4day_batch1) <- barcodes_TGFbeta1_4day_batch1$V1
# 
# rownames(matrix_TGFbeta1_8day_batch1) <- genes_TGFbeta1_8day_batch1$V1
# colnames(matrix_TGFbeta1_8day_batch1) <- barcodes_TGFbeta1_8day_batch1$V1


## Forma diversa.




# Criar objeto Seurat
# objetoSeurat_teste <- CreateSeuratObject(counts = data_teste, project = "MCF10A", min.cells = 3, min.features = 200)

# MCF10A_notreated_batch1 <- CreateSeuratObject(counts = matrix_notreated_batch1, project = "MCF10A_notreated_batch1", min.cells = 3, min.features = 200)
# MCF10A_notreated_batch2 <- CreateSeuratObject(counts = matrix_notreated_batch2, project = "MCF10A_notreated_batch2", min.cells = 3, min.features = 200)
# MCF10A_TGFbeta1_1day_batch2 <- CreateSeuratObject(counts = matrix_TGFbeta1_1day_batch2, project = "MCF10A_TGFbeta1_1day_batch2", min.cells = 3, min.features = 200)
# MCF10A_TGFbeta1_2day_batch2 <- CreateSeuratObject(counts = matrix_TGFbeta1_2day_batch2, project = "MCF10A_TGFbeta1_2day_batch2", min.cells = 3, min.features = 200)
# MCF10A_TGFbeta1_3day_batch2 <- CreateSeuratObject(counts = matrix_TGFbeta1_3day_batch2, project = "MCF10A_TGFbeta1_3day_batch2", min.cells = 3, min.features = 200)
# MCF10A_TGFbeta1_4day_batch1 <- CreateSeuratObject(counts = matrix_TGFbeta1_4day_batch1, project = "MCF10A_TGFbeta1_4day_batch1", min.cells = 3, min.features = 200)
# MCF10A_TGFbeta1_8day_batch1 <- CreateSeuratObject(counts = matrix_TGFbeta1_8day_batch1, project = "MCF10A_TGFbeta1_8day_batch1", min.cells = 3, min.features = 200)

MCF10A_notreated_batch1 <- CreateSeuratObject(counts = expression_data_notreated_batch1, project = "MCF10A_notreated_batch1")
MCF10A_notreated_batch2 <- CreateSeuratObject(counts = expression_data_notreated_batch2, project = "MCF10A_notreated_batch2")
MCF10A_TGFbeta1_1day_batch2 <- CreateSeuratObject(counts = expression_data_TGFbeta1_1day_batch2, project = "MCF10A_TGFbeta1_1day_batch2")
MCF10A_TGFbeta1_2day_batch2 <- CreateSeuratObject(counts = expression_data_TGFbeta1_2day_batch2, project = "MCF10A_TGFbeta1_2day_batch2")
MCF10A_TGFbeta1_3day_batch2 <- CreateSeuratObject(counts = expression_data_TGFbeta1_3day_batch2, project = "MCF10A_TGFbeta1_3day_batch2")
MCF10A_TGFbeta1_4day_batch1 <- CreateSeuratObject(counts = expression_data_TGFbeta1_4day_batch1, project = "MCF10A_TGFbeta1_4day_batch1")
MCF10A_TGFbeta1_8day_batch1 <- CreateSeuratObject(counts = expression_data_TGFbeta1_8day_batch1, project = "MCF10A_TGFbeta1_8day_batch1")

# Calcular a porcentagem de genes mitocondriais

MCF10A_notreated_batch1[["percent.mt"]] <- PercentageFeatureSet(MCF10A_notreated_batch1, pattern = "^MT-")  # Para dados humanos (genes mitocondriais começam com "MT-")
MCF10A_notreated_batch2[["percent.mt"]] <- PercentageFeatureSet(MCF10A_notreated_batch2, pattern = "^MT-")
MCF10A_TGFbeta1_1day_batch2[["percent.mt"]] <- PercentageFeatureSet(MCF10A_TGFbeta1_1day_batch2, pattern = "^MT-")
MCF10A_TGFbeta1_2day_batch2[["percent.mt"]] <- PercentageFeatureSet(MCF10A_TGFbeta1_2day_batch2, pattern = "^MT-")
MCF10A_TGFbeta1_3day_batch2[["percent.mt"]] <- PercentageFeatureSet(MCF10A_TGFbeta1_3day_batch2, pattern = "^MT-")
MCF10A_TGFbeta1_4day_batch1[["percent.mt"]] <- PercentageFeatureSet(MCF10A_TGFbeta1_4day_batch1, pattern = "^MT-")
MCF10A_TGFbeta1_8day_batch1[["percent.mt"]] <- PercentageFeatureSet(MCF10A_TGFbeta1_8day_batch1, pattern = "^MT-")

# Violin Plot das métricas de QC
VlnPlot(MCF10A_notreated_batch1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(MCF10A_notreated_batch2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(MCF10A_TGFbeta1_1day_batch2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(MCF10A_TGFbeta1_2day_batch2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(MCF10A_TGFbeta1_3day_batch2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(MCF10A_TGFbeta1_4day_batch1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(MCF10A_TGFbeta1_8day_batch1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Scatter Plots para identificar células problemáticas
FeatureScatter(MCF10A_notreated_batch1, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  FeatureScatter(MCF10A_notreated_batch1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

FeatureScatter(MCF10A_notreated_batch2, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  FeatureScatter(MCF10A_notreated_batch2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

FeatureScatter(MCF10A_TGFbeta1_1day_batch2, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  FeatureScatter(MCF10A_TGFbeta1_1day_batch2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

FeatureScatter(MCF10A_TGFbeta1_2day_batch2, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  FeatureScatter(MCF10A_TGFbeta1_2day_batch2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

FeatureScatter(MCF10A_TGFbeta1_3day_batch2, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  FeatureScatter(MCF10A_TGFbeta1_3day_batch2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

FeatureScatter(MCF10A_TGFbeta1_4day_batch1, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  FeatureScatter(MCF10A_TGFbeta1_4day_batch1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

FeatureScatter(MCF10A_TGFbeta1_8day_batch1, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  FeatureScatter(MCF10A_TGFbeta1_8day_batch1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Filtrar células com base nas métricas de QC
MCF10A_notreated_batch1 <- subset(MCF10A_notreated_batch1, subset = nFeature_RNA > 500 & percent.mt < 20)
MCF10A_notreated_batch2 <- subset(MCF10A_notreated_batch2, subset = nFeature_RNA > 500 & percent.mt < 20)
MCF10A_TGFbeta1_1day_batch2 <- subset(MCF10A_TGFbeta1_1day_batch2, subset = nFeature_RNA > 500 & percent.mt < 20)
MCF10A_TGFbeta1_2day_batch2 <- subset(MCF10A_TGFbeta1_2day_batch2, subset = nFeature_RNA > 500 & percent.mt < 20)
MCF10A_TGFbeta1_3day_batch2 <- subset(MCF10A_TGFbeta1_3day_batch2, subset = nFeature_RNA > 500 & percent.mt < 20)
MCF10A_TGFbeta1_4day_batch1 <- subset(MCF10A_TGFbeta1_4day_batch1, subset = nFeature_RNA > 500 & percent.mt < 20)
MCF10A_TGFbeta1_8day_batch1 <- subset(MCF10A_TGFbeta1_8day_batch1, subset = nFeature_RNA > 500 & percent.mt < 20)

# Normalização dos dados
MCF10A_notreated_batch1 <- NormalizeData(MCF10A_notreated_batch1)
MCF10A_notreated_batch2 <- NormalizeData(MCF10A_notreated_batch2)
MCF10A_TGFbeta1_1day_batch2 <- NormalizeData(MCF10A_TGFbeta1_1day_batch2)
MCF10A_TGFbeta1_2day_batch2 <- NormalizeData(MCF10A_TGFbeta1_2day_batch2)
MCF10A_TGFbeta1_3day_batch2 <- NormalizeData(MCF10A_TGFbeta1_3day_batch2)
MCF10A_TGFbeta1_4day_batch1 <- NormalizeData(MCF10A_TGFbeta1_4day_batch1)
MCF10A_TGFbeta1_8day_batch1 <- NormalizeData(MCF10A_TGFbeta1_8day_batch1)


# Identificar genes variáveis
MCF10A_notreated_batch1 <- FindVariableFeatures(MCF10A_notreated_batch1, selection.method = "vst", nfeatures = 2000)
MCF10A_notreated_batch2 <- FindVariableFeatures(MCF10A_notreated_batch2, selection.method = "vst", nfeatures = 2000)
MCF10A_TGFbeta1_1day_batch2 <- FindVariableFeatures(MCF10A_TGFbeta1_1day_batch2, selection.method = "vst", nfeatures = 2000)
MCF10A_TGFbeta1_2day_batch2 <- FindVariableFeatures(MCF10A_TGFbeta1_2day_batch2, selection.method = "vst", nfeatures = 2000)
MCF10A_TGFbeta1_3day_batch2 <- FindVariableFeatures(MCF10A_TGFbeta1_3day_batch2, selection.method = "vst", nfeatures = 2000)
MCF10A_TGFbeta1_4day_batch1 <- FindVariableFeatures(MCF10A_TGFbeta1_4day_batch1, selection.method = "vst", nfeatures = 2000)
MCF10A_TGFbeta1_8day_batch1 <- FindVariableFeatures(MCF10A_TGFbeta1_8day_batch1, selection.method = "vst", nfeatures = 2000)


# Escalonar os dados
MCF10A_notreated_batch1 <- ScaleData(MCF10A_notreated_batch1, features = rownames(MCF10A_notreated_batch1))
MCF10A_notreated_batch2 <- ScaleData(MCF10A_notreated_batch2, features = rownames(MCF10A_notreated_batch2))
MCF10A_TGFbeta1_1day_batch2 <- ScaleData(MCF10A_TGFbeta1_1day_batch2, features = rownames(MCF10A_TGFbeta1_1day_batch2))
MCF10A_TGFbeta1_2day_batch2 <- ScaleData(MCF10A_TGFbeta1_2day_batch2, features = rownames(MCF10A_TGFbeta1_2day_batch2))
MCF10A_TGFbeta1_3day_batch2 <- ScaleData(MCF10A_TGFbeta1_3day_batch2, features = rownames(MCF10A_TGFbeta1_3day_batch2))
MCF10A_TGFbeta1_4day_batch1 <- ScaleData(MCF10A_TGFbeta1_4day_batch1, features = rownames(MCF10A_TGFbeta1_4day_batch1))
MCF10A_TGFbeta1_8day_batch1 <- ScaleData(MCF10A_TGFbeta1_8day_batch1, features = rownames(MCF10A_TGFbeta1_8day_batch1))

# Executar PCA
MCF10A_notreated_batch1 <- RunPCA(MCF10A_notreated_batch1, features = VariableFeatures(object = MCF10A_notreated_batch1))
MCF10A_notreated_batch2 <- RunPCA(MCF10A_notreated_batch2, features = VariableFeatures(object = MCF10A_notreated_batch2))
MCF10A_TGFbeta1_1day_batch2 <- RunPCA(MCF10A_TGFbeta1_1day_batch2, features = VariableFeatures(object = MCF10A_TGFbeta1_1day_batch2))
MCF10A_TGFbeta1_2day_batch2 <- RunPCA(MCF10A_TGFbeta1_2day_batch2, features = VariableFeatures(object = MCF10A_TGFbeta1_2day_batch2))
MCF10A_TGFbeta1_3day_batch2 <- RunPCA(MCF10A_TGFbeta1_3day_batch2, features = VariableFeatures(object = MCF10A_TGFbeta1_3day_batch2))
MCF10A_TGFbeta1_4day_batch1 <- RunPCA(MCF10A_TGFbeta1_4day_batch1, features = VariableFeatures(object = MCF10A_TGFbeta1_4day_batch1))
MCF10A_TGFbeta1_8day_batch1 <- RunPCA(MCF10A_TGFbeta1_8day_batch1, features = VariableFeatures(object = MCF10A_TGFbeta1_8day_batch1))

# Visualizar os resultados do PCA
ElbowPlot(MCF10A_notreated_batch1)
ElbowPlot(MCF10A_notreated_batch2)
ElbowPlot(MCF10A_TGFbeta1_1day_batch2)
ElbowPlot(MCF10A_TGFbeta1_2day_batch2)
ElbowPlot(MCF10A_TGFbeta1_3day_batch2)
ElbowPlot(MCF10A_TGFbeta1_4day_batch1)
ElbowPlot(MCF10A_TGFbeta1_8day_batch1)

# Clustering
MCF10A_notreated_batch1 <- FindNeighbors(MCF10A_notreated_batch1, dims = 1:10)
MCF10A_notreated_batch2 <- FindNeighbors(MCF10A_notreated_batch2, dims = 1:10)
MCF10A_TGFbeta1_1day_batch2 <- FindNeighbors(MCF10A_TGFbeta1_1day_batch2, dims = 1:10)
MCF10A_TGFbeta1_2day_batch2 <- FindNeighbors(MCF10A_TGFbeta1_2day_batch2, dims = 1:10)
MCF10A_TGFbeta1_3day_batch2 <- FindNeighbors(MCF10A_TGFbeta1_3day_batch2, dims = 1:10)
MCF10A_TGFbeta1_4day_batch1 <- FindNeighbors(MCF10A_TGFbeta1_4day_batch1, dims = 1:10)
MCF10A_TGFbeta1_8day_batch1 <- FindNeighbors(MCF10A_TGFbeta1_8day_batch1, dims = 1:10)


MCF10A_notreated_batch1 <- FindClusters(MCF10A_notreated_batch1, resolution = 0.5)
MCF10A_notreated_batch2 <- FindClusters(MCF10A_notreated_batch2, resolution = 0.5)
MCF10A_TGFbeta1_1day_batch2 <- FindClusters(MCF10A_TGFbeta1_1day_batch2, resolution = 0.5)
MCF10A_TGFbeta1_2day_batch2 <- FindClusters(MCF10A_TGFbeta1_2day_batch2, resolution = 0.5)
MCF10A_TGFbeta1_3day_batch2 <- FindClusters(MCF10A_TGFbeta1_3day_batch2, resolution = 0.5)
MCF10A_TGFbeta1_4day_batch1 <- FindClusters(MCF10A_TGFbeta1_4day_batch1, resolution = 0.5)
MCF10A_TGFbeta1_8day_batch1 <- FindClusters(MCF10A_TGFbeta1_8day_batch1, resolution = 0.5)


# UMAP para visualização
MCF10A_notreated_batch1 <- RunUMAP(MCF10A_notreated_batch1, dims = 1:10)
MCF10A_notreated_batch2 <- RunUMAP(MCF10A_notreated_batch2, dims = 1:10)
MCF10A_TGFbeta1_1day_batch2 <- RunUMAP(MCF10A_TGFbeta1_1day_batch2, dims = 1:10)
MCF10A_TGFbeta1_2day_batch2 <- RunUMAP(MCF10A_TGFbeta1_2day_batch2, dims = 1:10)
MCF10A_TGFbeta1_3day_batch2 <- RunUMAP(MCF10A_TGFbeta1_3day_batch2, dims = 1:10)
MCF10A_TGFbeta1_4day_batch1 <- RunUMAP(MCF10A_TGFbeta1_4day_batch1, dims = 1:10)
MCF10A_TGFbeta1_8day_batch1 <- RunUMAP(MCF10A_TGFbeta1_8day_batch1, dims = 1:10)


DimPlot(MCF10A_notreated_batch1, reduction = "umap")
DimPlot(MCF10A_notreated_batch1, reduction = "pca")

DimPlot(MCF10A_notreated_batch2, reduction = "umap")
DimPlot(MCF10A_notreated_batch2, reduction = "pca")

DimPlot(MCF10A_TGFbeta1_1day_batch2, reduction = "umap")
DimPlot(MCF10A_TGFbeta1_1day_batch2, reduction = "pca")

DimPlot(MCF10A_TGFbeta1_2day_batch2, reduction = "umap")
DimPlot(MCF10A_TGFbeta1_2day_batch2, reduction = "pca")

DimPlot(MCF10A_TGFbeta1_3day_batch2, reduction = "umap")
DimPlot(MCF10A_TGFbeta1_3day_batch2, reduction = "pca")

DimPlot(MCF10A_TGFbeta1_4day_batch1, reduction = "umap")
DimPlot(MCF10A_TGFbeta1_4day_batch1, reduction = "pca")

DimPlot(MCF10A_TGFbeta1_8day_batch1, reduction = "umap")
DimPlot(MCF10A_TGFbeta1_8day_batch1, reduction = "pca")

# 2. Extração da Matriz de Expressão

# Extrair matriz de expressão
expression_matrix_MCF10A_notreated_batch1 <- as.data.frame(GetAssayData(MCF10A_notreated_batch1, layer = "data", assay = "RNA"))
expression_matrix_MCF10A_notreated_batch2 <- as.data.frame(GetAssayData(MCF10A_notreated_batch2, layer = "data", assay = "RNA"))
expression_matrix_MCF10A_TGFbeta1_1day_batch2 <- as.data.frame(GetAssayData(MCF10A_TGFbeta1_1day_batch2, layer = "data", assay = "RNA"))
expression_matrix_MCF10A_TGFbeta1_2day_batch2 <- as.data.frame(GetAssayData(MCF10A_TGFbeta1_2day_batch2, layer = "data", assay = "RNA"))
expression_matrix_MCF10A_TGFbeta1_3day_batch2 <- as.data.frame(GetAssayData(MCF10A_TGFbeta1_3day_batch2, layer = "data", assay = "RNA"))
expression_matrix_MCF10A_TGFbeta1_4day_batch1 <- as.data.frame(GetAssayData(MCF10A_TGFbeta1_4day_batch1, layer = "data", assay = "RNA"))
expression_matrix_MCF10A_TGFbeta1_8day_batch1 <- as.data.frame(GetAssayData(MCF10A_TGFbeta1_8day_batch1, layer = "data", assay = "RNA"))


genes_notreated_batch1 <- read.table("mtx_conversions/notreated_batch1/outs/filtered_feature_bc_matrix/genes.tsv", header = F)
genes_notreated_batch2 <- read.table("mtx_conversions/notreated_batch2/outs/filtered_feature_bc_matrix/genes.tsv", header = F)
genes_TGFbeta1_1day_batch2 <- read.table("mtx_conversions/TGFbeta1_1day_batch2/outs/filtered_feature_bc_matrix/genes.tsv", header = F)
genes_TGFbeta1_2day_batch2 <- read.table("mtx_conversions/TGFbeta1_2day_batch2/outs/filtered_feature_bc_matrix/genes.tsv", header = F)
genes_TGFbeta1_3day_batch2 <- read.table("mtx_conversions/TGFbeta1_3day_batch2/outs/filtered_feature_bc_matrix/genes.tsv", header = F)
genes_TGFbeta1_4day_batch1 <- read.table("mtx_conversions/TGFbeta1_4day_batch1/outs/filtered_feature_bc_matrix/genes.tsv", header = F)
genes_TGFbeta1_8day_batch1 <- read.table("mtx_conversions/TGFbeta1_8day_batch1/outs/filtered_feature_bc_matrix/genes.tsv", header = F)


barcodes_notreated_batch1 <- read.table("mtx_conversions/notreated_batch1/outs/filtered_feature_bc_matrix/barcodes.tsv", header = F)
barcodes_notreated_batch2 <- read.table("mtx_conversions/notreated_batch2/outs/filtered_feature_bc_matrix/barcodes.tsv", header = F)
barcodes_TGFbeta1_1day_batch2 <- read.table("mtx_conversions/TGFbeta1_1day_batch2/outs/filtered_feature_bc_matrix/barcodes.tsv", header = F)
barcodes_TGFbeta1_2day_batch2 <- read.table("mtx_conversions/TGFbeta1_2day_batch2/outs/filtered_feature_bc_matrix/barcodes.tsv", header = F)
barcodes_TGFbeta1_3day_batch2 <- read.table("mtx_conversions/TGFbeta1_3day_batch2/outs/filtered_feature_bc_matrix/barcodes.tsv", header = F)
barcodes_TGFbeta1_4day_batch1 <- read.table("mtx_conversions/TGFbeta1_4day_batch1/outs/filtered_feature_bc_matrix/barcodes.tsv", header = F)
barcodes_TGFbeta1_8day_batch1 <- read.table("mtx_conversions/TGFbeta1_8day_batch1/outs/filtered_feature_bc_matrix/barcodes.tsv", header = F)

# Verificar se os nomes dos genes no objeto `expression_matrix_MCF10A_notreated_batch1` correspondem aos símbolos do mapeamento
current_genes_expression_matrix_MCF10A_notreated_batch1 <- 
  rownames(expression_matrix_MCF10A_notreated_batch1)
# Substituir os símbolos pelos ENSG usando o mapeamento
new_gene_names_expression_matrix_MCF10A_notreated_batch1 <- 
  genes_notreated_batch1$V1[match(current_genes_expression_matrix_MCF10A_notreated_batch1,
                                  genes_notreated_batch1$V2)]
# Substituir os nomes das linhas da matriz de expressão
rownames(expression_matrix_MCF10A_notreated_batch1) <- 
  new_gene_names_expression_matrix_MCF10A_notreated_batch1

expression_matrix_MCF10A_notreated_batch1$RowNames <- rownames(expression_matrix_MCF10A_notreated_batch1)
expression_matrix_MCF10A_notreated_batch2$RowNames <- rownames(expression_matrix_MCF10A_notreated_batch2)
expression_matrix_MCF10A_TGFbeta1_1day_batch2$RowNames <- rownames(expression_matrix_MCF10A_TGFbeta1_1day_batch2)
expression_matrix_MCF10A_TGFbeta1_2day_batch2$RowNames <- rownames(expression_matrix_MCF10A_TGFbeta1_2day_batch2)
expression_matrix_MCF10A_TGFbeta1_3day_batch2$RowNames <- rownames(expression_matrix_MCF10A_TGFbeta1_3day_batch2)
expression_matrix_MCF10A_TGFbeta1_4day_batch1$RowNames <- rownames(expression_matrix_MCF10A_TGFbeta1_4day_batch1)
expression_matrix_MCF10A_TGFbeta1_8day_batch1$RowNames <- rownames(expression_matrix_MCF10A_TGFbeta1_8day_batch1)

expression_matrix_MCF10A <- full_join(expression_matrix_MCF10A_notreated_batch1, expression_matrix_MCF10A_notreated_batch2, by = "RowNames") %>%
  full_join(expression_matrix_MCF10A_TGFbeta1_1day_batch2, by = "RowNames") %>%
  full_join(expression_matrix_MCF10A_TGFbeta1_2day_batch2, by = "RowNames") %>%
  full_join(expression_matrix_MCF10A_TGFbeta1_3day_batch2, by = "RowNames") %>%
  full_join(expression_matrix_MCF10A_TGFbeta1_4day_batch1, by = "RowNames") %>%
  full_join(expression_matrix_MCF10A_TGFbeta1_8day_batch1, by = "RowNames")

# Preenchendo os NAs com 0
expression_matrix_MCF10A[is.na(expression_matrix_MCF10A)] <- 0

# Movendo os nomes das linhas de volta para rownames e removendo a coluna auxiliar
rownames(expression_matrix_MCF10A) <- expression_matrix_MCF10A$RowNames
expression_matrix_MCF10A$RowNames <- NULL
expression_matrix_MCF10A_notreated_batch1$RowNames <- NULL
expression_matrix_MCF10A_notreated_batch2$RowNames <- NULL
expression_matrix_MCF10A_TGFbeta1_1day_batch2$RowNames <- NULL
expression_matrix_MCF10A_TGFbeta1_2day_batch2$RowNames <- NULL
expression_matrix_MCF10A_TGFbeta1_3day_batch2$RowNames <- NULL
expression_matrix_MCF10A_TGFbeta1_4day_batch1$RowNames <- NULL
expression_matrix_MCF10A_TGFbeta1_8day_batch1$RowNames <- NULL


dim(expression_matrix_MCF10A_notreated_batch1)
dim(expression_matrix_MCF10A_notreated_batch2)
dim(expression_matrix_MCF10A_TGFbeta1_1day_batch2)
dim(expression_matrix_MCF10A_TGFbeta1_2day_batch2)
dim(expression_matrix_MCF10A_TGFbeta1_3day_batch2)
dim(expression_matrix_MCF10A_TGFbeta1_4day_batch1)
dim(expression_matrix_MCF10A_TGFbeta1_8day_batch1)
dim(expression_matrix_MCF10A)

sum(dim(expression_matrix_MCF10A_notreated_batch1)[2],
    dim(expression_matrix_MCF10A_notreated_batch2)[2],
    dim(expression_matrix_MCF10A_TGFbeta1_1day_batch2)[2],
    dim(expression_matrix_MCF10A_TGFbeta1_2day_batch2)[2],
    dim(expression_matrix_MCF10A_TGFbeta1_3day_batch2)[2],
    dim(expression_matrix_MCF10A_TGFbeta1_4day_batch1)[2],
    dim(expression_matrix_MCF10A_TGFbeta1_8day_batch1)[2])

# Criando matrizes de média e desvio padrão.

mean_expression_matrix_MCF10A_notreated_batch1 <- 
  rowMeans(expression_matrix_MCF10A_notreated_batch1[sapply(expression_matrix_MCF10A_notreated_batch1, is.numeric)])
sd_expression_matrix_MCF10A_notreated_batch1 <- 
  apply(expression_matrix_MCF10A_notreated_batch1[sapply(expression_matrix_MCF10A_notreated_batch1, is.numeric)], 1, sd)
mean_sd_expression_matrix_MCF10A_notreated_batch1 <- as.data.frame(cbind(MeanNotreatedBatch1 = mean_expression_matrix_MCF10A_notreated_batch1, 
                                                           SD_NotreatedBatch1= sd_expression_matrix_MCF10A_notreated_batch1))


mean_expression_matrix_MCF10A_notreated_batch2 <- 
  rowMeans(expression_matrix_MCF10A_notreated_batch2[sapply(expression_matrix_MCF10A_notreated_batch2, is.numeric)])
sd_expression_matrix_MCF10A_notreated_batch2 <- 
  apply(expression_matrix_MCF10A_notreated_batch2[sapply(expression_matrix_MCF10A_notreated_batch2, is.numeric)], 1, sd)
mean_sd_expression_matrix_MCF10A_notreated_batch2 <- as.data.frame(cbind(MeanNotreatedBatch2 = mean_expression_matrix_MCF10A_notreated_batch2, 
                                                           SD_NotreatedBatch2= sd_expression_matrix_MCF10A_notreated_batch2))

mean_expression_matrix_MCF10A_TGFbeta1_1day_batch2 <- 
  rowMeans(expression_matrix_MCF10A_TGFbeta1_1day_batch2[sapply(expression_matrix_MCF10A_TGFbeta1_1day_batch2, is.numeric)])
sd_expression_matrix_MCF10A_TGFbeta1_1day_batch2 <- 
  apply(expression_matrix_MCF10A_TGFbeta1_1day_batch2[sapply(expression_matrix_MCF10A_TGFbeta1_1day_batch2, is.numeric)], 1, sd)
mean_sd_expression_matrix_MCF10A_TGFbeta1_1day_batch2 <- as.data.frame(cbind(MeanTGFbeta1_1dayBatch2 = mean_expression_matrix_MCF10A_TGFbeta1_1day_batch2, 
                                                           SD_TGFbeta1_1dayBatch2= sd_expression_matrix_MCF10A_TGFbeta1_1day_batch2))

mean_expression_matrix_MCF10A_TGFbeta1_2day_batch2 <- 
  rowMeans(expression_matrix_MCF10A_TGFbeta1_2day_batch2[sapply(expression_matrix_MCF10A_TGFbeta1_2day_batch2, is.numeric)])
sd_expression_matrix_MCF10A_TGFbeta1_2day_batch2 <- 
  apply(expression_matrix_MCF10A_TGFbeta1_2day_batch2[sapply(expression_matrix_MCF10A_TGFbeta1_2day_batch2, is.numeric)], 1, sd)
mean_sd_expression_matrix_MCF10A_TGFbeta1_2day_batch2 <- as.data.frame(cbind(MeanTGFbeta1_2dayBatch2 = mean_expression_matrix_MCF10A_TGFbeta1_2day_batch2, 
                                                               SD_TGFbeta1_2dayBatch2= sd_expression_matrix_MCF10A_TGFbeta1_2day_batch2))

mean_expression_matrix_MCF10A_TGFbeta1_3day_batch2 <- 
  rowMeans(expression_matrix_MCF10A_TGFbeta1_3day_batch2[sapply(expression_matrix_MCF10A_TGFbeta1_3day_batch2, is.numeric)])
sd_expression_matrix_MCF10A_TGFbeta1_3day_batch2 <- 
  apply(expression_matrix_MCF10A_TGFbeta1_3day_batch2[sapply(expression_matrix_MCF10A_TGFbeta1_3day_batch2, is.numeric)], 1, sd)
mean_sd_expression_matrix_MCF10A_TGFbeta1_3day_batch2 <- as.data.frame(cbind(MeanTGFbeta1_3dayBatch2 = mean_expression_matrix_MCF10A_TGFbeta1_3day_batch2, 
                                                               SD_TGFbeta1_3dayBatch2= sd_expression_matrix_MCF10A_TGFbeta1_3day_batch2))

mean_expression_matrix_MCF10A_TGFbeta1_4day_batch1 <- 
  rowMeans(expression_matrix_MCF10A_TGFbeta1_4day_batch1[sapply(expression_matrix_MCF10A_TGFbeta1_4day_batch1, is.numeric)])
sd_expression_matrix_MCF10A_TGFbeta1_4day_batch1 <- 
  apply(expression_matrix_MCF10A_TGFbeta1_4day_batch1[sapply(expression_matrix_MCF10A_TGFbeta1_4day_batch1, is.numeric)], 1, sd)
mean_sd_expression_matrix_MCF10A_TGFbeta1_4day_batch1 <- as.data.frame(cbind(MeanTGFbeta1_4dayBatch1 = mean_expression_matrix_MCF10A_TGFbeta1_4day_batch1, 
                                                               SD_TGFbeta1_4dayBatch1= sd_expression_matrix_MCF10A_TGFbeta1_4day_batch1))

mean_expression_matrix_MCF10A_TGFbeta1_8day_batch1 <- 
  rowMeans(expression_matrix_MCF10A_TGFbeta1_8day_batch1[sapply(expression_matrix_MCF10A_TGFbeta1_8day_batch1, is.numeric)])
sd_expression_matrix_MCF10A_TGFbeta1_8day_batch1 <- 
  apply(expression_matrix_MCF10A_TGFbeta1_8day_batch1[sapply(expression_matrix_MCF10A_TGFbeta1_8day_batch1, is.numeric)], 1, sd)
mean_sd_expression_matrix_MCF10A_TGFbeta1_8day_batch1 <- as.data.frame(cbind(MeanTGFbeta1_8dayBatch1 = mean_expression_matrix_MCF10A_TGFbeta1_8day_batch1, 
                                                               SD_TGFbeta1_8dayBatch1= sd_expression_matrix_MCF10A_TGFbeta1_8day_batch1))

# Juntando as matrizes de médias e desvios.

mean_sd_expression_matrix_MCF10A_notreated_batch1$RowNames <- rownames(mean_sd_expression_matrix_MCF10A_notreated_batch1)
mean_sd_expression_matrix_MCF10A_notreated_batch2$RowNames <- rownames(mean_sd_expression_matrix_MCF10A_notreated_batch2)
mean_sd_expression_matrix_MCF10A_TGFbeta1_1day_batch2$RowNames <- rownames(mean_sd_expression_matrix_MCF10A_TGFbeta1_1day_batch2)
mean_sd_expression_matrix_MCF10A_TGFbeta1_2day_batch2$RowNames <- rownames(mean_sd_expression_matrix_MCF10A_TGFbeta1_2day_batch2)
mean_sd_expression_matrix_MCF10A_TGFbeta1_3day_batch2$RowNames <- rownames(mean_sd_expression_matrix_MCF10A_TGFbeta1_3day_batch2)
mean_sd_expression_matrix_MCF10A_TGFbeta1_4day_batch1$RowNames <- rownames(mean_sd_expression_matrix_MCF10A_TGFbeta1_4day_batch1)
mean_sd_expression_matrix_MCF10A_TGFbeta1_8day_batch1$RowNames <- rownames(mean_sd_expression_matrix_MCF10A_TGFbeta1_8day_batch1)

mean_sd_expression_matrix_MCF10A <- full_join(mean_sd_expression_matrix_MCF10A_notreated_batch1, 
                                              mean_sd_expression_matrix_MCF10A_notreated_batch2, by = "RowNames") %>%
  full_join(mean_sd_expression_matrix_MCF10A_TGFbeta1_1day_batch2, by = "RowNames") %>%
  full_join(mean_sd_expression_matrix_MCF10A_TGFbeta1_2day_batch2, by = "RowNames") %>%
  full_join(mean_sd_expression_matrix_MCF10A_TGFbeta1_3day_batch2, by = "RowNames") %>%
  full_join(mean_sd_expression_matrix_MCF10A_TGFbeta1_4day_batch1, by = "RowNames") %>%
  full_join(mean_sd_expression_matrix_MCF10A_TGFbeta1_8day_batch1, by = "RowNames")


# Preenchendo os NAs com 0
mean_sd_expression_matrix_MCF10A[is.na(mean_sd_expression_matrix_MCF10A)] <- 0

# Movendo os nomes das linhas de volta para rownames e removendo a coluna auxiliar
rownames(mean_sd_expression_matrix_MCF10A) <- mean_sd_expression_matrix_MCF10A$RowNames
mean_sd_expression_matrix_MCF10A$RowNames <- NULL
mean_sd_expression_matrix_MCF10A_notreated_batch1$RowNames <- NULL
mean_sd_expression_matrix_MCF10A_notreated_batch2$RowNames <- NULL
mean_sd_expression_matrix_MCF10A_TGFbeta1_1day_batch2$RowNames <- NULL
mean_sd_expression_matrix_MCF10A_TGFbeta1_2day_batch2$RowNames <- NULL
mean_sd_expression_matrix_MCF10A_TGFbeta1_3day_batch2$RowNames <- NULL
mean_sd_expression_matrix_MCF10A_TGFbeta1_4day_batch1$RowNames <- NULL
mean_sd_expression_matrix_MCF10A_TGFbeta1_8day_batch1$RowNames <- NULL



# Import expression matrix
# exp <- read.csv("expressionMatrixDayZero.csv")
# rownames(exp) <- exp$X
# exp$X <- NULL

# See transcriptogramerStep1 documentation
# First column of the dictionary MUST BE THE ENSEMBL PEPTIDE ID
# SECOND COLUMN MUST BE THE SAME ID AS ON THE GENE EXPRESSION MATRIX ROWNAMES
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
dictionary <- getBM(attributes = c("ensembl_peptide_id", "ensembl_gene_id"),
              mart = ensembl)
dictionary %>% 
  mutate(ensembl_peptide_id = ifelse(ensembl_peptide_id == "", NA, ensembl_peptide_id)) %>% 
  na.omit() -> dictionary

# Ordering
ord <- vroom("ordering_HomoSapiensScore800-2024-C.txt")

# Association matrix
assoc <- read.opj("Associationmatrix.opj")
assoc <- assoc$associationMa

inner_join(assoc, ord, by = c("A" = "dim1")) %>% 
  dplyr::rename(protein1 = Protein) %>% 
  dplyr::inner_join(ord, by = c("B" = "dim1")) %>% 
  dplyr::rename(protein2 = Protein) %>% 
  dplyr::select(protein1, protein2) -> assoc



# Run transcriptogram preprocess, radius = 30
t_notreated_batch1 <- transcriptogramPreprocess(association = assoc, ordering = ord$Protein, radius = 30)
t_notreated_batch1 <- transcriptogramStep1(object = t_notreated_batch1, expression = expression_matrix_MCF10A_notreated_batch1, dictionary = dictionary, nCores = 40)
t_notreated_batch1 <- transcriptogramStep2(object = t_notreated_batch1, nCores = 40)

t_notreated_batch2 <- transcriptogramPreprocess(association = assoc, ordering = ord$Protein, radius = 30)
t_notreated_batch2 <- transcriptogramStep1(object = t_notreated_batch2, expression = expression_matrix_MCF10A_notreated_batch2, dictionary = dictionary, nCores = 40)
t_notreated_batch2 <- transcriptogramStep2(object = t_notreated_batch2, nCores = 40)

t_TGFbeta1_1day_batch2 <- transcriptogramPreprocess(association = assoc, ordering = ord$Protein, radius = 30)
t_TGFbeta1_1day_batch2 <- transcriptogramStep1(object = t_TGFbeta1_1day_batch2, expression = expression_matrix_MCF10A_TGFbeta1_1day_batch2, dictionary = dictionary, nCores = 40)
t_TGFbeta1_1day_batch2 <- transcriptogramStep2(object = t_TGFbeta1_1day_batch2, nCores = 40)

t_TGFbeta1_2day_batch2 <- transcriptogramPreprocess(association = assoc, ordering = ord$Protein, radius = 30)
t_TGFbeta1_2day_batch2 <- transcriptogramStep1(object = t_TGFbeta1_2day_batch2, expression = expression_matrix_MCF10A_TGFbeta1_2day_batch2, dictionary = dictionary, nCores = 40)
t_TGFbeta1_2day_batch2 <- transcriptogramStep2(object = t_TGFbeta1_2day_batch2, nCores = 40)

t_TGFbeta1_3day_batch2 <- transcriptogramPreprocess(association = assoc, ordering = ord$Protein, radius = 30)
t_TGFbeta1_3day_batch2 <- transcriptogramStep1(object = t_TGFbeta1_3day_batch2, expression = expression_matrix_MCF10A_TGFbeta1_3day_batch2, dictionary = dictionary, nCores = 40)
t_TGFbeta1_3day_batch2 <- transcriptogramStep2(object = t_TGFbeta1_3day_batch2, nCores = 40)

t_TGFbeta1_4day_batch1 <- transcriptogramPreprocess(association = assoc, ordering = ord$Protein, radius = 30)
t_TGFbeta1_4day_batch1 <- transcriptogramStep1(object = t_TGFbeta1_4day_batch1, expression = expression_matrix_MCF10A_TGFbeta1_4day_batch1, dictionary = dictionary, nCores = 40)
t_TGFbeta1_4day_batch1 <- transcriptogramStep2(object = t_TGFbeta1_4day_batch1, nCores = 40)

t_TGFbeta1_8day_batch1 <- transcriptogramPreprocess(association = assoc, ordering = ord$Protein, radius = 30)
t_TGFbeta1_8day_batch1 <- transcriptogramStep1(object = t_TGFbeta1_8day_batch1, expression = expression_matrix_MCF10A_TGFbeta1_8day_batch1, dictionary = dictionary, nCores = 10)
t_TGFbeta1_8day_batch1 <- transcriptogramStep2(object = t_TGFbeta1_8day_batch1, nCores = 10)

save(t_notreated_batch1,
     t_notreated_batch2,
     t_TGFbeta1_1day_batch2,
     t_TGFbeta1_2day_batch2,
     t_TGFbeta1_3day_batch2,
     t_TGFbeta1_4day_batch1,
     t_TGFbeta1_8day_batch1, file = 
     "SomenteTranscriptogramas.RData")
# load("AposTranscriptogramasDiarios.RData")
# load("SomenteTranscriptogramas.RData")

# expression_matrix_MCF10A using colnames from days
expression_matrix_MCF10A_notreated_batch1_forDay <- expression_matrix_MCF10A_notreated_batch1
names_notreated_batch1 <- paste0("notreated_batch1_", seq_len(ncol(expression_matrix_MCF10A_notreated_batch1_forDay)))
colnames(expression_matrix_MCF10A_notreated_batch1_forDay) <- names_notreated_batch1

expression_matrix_MCF10A_notreated_batch2_forDay <- expression_matrix_MCF10A_notreated_batch2
names_notreated_batch2 <- paste0("notreated_batch2_", seq_len(ncol(expression_matrix_MCF10A_notreated_batch2_forDay)))
colnames(expression_matrix_MCF10A_notreated_batch2_forDay) <- names_notreated_batch2

expression_matrix_MCF10A_TGFbeta1_1day_batch2_forDay <- expression_matrix_MCF10A_TGFbeta1_1day_batch2
names_TGFbeta1_1day_batch2 <- paste0("TGFbeta1_1day_batch2_", seq_len(ncol(expression_matrix_MCF10A_TGFbeta1_1day_batch2_forDay)))
colnames(expression_matrix_MCF10A_TGFbeta1_1day_batch2_forDay) <- names_TGFbeta1_1day_batch2

expression_matrix_MCF10A_TGFbeta1_2day_batch2_forDay <- expression_matrix_MCF10A_TGFbeta1_2day_batch2
names_TGFbeta1_2day_batch2 <- paste0("TGFbeta1_2day_batch2_", seq_len(ncol(expression_matrix_MCF10A_TGFbeta1_2day_batch2_forDay)))
colnames(expression_matrix_MCF10A_TGFbeta1_2day_batch2_forDay) <- names_TGFbeta1_2day_batch2

expression_matrix_MCF10A_TGFbeta1_3day_batch2_forDay <- expression_matrix_MCF10A_TGFbeta1_3day_batch2
names_TGFbeta1_3day_batch2 <- paste0("TGFbeta1_3day_batch2_", seq_len(ncol(expression_matrix_MCF10A_TGFbeta1_3day_batch2_forDay)))
colnames(expression_matrix_MCF10A_TGFbeta1_3day_batch2_forDay) <- names_TGFbeta1_3day_batch2

expression_matrix_MCF10A_TGFbeta1_4day_batch1_forDay <- expression_matrix_MCF10A_TGFbeta1_4day_batch1
names_TGFbeta1_4day_batch1 <- paste0("TGFbeta1_4day_batch1_", seq_len(ncol(expression_matrix_MCF10A_TGFbeta1_4day_batch1_forDay)))
colnames(expression_matrix_MCF10A_TGFbeta1_4day_batch1_forDay) <- names_TGFbeta1_4day_batch1

expression_matrix_MCF10A_TGFbeta1_8day_batch1_forDay <- expression_matrix_MCF10A_TGFbeta1_8day_batch1
names_TGFbeta1_8day_batch1 <- paste0("TGFbeta1_8day_batch1_", seq_len(ncol(expression_matrix_MCF10A_TGFbeta1_8day_batch1_forDay)))
colnames(expression_matrix_MCF10A_TGFbeta1_8day_batch1_forDay) <- names_TGFbeta1_8day_batch1

expression_matrix_MCF10A_notreated_batch1_forDay$RowNames <- rownames(expression_matrix_MCF10A_notreated_batch1_forDay)
expression_matrix_MCF10A_notreated_batch2_forDay$RowNames <- rownames(expression_matrix_MCF10A_notreated_batch2_forDay)
expression_matrix_MCF10A_TGFbeta1_1day_batch2_forDay$RowNames <- rownames(expression_matrix_MCF10A_TGFbeta1_1day_batch2_forDay)
expression_matrix_MCF10A_TGFbeta1_2day_batch2_forDay$RowNames <- rownames(expression_matrix_MCF10A_TGFbeta1_2day_batch2_forDay)
expression_matrix_MCF10A_TGFbeta1_3day_batch2_forDay$RowNames <- rownames(expression_matrix_MCF10A_TGFbeta1_3day_batch2_forDay)
expression_matrix_MCF10A_TGFbeta1_4day_batch1_forDay$RowNames <- rownames(expression_matrix_MCF10A_TGFbeta1_4day_batch1_forDay)
expression_matrix_MCF10A_TGFbeta1_8day_batch1_forDay$RowNames <- rownames(expression_matrix_MCF10A_TGFbeta1_8day_batch1_forDay)

expression_matrix_MCF10A_forDay <- 
  full_join(expression_matrix_MCF10A_notreated_batch1_forDay, expression_matrix_MCF10A_notreated_batch2_forDay, by = "RowNames") %>%
  full_join(expression_matrix_MCF10A_TGFbeta1_1day_batch2_forDay, by = "RowNames") %>%
  full_join(expression_matrix_MCF10A_TGFbeta1_2day_batch2_forDay, by = "RowNames") %>%
  full_join(expression_matrix_MCF10A_TGFbeta1_3day_batch2_forDay, by = "RowNames") %>%
  full_join(expression_matrix_MCF10A_TGFbeta1_4day_batch1_forDay, by = "RowNames") %>%
  full_join(expression_matrix_MCF10A_TGFbeta1_8day_batch1_forDay, by = "RowNames")


# Preenchendo os NAs com 0
expression_matrix_MCF10A_forDay[is.na(expression_matrix_MCF10A_forDay)] <- 0


# Movendo os nomes das linhas de volta para rownames e removendo a coluna auxiliar
rownames(expression_matrix_MCF10A_forDay) <- expression_matrix_MCF10A_forDay$RowNames
expression_matrix_MCF10A_forDay$RowNames <- NULL

expression_matrix_MCF10A_notreated_batch1_forDay$RowNames <- NULL
expression_matrix_MCF10A_notreated_batch2_forDay$RowNames <- NULL
expression_matrix_MCF10A_TGFbeta1_1day_batch2_forDay$RowNames <- NULL
expression_matrix_MCF10A_TGFbeta1_2day_batch2_forDay$RowNames <- NULL
expression_matrix_MCF10A_TGFbeta1_3day_batch2_forDay$RowNames <- NULL
expression_matrix_MCF10A_TGFbeta1_4day_batch1_forDay$RowNames <- NULL
expression_matrix_MCF10A_TGFbeta1_8day_batch1_forDay$RowNames <- NULL

View(expression_matrix_MCF10A_forDay)
write.table(expression_matrix_MCF10A_forDay, 
            file = "expression_matrix_MCF10A_forDay_ENSG.txt", 
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
save.image("AposMatrizFullForDay.RData")

# Mudando os nomes das linhas da matriz completa de gene para proteína, usando BioMart




mean_t_notreated_batch1 <- rowMeans(t_notreated_batch1@transcriptogramS2[, -c(1,2)])
sd_t_notreated_batch1 <- apply(t_notreated_batch1@transcriptogramS2[, -c(1,2)], 1, sd)
mean_sd_t_notreated_batch1 <- cbind(mean_t_notreated_batch1, sd_t_notreated_batch1)
full_mean_sd_t_notreated_batch1 <- cbind(t_notreated_batch1@transcriptogramS2, mean_sd_t_notreated_batch1)
full_mean_sd_t_notreated_batch1 <- full_mean_sd_t_notreated_batch1[, c(1, 2, ncol(full_mean_sd_t_notreated_batch1)-1, ncol(full_mean_sd_t_notreated_batch1))]
write.table(full_mean_sd_t_notreated_batch1, file = "full_mean_sd_t_notreated_batch1.txt", row.names = F, col.names = T, sep = "\t")

mean_t_notreated_batch2 <- rowMeans(t_notreated_batch2@transcriptogramS2[, -c(1,2)])
sd_t_notreated_batch2 <- apply(t_notreated_batch2@transcriptogramS2[, -c(1,2)], 1, sd)
mean_sd_t_notreated_batch2 <- cbind(mean_t_notreated_batch2, sd_t_notreated_batch2)
full_mean_sd_t_notreated_batch2 <- cbind(t_notreated_batch2@transcriptogramS2, mean_sd_t_notreated_batch2)
full_mean_sd_t_notreated_batch2 <- full_mean_sd_t_notreated_batch2[, c(1, 2, ncol(full_mean_sd_t_notreated_batch2)-1, ncol(full_mean_sd_t_notreated_batch2))]
write.table(full_mean_sd_t_notreated_batch2, file = "full_mean_sd_t_notreated_batch2.txt", row.names = F, col.names = T, sep = "\t")

mean_t_TGFbeta1_1day_batch2 <- rowMeans(t_TGFbeta1_1day_batch2@transcriptogramS2[, -c(1,2)])
sd_t_TGFbeta1_1day_batch2 <- apply(t_TGFbeta1_1day_batch2@transcriptogramS2[, -c(1,2)], 1, sd)
mean_sd_t_TGFbeta1_1day_batch2 <- cbind(mean_t_TGFbeta1_1day_batch2, sd_t_TGFbeta1_1day_batch2)
full_mean_sd_t_TGFbeta1_1day_batch2 <- cbind(t_TGFbeta1_1day_batch2@transcriptogramS2, mean_sd_t_TGFbeta1_1day_batch2)
full_mean_sd_t_TGFbeta1_1day_batch2 <- full_mean_sd_t_TGFbeta1_1day_batch2[, c(1, 2, ncol(full_mean_sd_t_TGFbeta1_1day_batch2)-1, ncol(full_mean_sd_t_TGFbeta1_1day_batch2))]
write.table(full_mean_sd_t_TGFbeta1_1day_batch2, file = "full_mean_sd_t_TGFbeta1_1day_batch2.txt", row.names = F, col.names = T, sep = "\t")

mean_t_TGFbeta1_2day_batch2 <- rowMeans(t_TGFbeta1_2day_batch2@transcriptogramS2[, -c(1,2)])
sd_t_TGFbeta1_2day_batch2 <- apply(t_TGFbeta1_2day_batch2@transcriptogramS2[, -c(1,2)], 1, sd)
mean_sd_t_TGFbeta1_2day_batch2 <- cbind(mean_t_TGFbeta1_2day_batch2, sd_t_TGFbeta1_2day_batch2)
full_mean_sd_t_TGFbeta1_2day_batch2 <- cbind(t_TGFbeta1_2day_batch2@transcriptogramS2, mean_sd_t_TGFbeta1_2day_batch2)
full_mean_sd_t_TGFbeta1_2day_batch2 <- full_mean_sd_t_TGFbeta1_2day_batch2[, c(1, 2, ncol(full_mean_sd_t_TGFbeta1_2day_batch2)-1, ncol(full_mean_sd_t_TGFbeta1_2day_batch2))]
write.table(full_mean_sd_t_TGFbeta1_2day_batch2, file = "full_mean_sd_t_TGFbeta1_2day_batch2.txt", row.names = F, col.names = T, sep = "\t")

mean_t_TGFbeta1_3day_batch2 <- rowMeans(t_TGFbeta1_3day_batch2@transcriptogramS2[, -c(1,2)])
sd_t_TGFbeta1_3day_batch2 <- apply(t_TGFbeta1_3day_batch2@transcriptogramS2[, -c(1,2)], 1, sd)
mean_sd_t_TGFbeta1_3day_batch2 <- cbind(mean_t_TGFbeta1_3day_batch2, sd_t_TGFbeta1_3day_batch2)
full_mean_sd_t_TGFbeta1_3day_batch2 <- cbind(t_TGFbeta1_3day_batch2@transcriptogramS2, mean_sd_t_TGFbeta1_3day_batch2)
full_mean_sd_t_TGFbeta1_3day_batch2 <- full_mean_sd_t_TGFbeta1_3day_batch2[, c(1, 2, ncol(full_mean_sd_t_TGFbeta1_3day_batch2)-1, ncol(full_mean_sd_t_TGFbeta1_3day_batch2))]
write.table(full_mean_sd_t_TGFbeta1_3day_batch2, file = "full_mean_sd_t_TGFbeta1_3day_batch2.txt", row.names = F, col.names = T, sep = "\t")

mean_t_TGFbeta1_4day_batch1 <- rowMeans(t_TGFbeta1_4day_batch1@transcriptogramS2[, -c(1,2)])
sd_t_TGFbeta1_4day_batch1 <- apply(t_TGFbeta1_4day_batch1@transcriptogramS2[, -c(1,2)], 1, sd)
mean_sd_t_TGFbeta1_4day_batch1 <- cbind(mean_t_TGFbeta1_4day_batch1, sd_t_TGFbeta1_4day_batch1)
full_mean_sd_t_TGFbeta1_4day_batch1 <- cbind(t_TGFbeta1_4day_batch1@transcriptogramS2, mean_sd_t_TGFbeta1_4day_batch1)
full_mean_sd_t_TGFbeta1_4day_batch1 <- full_mean_sd_t_TGFbeta1_4day_batch1[, c(1, 2, ncol(full_mean_sd_t_TGFbeta1_4day_batch1)-1, ncol(full_mean_sd_t_TGFbeta1_4day_batch1))]
write.table(full_mean_sd_t_TGFbeta1_4day_batch1, file = "full_mean_sd_t_TGFbeta1_4day_batch1.txt", row.names = F, col.names = T, sep = "\t")

mean_t_TGFbeta1_8day_batch1 <- rowMeans(t_TGFbeta1_8day_batch1@transcriptogramS2[, -c(1,2)])
sd_t_TGFbeta1_8day_batch1 <- apply(t_TGFbeta1_8day_batch1@transcriptogramS2[, -c(1,2)], 1, sd)
mean_sd_t_TGFbeta1_8day_batch1 <- cbind(mean_t_TGFbeta1_8day_batch1, sd_t_TGFbeta1_8day_batch1)
full_mean_sd_t_TGFbeta1_8day_batch1 <- cbind(t_TGFbeta1_8day_batch1@transcriptogramS2, mean_sd_t_TGFbeta1_8day_batch1)
full_mean_sd_t_TGFbeta1_8day_batch1 <- full_mean_sd_t_TGFbeta1_8day_batch1[, c(1, 2, ncol(full_mean_sd_t_TGFbeta1_8day_batch1)-1, ncol(full_mean_sd_t_TGFbeta1_8day_batch1))]
write.table(full_mean_sd_t_TGFbeta1_8day_batch1, file = "full_mean_sd_t_TGFbeta1_8day_batch1.txt", row.names = F, col.names = T, sep = "\t")


# Somente os transcriptogramas.

png(file = "images/t_notreated_batch1_HeatMap.png", width = 1920, height = 1080, res = 150)
Heatmap(t_notreated_batch1@transcriptogramS2)
dev.off()



png(file = "images/t_notreated_batch2_HeatMap.png", width = 1920, height = 1080, res = 150)
Heatmap(t_notreated_batch2@transcriptogramS2, parallel = 40)
dev.off()



png(file = "images/t_TGFbeta1_1day_batch2_HeatMap.png", width = 1920, height = 1080, res = 150)
Heatmap(t_TGFbeta1_1day_batch2@transcriptogramS2, parallel = 40)
dev.off()


png(file = "images/t_TGFbeta1_2day_batch2_HeatMap.png", width = 1920, height = 1080, res = 150)
Heatmap(t_TGFbeta1_2day_batch2@transcriptogramS2, parallel = 40)
dev.off()



png(file = "images/t_TGFbeta1_3day_batch2_HeatMap.png", width = 1920, height = 1080, res = 150)
Heatmap(t_TGFbeta1_3day_batch2@transcriptogramS2, parallel = 40)
dev.off()



png(file = "images/t_TGFbeta1_4day_batch1_HeatMap.png", width = 1920, height = 1080, res = 150)
Heatmap(t_TGFbeta1_4day_batch1@transcriptogramS2, parallel = 40)
dev.off()



png(file = "images/t_TGFbeta1_8day_batch1_HeatMap.png", width = 1920, height = 1080, res = 150)
Heatmap(t_TGFbeta1_8day_batch1@transcriptogramS2, parallel = 40)
dev.off()





# t_notreated_batch1 <- differentiallyExpressed(object = t_notreated_batch1, levels = levels, pValue = 0.01,
#                                               trend = FALSE, title = "radius 80")
# 
# t_notreated_batch2 <- differentiallyExpressed(object = t_notreated_batch2, pValue = 0.01,
#                                               trend = FALSE, title = "radius 80")
# 
# t_TGFbeta1_1day_batch2 <- differentiallyExpressed(object = t_TGFbeta1_1day_batch2, pValue = 0.01,
#                                               trend = FALSE, title = "radius 80")
# 
# t_TGFbeta1_2day_batch2 <- differentiallyExpressed(object = t_TGFbeta1_2day_batch2, pValue = 0.01,
#                                               trend = FALSE, title = "radius 80")
# 
# t_TGFbeta1_3day_batch2 <- differentiallyExpressed(object = t_TGFbeta1_3day_batch2, pValue = 0.01,
#                                               trend = FALSE, title = "radius 80")
# 
# t_TGFbeta1_4day_batch1 <- differentiallyExpressed(object = t_TGFbeta1_4day_batch1, pValue = 0.01,
#                                               trend = FALSE, title = "radius 80")
# 
# t_TGFbeta1_8day_batch1 <- differentiallyExpressed(object = t_TGFbeta1_8day_batch1, pValue = 0.01,
#                                               trend = FALSE, title = "radius 80")

# Transcriptogramer de todos os dias.

# t <- transcriptogramPreprocess(association = assoc, ordering = ord$Protein, radius = 80)
# t <- transcriptogramStep1(object = t, expression = expression_matrix_MCF10A, dictionary = dictionary, nCores = 10)
# t <- transcriptogramStep2(object = t, nCores = 10)

t_mean_sd <- transcriptogramPreprocess(association = assoc, ordering = ord$Protein, radius = 80)
t_mean_sd <- transcriptogramStep1(object = t_mean_sd, expression = mean_sd_expression_matrix_MCF10A, dictionary = dictionary, nCores = 10)
t_mean_sd <- transcriptogramStep2(object = t_mean_sd, nCores = 10)


# Transcriptogramer com todos os dias juntos.
# save.image("AposTranscriptogramasDiarios.RData")
# load("AposTranscriptogramasDiarios.RData")


# subset <- expression_matrix[, c(1:10)]
# write.csv(subset, file= "~/mestrado-single-cell/expressionMatrixDayZero.csv", row.names = T)
# Subamostragem ou agregação pode ser necessária devido ao grande volume de dados
# Aqui está um exemplo simples de subamostragem (selecionando um subconjunto de células)
# set.seed(123)
# subsampled_matrix <- expression_matrix[, sample(ncol(expression_matrix), 500)]

# 3. Análise com transcriptogramer

# Carregar dados de genes.
# genes <- read.table("mtx_conversions/notreated_batch1/outs/filtered_feature_bc_matrix/genes.tsv", header = FALSE, sep = "\t")
# barcodes <- read.table("mtx_conversions/notreated_batch1/outs/filtered_feature_bc_matrix/barcodes.tsv", header = FALSE, sep = "\t")

# Conectando ao Ensembl BioMart
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Consultando informações de genes
# dictionary <- getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id"),
#                     filters = "ensembl_gene_id", values = genes, mart = ensembl)


# dictionary <- dictionary[!apply(dictionary, 1, function(row) any(is.na(row) | row == "")), ]

# colnames(dictionary) <- NULL
# 
# head(dictionary)
# 
# class(dictionary)
# 
# dim(dictionary)
# 
# proteinas <- dictionary[,2]
# 
# head(proteinas)
# 
# write.csv(proteinas, file = "proteinasLista.csv", row.names = F, quote = F)
# 
# head(expression_matrix)
# 
# class(expression_matrix)
# 
# head(expression_matrix[, 'AAACCTGCATTGGTAC-1'], digits = 5, 20)
# 
# AAACCTGCATTGGTAC <- expression_matrix[, 'AAACCTGCATTGGTAC-1', drop = FALSE]
# 
# AAACCTGCATTGGTAC <- as.data.frame(AAACCTGCATTGGTAC)
# 
# head(AAACCTGCATTGGTAC, digits = 5, 20)
# 
# class(AAACCTGCATTGGTAC)
# 
# colnames(dictionary) <- c("V1","V2")
# 
# head(genes)
# 
# head(dictionary)
# 
# merged <- merge(dictionary, genes, by = "V1",all = TRUE)
# 
# head(merged)
# 
# dictionaryP <- merged[,2:3]
# 
# colnames(dictionaryP) <- c("ensembl_peptide_id","hgnc_id")
# 
# head(dictionaryP)
# 
# dim(dictionaryP)
# 
# head(AAACCTGCATTGGTAC, nsmall = 5, 20)
# 
# dim(AAACCTGCATTGGTAC)
# 
# head(rownames(AAACCTGCATTGGTAC))
# 
# sum(rownames(AAACCTGCATTGGTAC) %in% dictionaryP[,2])
# 
# # Criar um data frame de exemplo
# 
# # Exibir o data frame original
# 
# 
# # Adicionar o termo "9606." à coluna 'ensembl_peptide_id'
# dictionaryP$ensembl_peptide_id <- paste("9606.", dictionaryP$ensembl_peptide_id, sep = "")
# 
# head(dictionaryP)
# 
# 
# ### TRANSCRIPTOGRAMER
# 
# ## creating the object and setting the radius as 80
# t <- transcriptogramPreprocess(association = association, ordering = Hs900,
#                                radius = 80)
## after the preprocessing

## modifying the radius of an existing Transcriptogram object
# radius(object = t) <- 50
# 
# ## getting the radius of an existing Transcriptogram object
# r <- radius(object = t)
# 
# oPropertiesR50 <- orderingProperties(object = t, nCores = 5)
# 
# ## slight change of radius
# radius(object = t) <- 80

## this output is partially different comparing to oPropertiesR50
# oPropertiesR80 <- orderingProperties(object = t, nCores = 10)
# 
# # sum(oPropertiesR50$windowModularity)
# 
# sum(oPropertiesR80$windowModularity)
# 
# cProperties <- connectivityProperties(object = t)
# 
# t <- transcriptogramStep1(object = t, expression = AAACCTGCATTGGTAC,
#                           dictionary = dictionaryP, nCores = 10)
# # t2 <- t
# 
# head(t@transcriptogramS1)
# 
# t <- transcriptogramStep2(object = t)

# # radius(object = t2) <- 50
# 
# t2 <- transcriptogramStep2(object = t2, nCores = 5)

## trend = FALSE for microarray data or voom log2-counts-per-million
## the default value for trend is FALSE

# controle <- dim(expression_matrix_MCF10A_notreated_batch1)[2] + dim(expression_matrix_MCF10A_notreated_batch2)[2]
# caso <- dim(expression_matrix_MCF10A_TGFbeta1_1day_batch2)[2] + 
#   dim(expression_matrix_MCF10A_TGFbeta1_2day_batch2)[2] +
#   dim(expression_matrix_MCF10A_TGFbeta1_3day_batch2)[2] +
#   dim(expression_matrix_MCF10A_TGFbeta1_4day_batch1)[2] +
#   dim(expression_matrix_MCF10A_TGFbeta1_8day_batch1)[2]







#### A PARTIR DAQUI, ESTÁ DANDO ERRO, AFIRMNADO QUE NÃO HÁ DIFERENÇA DE EXPRESSÃO.
levels <- c(rep(FALSE, 4),
            rep(TRUE, 10))


t_mean_sd <- differentiallyExpressed(object = t_mean_sd, levels = levels, pValue = 0.5,
                             trend = FALSE, title = "radius 80")

# ## the radius 50 will affect the output significantly
# t2 <- differentiallyExpressed(object = t2, levels = levels, pValue = 0.01,
#                               species = DEsymbols, title = "radius 50")

## using the species argument to translate ENSEMBL Peptide IDs to Symbols
## Internet connection is required for this command
t_mean_sd <- differentiallyExpressed(object = t_mean_sd, levels = levels, pValue = 0.01,
                             species = "Homo sapiens", title = "radius 80")

## translating ENSEMBL Peptide IDs to Symbols using the DEsymbols dataset
t_mean_sd <- differentiallyExpressed(object = t_mean_sd, levels = levels, pValue = 0.01,
                             species = DEsymbols, title = "radius 80")

DE <- DE(object = t_mean_sd)
# DE2 <- DE(object = t2)
nrow(DE)

nrow(DE)/length(unique(DE$ClusterNumber))

# nrow(DE2)

# nrow(DE2)/length(unique(DE2$ClusterNumber))

rdp <- clusterVisualization(object = t_mean_sd)

## using the HsBPTerms dataset to create the Protein2GO data.frame
t_mean_sd <- clusterEnrichment(object = t_mean_sd, species = HsBPTerms,
                       pValue = 0.005, nCores = 5)


## using the species argument to create the Protein2GO data.frame
## Internet connection is required for this command
t_mean_sd <- clusterEnrichment(object = t_mean_sd, species = "Homo sapiens",
                       pValue = 0.005, nCores = 5)

