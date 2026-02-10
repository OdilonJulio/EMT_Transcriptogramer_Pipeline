# load_data.R
# Carrega os dados e cria os objetos Seurat.

library(Seurat)
library(SeuratObject)
library(scDblFinder)
library(SingleCellExperiment)

# Caminhos para as amostras
paths <- c(
  "mtx_conversions/notreated_batch1/outs/filtered_feature_bc_matrix",
  "mtx_conversions/notreated_batch2/outs/filtered_feature_bc_matrix/",
  "mtx_conversions/TGFbeta1_1day_batch2/outs/filtered_feature_bc_matrix/",
  "mtx_conversions/TGFbeta1_2day_batch2/outs/filtered_feature_bc_matrix/",
  "mtx_conversions/TGFbeta1_3day_batch2/outs/filtered_feature_bc_matrix/",
  "mtx_conversions/TGFbeta1_4day_batch1/outs/filtered_feature_bc_matrix/",
  "mtx_conversions/TGFbeta1_8day_batch1/outs/filtered_feature_bc_matrix/"
)

# Nomes para os objetos
sample_names <- c(
  "notreated-batch1", "notreated-batch2", "TGFbeta1-1day-batch2",
  "TGFbeta1-2day-batch2", "TGFbeta1-3day-batch2", "TGFbeta1-4day-batch1", "TGFbeta1-8day-batch1"
)

# Criar objetos Seurat e adicionar anotações
seurat_objects <- lapply(1:length(paths), function(i) {
  data <- Read10X(paths[i])
  seurat <- CreateSeuratObject(counts = data, project = sample_names[i])
  seurat$batch <- ifelse(i %in% c(1, 6, 7), "Batch1", "Batch2")
  seurat$day <- c("0", "0", "1", "2", "3", "4", "8")[i]
  return(seurat)
})

# Salvar objetos Seurat
save(seurat_objects, file = "seurat_objects.RData")
