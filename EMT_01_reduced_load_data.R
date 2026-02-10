# load_data_subset.R
# Carrega um subconjunto de 70 células (10 por amostra) para testes rápidos

library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)

# 1. Definir caminhos e metadados
paths <- c(
  "mtx_conversions/notreated_batch1/outs/filtered_feature_bc_matrix",
  "mtx_conversions/notreated_batch2/outs/filtered_feature_bc_matrix/",
  "mtx_conversions/TGFbeta1_1day_batch2/outs/filtered_feature_bc_matrix/",
  "mtx_conversions/TGFbeta1_2day_batch2/outs/filtered_feature_bc_matrix/",
  "mtx_conversions/TGFbeta1_3day_batch2/outs/filtered_feature_bc_matrix/",
  "mtx_conversions/TGFbeta1_4day_batch1/outs/filtered_feature_bc_matrix/",
  "mtx_conversions/TGFbeta1_8day_batch1/outs/filtered_feature_bc_matrix/"
)

sample_names <- c(
  "notreated-batch1", "notreated-batch2", "TGFbeta1-1day-batch2",
  "TGFbeta1-2day-batch2", "TGFbeta1-3day-batch2", "TGFbeta1-4day-batch1", "TGFbeta1-8day-batch1"
)

# 2. Função para carregar subconjunto
load_subset <- function(path, n_cells = 20) {
  data <- Read10X(path)
  set.seed(123) # Reprodutibilidade
  if(ncol(data) < n_cells) {
    message(paste("Amostra com apenas", ncol(data), "células. Usando todas."))
    return(data)
  }
  return(data[, sample(colnames(data), n_cells)])
}

# 3. Criar objetos Seurat com 10 células cada
seurat_list <- lapply(1:7, function(i) {
  subset <- load_subset(paths[i])
  obj <- CreateSeuratObject(
    counts = subset,
    project = sample_names[i],
    min.cells = 3,
    min.features = 200
  )
  obj$batch <- ifelse(i %in% c(1, 6, 7), "Batch1", "Batch2")
  obj$day <- c("0", "0", "1", "2", "3", "4", "8")[i]
  return(obj)
})

# 4. Combinar e salvar
combined <- merge(seurat_list[[1]], y = seurat_list[-1], 
                  add.cell.ids = sample_names)
save(combined, file = "seurat_140cells.RData")

# 5. Verificação
cat("\nObjeto criado com:\n",
    ncol(combined), "células (20 por amostra)\n",
    nrow(combined), "genes\n",
    "Metadados: batch e day\n",
    "Salvo como 'seurat_140cells.RData'\n")
