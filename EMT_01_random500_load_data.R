# load_data_subset_v5_fixed.R
# Versão final corrigida para Seurat v5

library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)

## 1. Configuração das amostras ----
sample_data <- data.frame(
  path = c(
    "mtx_conversions/notreated_batch1/outs/filtered_feature_bc_matrix",
    "mtx_conversions/notreated_batch2/outs/filtered_feature_bc_matrix",
    "mtx_conversions/TGFbeta1_1day_batch2/outs/filtered_feature_bc_matrix",
    "mtx_conversions/TGFbeta1_2day_batch2/outs/filtered_feature_bc_matrix",
    "mtx_conversions/TGFbeta1_3day_batch2/outs/filtered_feature_bc_matrix",
    "mtx_conversions/TGFbeta1_4day_batch1/outs/filtered_feature_bc_matrix",
    "mtx_conversions/TGFbeta1_8day_batch1/outs/filtered_feature_bc_matrix"
  ),
  sample_name = c(
    "notreated-batch1", "notreated-batch2", "TGFbeta1-1day-batch2",
    "TGFbeta1-2day-batch2", "TGFbeta1-3day-batch2", "TGFbeta1-4day-batch1",
    "TGFbeta1-8day-batch1"
  ),
  batch = c("Batch1", "Batch2", "Batch2", "Batch2", "Batch2", "Batch1", "Batch1"),
  day = c("0", "0", "1", "2", "3", "4", "8"),
  stringsAsFactors = FALSE
)

## 2. Função para carregar subconjuntos ----
load_sample <- function(path, sample_name, n_cells = 500) {
  # Carrega os dados
  data <- Read10X(data.dir = path)
  
  # Verifica células disponíveis
  if (ncol(data) < n_cells) {
    message(sprintf("%s: usando todas as %d células", sample_name, ncol(data)))
    cell_ids <- colnames(data)
  } else {
    set.seed(123) # Reprodutibilidade
    cell_ids <- sample(colnames(data), n_cells)
  }
  
  # Cria objeto Seurat
  CreateSeuratObject(
    counts = data[, cell_ids, drop = FALSE],
    project = sample_name,
    min.cells = 3,
    min.features = 200
  )
}

## 3. Pipeline principal ----
# Carrega e processa todas as amostras
seurat_list <- lapply(1:nrow(sample_data), function(i) {
  obj <- load_sample(
    sample_data$path[i],
    sample_data$sample_name[i]
  )
  
  # Adiciona metadados de forma segura
  obj$batch <- sample_data$batch[i]
  obj$day <- sample_data$day[i]
  
  return(obj)
})

# Combina as amostras
combined <- merge(
  x = seurat_list[[1]],
  y = seurat_list[2:length(seurat_list)],
  add.cell.ids = sample_data$sample_name
)

## 4. Detecção de dupletos (método robusto para v5) ----
# Extrai todas as camadas de contagem
all_counts <- lapply(Layers(combined), function(layer) {
  LayerData(combined, layer = layer)
})

# Combina todas as camadas
combined_counts <- do.call(cbind, all_counts)

# Cria objeto SCE garantindo correspondência com metadados
sce <- SingleCellExperiment(
  assays = list(counts = combined_counts),
  colData = combined[[]][colnames(combined_counts), , drop = FALSE]
)

# Detecta dupletos
sce <- scDblFinder(sce)

# Adiciona resultados ao objeto Seurat
combined$scDblFinder_score <- sce$scDblFinder.score[match(colnames(combined), colnames(sce))]
combined$scDblFinder_class <- sce$scDblFinder.class[match(colnames(combined), colnames(sce))]

## 5. Salvar resultados ----
output_file <- "seurat_3500cells_v5_fixed.RData"
save(combined, file = output_file)

## 6. Relatório final ----
cat("=== RESUMO FINAL ===\n")
cat("Total de células:", ncol(combined), "\n")
cat("Total de genes:", nrow(combined), "\n")

cat("\nDistribuição por amostra:\n")
print(table(combined$orig.ident))

cat("\nResumo de dupletos:\n")
print(table(combined$scDblFinder_class))

cat("\nArquivo salvo como:", output_file, "\n")