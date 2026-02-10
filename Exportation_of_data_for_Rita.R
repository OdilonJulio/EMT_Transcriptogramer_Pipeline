## ==========================================
## Script para exportação de dados tabulados
## ==========================================

# Carregar bibliotecas necessárias
library(data.table)

## ==========================================
## 1. Função para exportar matrizes
## ==========================================
export_matrix <- function(mat, filename, rownames_col = "Gene") {
  # Converter para data.table mantendo os nomes das linhas
  dt <- as.data.table(mat, keep.rownames = rownames_col)
  
  # Definir nome do arquivo de saída
  output_file <- paste0(filename, ".tsv")
  
  # Exportar com tabulação
  fwrite(dt, file = output_file, sep = "\t", quote = FALSE, na = "NA")
  
  message(paste("Matriz exportada para:", output_file))
}

## ==========================================
## 2. Função para exportar resultados de análise
## ==========================================
export_analysis_results <- function(analysis_obj, prefix) {
  # Exportar dados de R² por número de PCs
  r2_data <- analysis_obj$results
  fwrite(r2_data, file = paste0(prefix, "_r2_by_pcs.tsv"), sep = "\t")
  
  message(paste("Resultados de análise exportados para:", prefix, "*"))
}

## ==========================================
## 3. Função para exportar dados de reconstrução completa
## ==========================================
export_full_reconstruction <- function(recon_obj, prefix) {
  # Criar lista de resultados para exportação
  results <- list(
    global_r2 = data.frame(Global_R2 = recon_obj$global),
    by_gene = data.frame(Gene = names(recon_obj$by_gene), R2 = recon_obj$by_gene),
    by_cell = data.frame(Cell = names(recon_obj$by_cell), R2 = recon_obj$by_cell),
    zero_adjustment = data.frame(
      Adjusted_Zero_Genes = recon_obj$adjusted_zero_genes,
      Adjustment_Value = recon_obj$adjustment_value
    )
  )
  
  # Exportar cada componente
  fwrite(results$global_r2, file = paste0(prefix, "_global_r2.tsv"), sep = "\t")
  fwrite(results$by_gene, file = paste0(prefix, "_r2_by_gene.tsv"), sep = "\t")
  fwrite(results$by_cell, file = paste0(prefix, "_r2_by_cell.tsv"), sep = "\t")
  fwrite(results$zero_adjustment, file = paste0(prefix, "_zero_adjustment.tsv"), sep = "\t")
  
  message(paste("Dados completos de reconstrução exportados para:", prefix, "*"))
}

prepare_original_matrix <- function(obj) {
  mat <- as.matrix(obj@transcriptogramS2[, !names(obj@transcriptogramS2) %in% c("Protein", "Position")])
  rownames(mat) <- obj@transcriptogramS2$Protein
  gc()
  t(mat)
}

## ==========================================
## 4. Execução principal
## ==========================================

# Carregar os dados
load_data()
load("~/mestrado/analysis_R0.RData")
load("~/mestrado/analysis_R30.RData")

# Preparar matrizes originais

original_matrix_R0 <- prepare_original_matrix(t_matrix_R0)
original_matrix_R30 <- prepare_original_matrix(t_matrix_R30)

# Exportar matrizes originais
export_matrix(original_matrix_R0, "original_matrix_R0")
export_matrix(original_matrix_R30, "original_matrix_R30")

# Exportar matrizes reconstruídas
export_matrix(reconstructed_matrix_R0, "reconstructed_matrix_R0")
export_matrix(reconstructed_matrix_R30, "reconstructed_matrix_R30")

# Exportar resultados das análises
export_analysis_results(analysis_R0, "analysis_R0")
export_analysis_results(analysis_R30, "analysis_R30")

# # Exportar dados de reconstrução completa
# export_full_reconstruction(full_recon_r2_R0, "full_recon_R0")
# export_full_reconstruction(full_recon_r2_R30, "full_recon_R30")

# # Exportar componentes principais (opcional)
# export_matrix(pca_result_R0$pca_result$x, "pca_scores_R0", "Cell")
# export_matrix(pca_result_R30$pca_result$x, "pca_scores_R30", "Cell")
# export_matrix(pca_result_R0$pca_result$rotation, "pca_loadings_R0", "Gene")
# export_matrix(pca_result_R30$pca_result$rotation, "pca_loadings_R30", "Gene")

# Exportar variância explicada
export_pca_variance <- function(pca_res, prefix) {
  variance_data <- data.frame(
    PC = 1:length(pca_res$pca_result$sdev),
    Variance = pca_res$pca_result$sdev^2,
    Cumulative = cumsum(pca_res$pca_result$sdev^2)/sum(pca_res$pca_result$sdev^2)
  )
  fwrite(variance_data, file = paste0(prefix, "_pca_variance.tsv"), sep = "\t")
}

export_pca_variance(pca_result_R0, "R0")
export_pca_variance(pca_result_R30, "R30")

message("Exportação de dados concluída com sucesso!")