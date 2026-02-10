# reconstruct_transcriptogram_final.R
# Reconstrução do transcriptograma com correções de transposição e estrutura

library(transcriptogramer)
library(dplyr)

# 1. Carregar dados necessários
load("pca_result_R0_reduced.RData")  # Objeto pca_result_R0
load("pca_result_R30_reduced.RData") # Objeto pca_result_R30
load("t_matrix_R0_reduced.RData")    # Objeto t_matrix_R0 original
load("t_matrix_R30_reduced.RData")   # Objeto t_matrix_R30 original

# 2. Função para preparar a matriz do transcriptograma original
prepare_original_matrix <- function(transcriptogram_obj) {
  # Extrair dados e definir nomes das linhas
  df <- transcriptogram_obj@transcriptogramS2
  rownames(df) <- df$Protein  # Usar coluna Protein como nomes das linhas
  
  # Remover colunas Protein e Position
  df <- df[, !colnames(df) %in% c("Protein", "Position")]
  
  # Converter para matriz e transpor (genes nas colunas, amostras nas linhas)
  t(as.matrix(df))
}

# 3. Função para reconstruir o transcriptograma
reconstruct_transcriptogram <- function(pca_result, original_matrix) {
  # Extrair componentes e parâmetros da PCA
  rotation <- pca_result$pca_result$rotation
  components <- pca_result$pca_result$x
  center <- pca_result$pca_result$center
  
  # Reconstruir a matriz (transposta para genes x amostras)
  reconstructed <- t(components %*% t(rotation))
  
  # Adicionar o centro de volta (se a PCA foi centralizada)
  if (!is.null(center)) {
    reconstructed <- sweep(reconstructed, 1, center, "+")
  }
  
  # Transpor para amostras x genes (como no original)
  reconstructed <- t(reconstructed)
  
  # Verificar dimensões
  if (!all(dim(reconstructed) == dim(original_matrix))) {
    warning(paste("Dimensões não coincidem: Reconstruído", 
                  paste(dim(reconstructed), collapse = "x"),
                  "vs Original", paste(dim(original_matrix), collapse = "x")))
  }
  
  return(reconstructed)
}

# 4. Preparar matrizes originais
original_matrix_R0 <- prepare_original_matrix(t_matrix_R0)
original_matrix_R30 <- prepare_original_matrix(t_matrix_R30)

# 5. Reconstruir transcriptogramas
reconstructed_matrix_R0 <- reconstruct_transcriptogram(pca_result_R0, original_matrix_R0)
reconstructed_matrix_R30 <- reconstruct_transcriptogram(pca_result_R30, original_matrix_R30)

# 6. Criar novos objetos transcriptograma com a estrutura correta
create_transcriptogram_object <- function(reconstructed_matrix, original_obj) {
  # Converter para data frame
  df <- as.data.frame(t(reconstructed_matrix))  # Transpor para genes x amostras
  
  # Adicionar colunas Protein e Position do objeto original
  df$Protein <- original_obj@transcriptogramS2$Protein
  df$Position <- original_obj@transcriptogramS2$Position
  
  # Reordenar colunas
  df <- df[, c("Protein", "Position", setdiff(colnames(df), c("Protein", "Position")))]
  
  # Criar novo objeto mantendo outros slots do original
  new_obj <- original_obj
  new_obj@transcriptogramS2 <- df
  
  return(new_obj)
}

reconstructed_R0 <- create_transcriptogram_object(reconstructed_matrix_R0, t_matrix_R0)
reconstructed_R30 <- create_transcriptogram_object(reconstructed_matrix_R30, t_matrix_R30)

# 7. Verificação dimensional
cat("=== Verificação R0 ===\n")
cat("Original:", dim(original_matrix_R0), "\n")
cat("Reconstruído:", dim(reconstructed_matrix_R0), "\n\n")

cat("=== Verificação R30 ===\n")
cat("Original:", dim(original_matrix_R30), "\n")
cat("Reconstruído:", dim(reconstructed_matrix_R30), "\n")


# 9. Salvar objetos reconstruídos
save(reconstructed_R0, file = "reconstructed_transcriptogram_R0_final.RData")
save(reconstructed_R30, file = "reconstructed_transcriptogram_R30_final.RData")

# 10. Salvar matrizes reconstruídas (opcional)
save(reconstructed_matrix_R0, file = "reconstructed_matrix_R0.RData")
save(reconstructed_matrix_R30, file = "reconstructed_matrix_R30.RData")

cat("\nProcesso concluído. Objetos salvos com sucesso.\n")
View(reconstructed_matrix_R0)
View(reconstructed_matrix_R30)
View(original_matrix_R0)
View(original_matrix_30)
