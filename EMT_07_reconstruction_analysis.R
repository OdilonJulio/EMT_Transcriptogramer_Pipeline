# reconstruct_transcriptogram_final.R
# Reconstrução do transcriptograma com correções de transposição e estrutura

library(transcriptogramer)
library(dplyr)

# 1. Carregar dados necessários
load("pca_result_R0.RData")  # Objeto pca_result_R0
load("pca_result_R30.RData") # Objeto pca_result_R30
load("t_matrix_R0.RData")    # Objeto t_matrix_R0 original
load("t_matrix_R30.RData")   # Objeto t_matrix_R30 original

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


#### GRÁFICOS

# === 11. Visualizações de validação (PDF alta qualidade) ===
library(ggplot2)
library(pheatmap)

# --- 11.1 Scatter plot Original vs Reconstruído (R0) ---
df_compare_R0 <- data.frame(
  Original = as.vector(original_matrix_R0),
  Reconstruido = as.vector(reconstructed_matrix_R0)
)

p_R0 <- ggplot(df_compare_R0, aes(x = Original, y = Reconstruido)) +
  geom_point(alpha = 0.3, size = 0.6, color = "#1f78b4") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  theme_classic(base_size = 14) +
  labs(title = "R0: Original vs Reconstruído",
       x = "Expressão Original", y = "Expressão Reconstruída")

ggsave("images/scatter_reconstruction_analysis_R0.pdf", p_R0, width = 6, height = 5, dpi = 300)

# --- 11.2 Scatter plot Original vs Reconstruído (R30) ---
df_compare_R30 <- data.frame(
  Original = as.vector(original_matrix_R30),
  Reconstruido = as.vector(reconstructed_matrix_R30)
)

p_R30 <- ggplot(df_compare_R30, aes(x = Original, y = Reconstruido)) +
  geom_point(alpha = 0.3, size = 0.6, color = "#33a02c") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  theme_classic(base_size = 14) +
  labs(title = "R30: Original vs Reconstruído",
       x = "Expressão Original", y = "Expressão Reconstruída")

ggsave("images/scatter_reconstruction_analysis_R30.pdf", p_R30, width = 6, height = 5, dpi = 300)

# --- 11.3 Heatmap das diferenças (R0) ---
diff_matrix_R0 <- original_matrix_R0 - reconstructed_matrix_R0
pdf("images/heatmap_diff_reconstruction_analysis_R0.pdf", width = 7, height = 6)
pheatmap(diff_matrix_R0,
         color = colorRampPalette(c("blue","white","red"))(50),
         show_rownames = FALSE, show_colnames = FALSE,
         main = "Diferença R0 (Original - Reconstruído)")
dev.off()

# --- 11.4 Heatmap das diferenças (R30) ---
diff_matrix_R30 <- original_matrix_R30 - reconstructed_matrix_R30
pdf("images/heatmap_diff_reconstruction_analysis_R30.pdf", width = 7, height = 6)
pheatmap(diff_matrix_R30,
         color = colorRampPalette(c("blue","white","red"))(50),
         show_rownames = FALSE, show_colnames = FALSE,
         main = "Diferença R30 (Original - Reconstruído)")
dev.off()

# --- 11.5 Curva de variância explicada (R0) ---
var_expl_R0 <- (pca_result_R0$pca_result$sdev^2) / sum(pca_result_R0$pca_result$sdev^2)
df_var_R0 <- data.frame(PC = seq_along(var_expl_R0),
                        Variancia = cumsum(var_expl_R0))

p_var_R0 <- ggplot(df_var_R0, aes(x = PC, y = Variancia)) +
  geom_line(color = "#1f78b4", size = 1) +
  geom_point(color = "#1f78b4") +
  theme_classic(base_size = 14) +
  labs(title = "Curva de Variância Explicada (R0)",
       x = "Componentes Principais", y = "Variância Acumulada")

ggsave("images/variancia_reconstruction_analysis_R0.pdf", p_var_R0, width = 6, height = 4.5, dpi = 300)

# --- 11.6 R² global ---
R2_R0 <- cor(as.vector(original_matrix_R0), as.vector(reconstructed_matrix_R0))^2
R2_R30 <- cor(as.vector(original_matrix_R30), as.vector(reconstructed_matrix_R30))^2
cat("R² R0:", R2_R0, "\n")
cat("R² R30:", R2_R30, "\n")
