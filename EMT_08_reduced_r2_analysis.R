# transcriptogram_reconstruction_validation.R
# Cálculo do R² para avaliar qualidade da reconstrução

library(transcriptogramer)
library(dplyr)
library(ggplot2)

# 1. Carregar dados necessários
load("reconstructed_matrix_R0.RData")  # Matriz reconstruída R0
load("reconstructed_matrix_R30.RData") # Matriz reconstruída R30
load("t_matrix_R0_reduced.RData")      # Objeto original R0
load("t_matrix_R30_reduced.RData")     # Objeto original R30
load("pca_result_R0_reduced.RData")    # Resultados PCA R0
load("pca_result_R30_reduced.RData")   # Resultados PCA R30

# 2. Função para preparar matriz original (igual à anterior)
prepare_original_matrix <- function(transcriptogram_obj) {
  df <- transcriptogram_obj@transcriptogramS2
  rownames(df) <- df$Protein
  df <- df[, !colnames(df) %in% c("Protein", "Position")]
  t(as.matrix(df))
}

# 3. Preparar matrizes originais
original_matrix_R0 <- prepare_original_matrix(t_matrix_R0)
original_matrix_R30 <- prepare_original_matrix(t_matrix_R30)

# 4. Função para calcular R² conforme fórmula especificada

calculate_r_squared <- function(original, reconstructed) {
  # Calcular médias dos genes (Tmed_i)
  gene_means <- rowMeans(original)
  
  # Calcular valor típico (moda) para proteção contra zero
  mode_value <- as.numeric(names(sort(table(original), decreasing = TRUE)[1])) / 1000
  
  # Substituir médias zero pela moda/1000
  zero_means <- gene_means == 0
  if (any(zero_means)) {
    gene_means[zero_means] <- mode_value
  }
  
  # Calcular diferenças quadráticas normalizadas
  squared_diff <- (original - reconstructed)^2
  normalized_diff <- squared_diff / (gene_means^2)
  
  # Retornar resultados
  list(
    global = mean(normalized_diff),
    by_gene = rowMeans(normalized_diff),
    by_cell = colMeans(normalized_diff),
    adjusted_zero_genes = sum(zero_means),
    adjustment_value = mode_value
  )
}


# 5. Função para reconstruir com n PCs específico
reconstruct_with_n_pcs <- function(pca_result, original_matrix, n_pcs) {
  rotation <- pca_result$pca_result$rotation[, 1:n_pcs, drop = FALSE]
  components <- pca_result$pca_result$x[, 1:n_pcs, drop = FALSE]
  center <- pca_result$pca_result$center
  
  reconstructed <- t(components %*% t(rotation))
  
  if (!is.null(center)) {
    reconstructed <- sweep(reconstructed, 1, center, "+")
  }
  
  t(reconstructed)  # Retornar no formato amostras x genes
}

# 6. Análise para diferentes números de PCs
analyze_pc_components <- function(pca_result, original_matrix, radius_label) {
  max_pcs <- ncol(pca_result$pca_result$x)
  results <- data.frame(
    n_pcs = 1:max_pcs,
    r_squared = numeric(max_pcs)
  )
  
  for (n in 1:max_pcs) {
    reconstructed <- reconstruct_with_n_pcs(pca_result, original_matrix, n)
    r2 <- calculate_r_squared(original_matrix, reconstructed)
    results$r_squared[n] <- r2$global
  }
  
  # Plotar curva de erro
  plot <- ggplot(results, aes(x = n_pcs, y = r_squared)) +
    geom_line(color = "steelblue") +
    geom_point(color = "steelblue") +
    labs(
      title = paste("Erro de reconstrução por número de PCs -", radius_label),
      x = "Número de Componentes Principais",
      y = "Erro R²",
      caption = paste("Análise para", radius_label)
    ) +
    theme_minimal()
  
  list(
    results = results,
    plot = plot
  )
}

# 7. Executar análises para R0 e R30
analysis_R0 <- analyze_pc_components(pca_result_R0, original_matrix_R0, "R0")
analysis_R30 <- analyze_pc_components(pca_result_R30, original_matrix_R30, "R30")

# 8. Resultados com todos os PCs (reconstrução completa)
full_recon_r2_R0 <- calculate_r_squared(original_matrix_R0, reconstructed_matrix_R0)
full_recon_r2_R30 <- calculate_r_squared(original_matrix_R30, reconstructed_matrix_R30)

# 9. Exibir resultados
cat("=== Resultados para R0 ===\n")
cat("Erro R² global (todos PCs):", full_recon_r2_R0$global, "\n")
print(head(analysis_R0$results))

cat("\n=== Resultados para R30 ===\n")
cat("Erro R² global (todos PCs):", full_recon_r2_R30$global, "\n")
print(head(analysis_R30$results))

# 10. Mostrar gráficos
print(analysis_R0$plot)
print(analysis_R30$plot)

# 11. Salvar resultados
save(analysis_R0, analysis_R30, full_recon_r2_R0, full_recon_r2_R30,
     file = "transcriptogram_reconstruction_validation.RData")

# 12. Salvar gráficos
ggsave("reconstruction_error_R0.png", analysis_R0$plot, width = 8, height = 6)
ggsave("reconstruction_error_R30.png", analysis_R30$plot, width = 8, height = 6)

cat("\nAnálise concluída. Resultados salvos em:\n",
    "- transcriptogram_reconstruction_validation.RData\n",
    "- reconstruction_error_R0.png\n",
    "- reconstruction_error_R30.png\n")


generate_elbow_plot_with_cumulative <- function(pca_result, radius_label, top_n = NULL) {
  # Variância explicada
  explained_var <- pca_result$pca_result$sdev^2
  explained_var_ratio <- explained_var / sum(explained_var)
  cumulative_var <- cumsum(explained_var_ratio)
  
  df <- data.frame(
    PC = 1:length(explained_var_ratio),
    VarianceExplained = explained_var_ratio,
    Cumulative = cumulative_var
  )
  
  if (!is.null(top_n)) {
    df <- df[1:top_n, ]
  }
  
  plot <- ggplot(df, aes(x = PC)) +
    geom_col(aes(y = VarianceExplained), fill = "steelblue", alpha = 0.7) +
    geom_line(aes(y = Cumulative), color = "darkorange", size = 1) +
    geom_point(aes(y = Cumulative), color = "darkorange", size = 2) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1),
      sec.axis = dup_axis(name = "Variância Acumulada")
    ) +
    labs(
      title = paste("Curva de Cotovelo com Acumulado -", radius_label),
      x = "Componente Principal",
      y = "Variância Explicada"
    ) +
    theme_minimal()
  
  return(plot)
}


# Top 10 PCs
elbow_top10_R0 <- generate_elbow_plot_with_cumulative(pca_result_R0, "R0", top_n = 10)
elbow_top10_R30 <- generate_elbow_plot_with_cumulative(pca_result_R30, "R30", top_n = 10)

# Todos os PCs
elbow_all_R0 <- generate_elbow_plot_with_cumulative(pca_result_R0, "R0")
elbow_all_R30 <- generate_elbow_plot_with_cumulative(pca_result_R30, "R30")



ggsave("elbow_plot_with_cumulative_top10_R0.png", elbow_top10_R0, width = 7, height = 5)
ggsave("elbow_plot_with_cumulative_top10_R30.png", elbow_top10_R30, width = 7, height = 5)
ggsave("elbow_plot_with_cumulative_all_R0.png", elbow_all_R0, width = 8, height = 6)
ggsave("elbow_plot_with_cumulative_all_R30.png", elbow_all_R30, width = 8, height = 6)


print(elbow_top10_R0)
print(elbow_top10_R30)
print(elbow_all_R0)
print(elbow_all_R30)

generate_elbow_plot <- function(pca_result, r2_results, radius_label, n_fit = 2) {
  # Preparar dados com escalas normalizadas
  max_r2 <- max(r2_results$r_squared)
  df <- data.frame(
    PC = r2_results$n_pcs,
    R2_Error = r2_results$r_squared / max_r2,  # Normalizado para [0,1]
    CumulativeVar = cumsum(pca_result$pca_result$sdev[1:nrow(r2_results)]^2) / 
      sum(pca_result$pca_result$sdev^2)
  )
  
  # Ajustar retas para o método do cotovelo
  first_fit <- lm(CumulativeVar ~ PC, data = df[1:n_fit, ])
  last_fit <- lm(CumulativeVar ~ PC, data = tail(df, n_fit))
  
  # Calcular ponto de interseção (cotovelo)
  elbow_x <- (coef(last_fit)[[1]] - coef(first_fit)[[1]]) / 
    (coef(first_fit)[[2]] - coef(last_fit)[[2]])
  elbow_y <- coef(first_fit)[[2]] * elbow_x + coef(first_fit)[[1]]
  
  # Criar o gráfico com eixos sincronizados
  p <- ggplot(df, aes(x = factor(PC))) +  # Fator para remover números das barras
    
    # Barras de erro R² (sem números no eixo x)
    geom_col(aes(y = R2_Error * df$CumulativeVar[1]), fill = "#1f77b4", alpha = 0.6, width = 0.7) +
    
    # Linha de variância acumulada (mesma escala)
    geom_line(aes(y = CumulativeVar, group = 1), color = "#ff7f0e", size = 1.2) +
    geom_point(aes(y = CumulativeVar), color = "#ff7f0e", size = 3) +
    
    # Retas do método do cotovelo
    geom_segment(
      aes(
        x = 1, 
        xend = elbow_x, 
        y = predict(first_fit, newdata = data.frame(PC = 1)),
        yend = elbow_y
      ),
      color = "red", linetype = "dashed", size = 0.8
    ) +
    geom_segment(
      aes(
        x = elbow_x, 
        xend = nrow(df),
        y = elbow_y, 
        yend = predict(last_fit, newdata = data.frame(PC = nrow(df)))
      ),
      color = "red", linetype = "dashed", size = 0.8
    ) +
    
    # Ponto de cotovelo
    geom_point(aes(x = elbow_x, y = elbow_y), color = "#2ca02c", size = 4) +
    annotate(
      "text", 
      x = elbow_x, 
      y = elbow_y + 0.05, 
      label = sprintf("PC %.1f", elbow_x),
      color = "#2ca02c", 
      vjust = 0
    ) +
    
    # Formatação
    labs(
      title = paste("Análise de Componentes -", radius_label),
      x = "Número de Componentes Principais",
      y = "Erro R² (normalizado) / Variância Acumulada"
    ) +
    scale_y_continuous(limits = c(0, 1.05)) +
    scale_x_discrete() +  # Remove números das barras
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_blank(),  # Remove rótulos do eixo x
      axis.ticks.x = element_blank()  # Remove ticks do eixo x
    )
  
  return(list(plot = p, elbow_point = c(elbow_x, elbow_y)))
}


# Exemplo de uso:
result_R0 <- generate_elbow_plot(pca_result_R0, analysis_R0$results, "R0")
result_R30 <- generate_elbow_plot(pca_result_R30, analysis_R30$results, "R30")
print(result_R0$plot)
print(result_R30$plot)


# Salvar
ggsave("images/r2_vs_pcs_with_elbow_R0.png", r2_elbow_R0$plot, width = 8, height = 6)
ggsave("images/r2_vs_pcs_with_elbow_R30.png", r2_elbow_R30$plot, width = 8, height = 6)





