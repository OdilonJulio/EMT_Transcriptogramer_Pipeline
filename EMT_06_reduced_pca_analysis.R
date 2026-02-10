# pca_analysis_enhanced.R
# Versão aprimorada do Código 1 com todas as funcionalidades do Código 2

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# 1. Carregar objetos de transcriptograma
load("t_matrix_R0_reduced.RData")
load("t_matrix_R30_reduced.RData")

# 2. Definir funções equivalentes às do Código 2

# Função para preparar os dados para PCA (equivalente à prepare_pca_data)
prepare_pca_data <- function(transcriptogram_obj) {
  df <- transcriptogram_obj@transcriptogramS2
  
  # Remover colunas desnecessárias (equivalente a [,-2] no Código 1)
  df <- df[, !colnames(df) %in% c("Protein", "Position")]
  
  # Verificar NAs (novo)
  if (any(is.na(df))) {
    stop("O dataset contém valores NA. Verifique e limpe antes de prosseguir.")
  }
  
  # Transpor a matriz (equivalente ao Código 1)
  t_df <- t(df)
  
  return(t_df)
}

# Função para rodar PCA com verificações (equivalente à run_pca)
run_pca_enhanced <- function(pca_input) {
  # Remover colunas com variância zero (novo)
  variances <- apply(pca_input, 2, var)
  #pca_input_cleaned <- pca_input[, variances > 0]
  
  if (ncol(pca_input_cleaned) == 0) {
    stop("Todas as colunas têm variância zero.")
  }
  
  # Rodar PCA (equivalente ao Código 1)
  pca_result <- prcomp(pca_input_cleaned, center = TRUE, scale. = FALSE)
  
  # Cálculos adicionais (novo)
  eigenvalues <- pca_result$sdev^2
  explained_variance <- eigenvalues / sum(eigenvalues)
  cumulative_variance <- cumsum(explained_variance)
  
  cat("Eigenvalues (autovalores):\n", eigenvalues, "\n")
  cat("Explained Variance:\n", explained_variance, "\n")
  cat("Sum of Explained Variance:\n", sum(explained_variance), "\n")
  cat("Cumulative Variance:\n", cumulative_variance, "\n")
  
  return(list(
    pca_result = pca_result,
    eigenvalues = eigenvalues,
    explained_variance = explained_variance,
    cumulative_variance = cumulative_variance
  ))
}

# Função para extrair condições (equivalente à extract_condition)
extract_condition <- function(colnames) {
  conditions <- c("notreated-batch1", "notreated-batch2", "TGFbeta1-1day-batch2",
                  "TGFbeta1-2day-batch2", "TGFbeta1-3day-batch2",
                  "TGFbeta1-4day-batch1", "TGFbeta1-8day-batch1")
  
  sapply(colnames, function(col) {
    match <- grep(paste(conditions, collapse = "|"), col, value = TRUE)
    if (length(match) > 0) {
      return(sub("_.*", "", match[1]))
    } else {
      return("Unknown")
    }
  })
}

# 3. Processamento dos dados (equivalente ao Código 2)

# Preparar dados para PCA
df_R0 <- prepare_pca_data(t_matrix_R0)
df_R30 <- prepare_pca_data(t_matrix_R30)

# Extrair condições
condition_labels_R0 <- extract_condition(rownames(df_R0))
condition_labels_R30 <- extract_condition(rownames(df_R30))

# Definir paleta de cores (equivalente ao Código 2)
color_palette <- c(
  "notreated-batch1" = "#984EA3",
  "notreated-batch2" = "#377EB8",
  "TGFbeta1-1day-batch2" = "#4DAF4A",
  "TGFbeta1-2day-batch2" = "#3A5400",
  "TGFbeta1-3day-batch2" = "#A65628",
  "TGFbeta1-4day-batch1" = "#FF7F00",
  "TGFbeta1-8day-batch1" = "#E41A1C"
)

# 4. Executar PCA com as novas funções
pca_result_R0 <- run_pca_enhanced(df_R0)
pca_result_R30 <- run_pca_enhanced(df_R30)

# 5. Adicionar metadados (equivalente à add_day_column)
add_day_column <- function(pca_df, sample_names) {
  pca_df$Day <- case_when(
    grepl("notreated-batch1", sample_names) ~ "notreated-batch1",
    grepl("notreated-batch2", sample_names) ~ "notreated-batch2",
    grepl("TGFbeta1-1day-batch2", sample_names) ~ "TGFbeta1-1day-batch2",
    grepl("TGFbeta1-2day-batch2", sample_names) ~ "TGFbeta1-2day-batch2",
    grepl("TGFbeta1-3day-batch2", sample_names) ~ "TGFbeta1-3day-batch2",
    grepl("TGFbeta1-4day-batch1", sample_names) ~ "TGFbeta1-4day-batch1",
    grepl("TGFbeta1-8day-batch1", sample_names) ~ "TGFbeta1-8day-batch1",
    TRUE ~ "Unknown"
  )
  
  if (any(pca_df$Day == "Unknown")) {
    warning(sum(pca_df$Day == "Unknown"), " amostras não mapeadas")
  }
  
  return(pca_df)
}

# Criar data frames com componentes principais
pca_r0_df <- as.data.frame(pca_result_R0$pca_result$x)
pca_r30_df <- as.data.frame(pca_result_R30$pca_result$x)

# Adicionar informações de dia
pca_r0_df <- add_day_column(pca_r0_df, rownames(pca_r0_df))
pca_r30_df <- add_day_column(pca_r30_df, rownames(pca_r30_df))

# 6. Salvar resultados (equivalente ao Código 1 + adicional)
save(pca_result_R0, file = "pca_result_R0_reduced.RData")
save(pca_result_R30, file = "pca_result_R30_reduced.RData")
save(pca_r0_df, pca_r30_df, file = "pca_results_with_metadata_reduced.RData")

# 7. Gerar gráficos (equivalente ao plot_pca.R)
generate_pca_plot <- function(pca_df, title) {
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Day)) +
    geom_point(size = 0.5) +
    scale_color_manual(values = color_palette) +
    ggtitle(title) +
    theme_minimal()
}

plot_R0 <- generate_pca_plot(pca_r0_df, "PCA - Condição R0")
plot_R30 <- generate_pca_plot(pca_r30_df, "PCA - Condição R30")

# Combinar gráficos
combined_plots <- plot_R0 + plot_R30 + plot_layout(guides = 'collect')

# Salvar gráficos
ggsave("pca_combined_plot_reduced.png", combined_plots, width = 12, height = 6)
