# Carregar pacotes
library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)
library(magrittr)
library(transcriptogramer)
library(biomaRt)
library(dplyr)
library(vroom)
library(RcppCWB)
library(Ropj)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(RColorBrewer)

# Carregar os dataframes
load("t_matrix_batch_effect_correction_R0.RData")
load("t_matrix_batch_effect_correction_R30.RData")

##############
#### Funções de PCA
##############

# Função para preparar os dados para PCA
prepare_pca_data <- function(df) {
  if (any(is.na(df))) {
    stop("O dataset contém valores NA. Verifique e limpe antes de prosseguir.")
  }
  df <- df[, !is.na(colnames(df))]  # Remover colunas com rótulos NA
  rownames(df) <- df$Protein        # Usar a coluna "Protein" como nomes das linhas
  df <- df[, -1]                    # Remover a coluna "Protein"
  return(t(df))                     # Transpor a matriz para PCA
}

run_pca <- function(pca_input) {
  # Remover colunas com variância zero
  variances <- apply(pca_input, 2, var)
  pca_input_cleaned <- pca_input[, variances > 0]
  
  if (ncol(pca_input_cleaned) == 0) {
    stop("Todas as colunas têm variância zero.")
  }
  head(pca_result_R0$pca_result$x)
colnames(pca_result_R0$pca_result$x)

  # Rodar PCA
  pca_result <- prcomp(pca_input_cleaned, center = TRUE, scale. = FALSE)
  
  # Calcular autovalores
  eigenvalues <- pca_result$sdev^2
  cat("Eigenvalues (autovalores):\n", eigenvalues, "\n")
  
  # Calcular variância explicada
  explained_variance <- eigenvalues / sum(eigenvalues)
  cat("Explained Variance:\n", explained_variance, "\n")
  cat("Sum of Explained Variance:\n", sum(explained_variance), "\n")
  
  # Calcular variância cumulativa
  cumulative_variance <- cumsum(explained_variance)
  cat("Cumulative Variance:\n", cumulative_variance, "\n")
  
  return(list(
    pca_result = pca_result,
    eigenvalues = eigenvalues,
    explained_variance = explained_variance,
    cumulative_variance = cumulative_variance
  ))
}

# Função para gerar gráficos pairwise
plot_pairwise <- function(pca_matrix, condition_labels, pcs_to_plot = 1:10, color_palette, file_name) {
  pairwise_plots <- list()
  pca_df <- as.data.frame(pca_matrix)
  pca_df$Condition <- factor(condition_labels[rownames(pca_df)])
  
  for (i in pcs_to_plot) {
    pairwise_plots[[i]] <- ggplot(pca_df, aes_string(x = "PC1", y = paste0("PC", i), color = "Condition")) +
      geom_point(size = 1.5, alpha = 0.8) +
      labs(title = paste0("PC1 vs PC", i), x = "PC1", y = paste0("PC", i)) +
      theme_minimal() +
      scale_color_manual(values = color_palette) +
      theme(legend.position = "right")
  }
  
  combined_plot <- patchwork::wrap_plots(pairwise_plots[pcs_to_plot], ncol = 3)
  ggsave(filename = file_name, plot = combined_plot, dpi = 300, width = 15, height = 10)
  
  if (file.exists(file_name)) {
    message("Gráfico salvo com sucesso em ", file_name)
  } else {
    stop("Falha ao salvar o gráfico.")
  }
}

##############
#### Processo de PCA
##############

# Lista de condições para filtrar
conditions <- c("notreated-batch1", "notreated-batch2", "TGFbeta1-1day-batch2", 
                "TGFbeta1-2day-batch2", "TGFbeta1-3day-batch2", 
                "TGFbeta1-4day-batch1", "TGFbeta1-8day-batch1")

# Preparar os dados para PCA
df_t_matrix_batch_effect_correction_R0_for_PCA <- prepare_pca_data(t_matrix_batch_effect_correction_R0@transcriptogramS2)
df_t_matrix_batch_effect_correction_R30_for_PCA <- prepare_pca_data(t_matrix_batch_effect_correction_R30@transcriptogramS2)

# Função para extrair as condições dos nomes das colunas
extract_condition <- function(colnames, conditions) {
  sapply(colnames, function(col) {
    match <- grep(paste(conditions, collapse = "|"), col, value = TRUE)
    if (length(match) > 0) {
      return(sub("_.*", "", match[1]))  # Extrai apenas o prefixo da condição
    } else {
      return("Unknown")  # Caso nenhuma condição seja encontrada
    }
  })
}

# Gerar labels de condições
condition_labels_R0 <- extract_condition(rownames(df_t_matrix_batch_effect_correction_R0_for_PCA), conditions)
condition_labels_R30 <- extract_condition(rownames(df_t_matrix_batch_effect_correction_R30_for_PCA), conditions)

# Garantir que os labels sejam fatores para um mapeamento correto de cores
condition_labels_R0 <- factor(condition_labels_R0, levels = conditions)
condition_labels_R30 <- factor(condition_labels_R30, levels = conditions)

# Definir paleta de cores fixa para 7 condições
color_palette <- RColorBrewer::brewer.pal(7, "Set1")
names(color_palette) <- conditions  # Mapear cores para condições

# Rodar PCA
pca_result_R0 <- run_pca(df_t_matrix_batch_effect_correction_R0_for_PCA)
pca_result_R30 <- run_pca(df_t_matrix_batch_effect_correction_R30_for_PCA)

# Gráficos pairwise para R0
plot_pairwise(
  pca_result_R0$pca_result$x, 
  condition_labels_R0, 
  pcs_to_plot = 2:10, 
  color_palette = color_palette, 
  file_name = "images/R0_pairwise_plots_patchwork.png"
)

# Gráficos pairwise para R30
plot_pairwise(
  pca_result_R30$pca_result$x, 
  condition_labels_R30, 
  pcs_to_plot = 2:10, 
  color_palette = color_palette, 
  file_name = "images/R30_pairwise_plots_patchwork.png"
)

# Função para adicionar a coluna "Day"
add_day_column <- function(pca_df, colnames_data) {
  pca_df$Day <- case_when(
    grepl("notreated-batch1", colnames_data) ~ "Day0B",
    grepl("notreated-batch2", colnames_data) ~ "Day0",
    grepl("TGFbeta1-1day-batch2", colnames_data) ~ "Day1",
    grepl("TGFbeta1-2day-batch2", colnames_data) ~ "Day2",
    grepl("TGFbeta1-3day-batch2", colnames_data) ~ "Day3",
    grepl("TGFbeta1-4day-batch1", colnames_data) ~ "Day4",
    grepl("TGFbeta1-8day-batch1", colnames_data) ~ "Day8",
    TRUE ~ "Unknown"
  )
  unmatched <- sum(pca_df$Day == "Unknown")
  if (unmatched > 0) {
    warning(unmatched, " amostras não foram mapeadas para nenhuma condição conhecida.")
  }
  return(pca_df)
}

# Adicionar coluna "Day"
pca_R0 <- pca_result_R0$pca_result$x
pca_R30 <- pca_result_R30$pca_result$x

pca_r0 <- as.data.frame(pca_R0)
pca_r0 <- add_day_column(pca_r0, rownames(df_t_matrix_batch_effect_correction_R0_for_PCA))

pca_r30 <- as.data.frame(pca_R30)
pca_r30 <- add_day_column(pca_r30, rownames(df_t_matrix_batch_effect_correction_R30_for_PCA))

# Função para calcular médias e desvios padrão por dia
calculate_mean_sd_by_day <- function(pca_df) {
  pca_df %>%
    group_by(Day) %>%
    summarise(
      across(starts_with("PC"), list(mean = ~mean(. , na.rm = TRUE), sd = ~sd(. , na.rm = TRUE))),
      .groups = 'drop'
    ) %>%
    arrange(factor(Day, levels = c("Day0B", "Day0", "Day1", "Day2", "Day3", "Day4", "Day8")))
}


mean_sd_R0 <- calculate_mean_sd_by_day(pca_r0)
mean_sd_R30 <- calculate_mean_sd_by_day(pca_r30)

# Exportar médias
write.table(mean_sd_R0, "mean_sd_R0_PCA.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(mean_sd_R30, "mean_sd_R30_PCA.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

mean_sd_R0_top10 <- mean_sd_R0[,1:21]
mean_sd_R30_top10 <- mean_sd_R30[,1:21]


##### RASTROS DE MÉDIAS E DESVIOS
color_palette <- c(
  "notreated_batch1" = "#1f77b4",
  "notreated_batch2" = "#ff7f0e",
  "TGFbeta1_1day_batch2" = "#2ca02c",
  "TGFbeta1_2day_batch2" = "#d62728",
  "TGFbeta1_3day_batch2" = "#9467bd",
  "TGFbeta1_4day_batch1" = "#8c564b",
  "TGFbeta1_8day_batch1" = "#e377c2"
)

plot_combined_pairwise_mean_sd <- function(pca_matrix, mean_sd_data, condition_labels, pcs_to_plot, color_palette, file_name) {
  library(ggplot2)
  library(patchwork)
  
  # Converter PCA matrix para data frame
  pca_df <- as.data.frame(pca_matrix)
  pca_df$Condition <- factor(condition_labels[rownames(pca_df)]) # Mapeia condições
  
  combined_plots <- list()
  
  for (i in pcs_to_plot) {
    pc_name <- paste0("PC", i)
    pc_mean_col <- paste0(pc_name, "_mean")
    pc_sd_col <- paste0(pc_name, "_sd")
    
    # Verifica se as colunas necessárias existem
    if (!all(c("PC1_mean", pc_mean_col) %in% colnames(mean_sd_data))) {
      stop(paste("Colunas necessárias para PC1 e", pc_name, "não estão presentes em mean_sd_data."))
    }
    
    combined_plots[[i]] <- ggplot() +
      # Adiciona pontos do PCA
      geom_point(data = pca_df, aes(x = PC1, y = !!sym(pc_name), color = Condition),
                 size = 1, alpha = 0.5) +
      # Adiciona pontos de média
      geom_point(data = mean_sd_data, 
                 aes(x = PC1_mean, y = !!sym(pc_mean_col), fill = Day), 
                 size = 4, shape = 21, color = "black", inherit.aes = FALSE) +
      # Adiciona barras de erro
      geom_errorbar(data = mean_sd_data, 
                    aes(x = PC1_mean, ymin = !!sym(pc_mean_col) - !!sym(pc_sd_col), 
                        ymax = !!sym(pc_mean_col) + !!sym(pc_sd_col), color = Day), 
                    width = 0.02, inherit.aes = FALSE) +
      labs(
        title = paste0("PC1 vs ", pc_name),
        x = "PC1",
        y = pc_name,
        subtitle = "Pontos menores: distribuição; Pontos maiores: médias e desvios"
      ) +
      theme_minimal() +
      scale_color_manual(values = color_palette) +
      scale_fill_manual(values = color_palette) +  # Corrige cores de preenchimento
      theme(legend.position = "right")
  }
  
  combined_plot <- patchwork::wrap_plots(combined_plots[pcs_to_plot], ncol = 3)
  ggsave(filename = file_name, plot = combined_plot, dpi = 300, width = 15, height = 10)
  message("Gráfico salvo com sucesso em: ", file_name)
}

# Exemplo de uso corrigido
plot_combined_pairwise_mean_sd(
  pca_matrix = pca_result_R0$pca_result$x,
  mean_sd_data = mean_sd_R0_top10,
  condition_labels = condition_labels_R0,
  pcs_to_plot = 2:10,
  color_palette = color_palette,
  file_name = "images/R0_combined_pairwise_mean_sd.png"
)

plot_combined_pairwise_mean_sd(
  pca_matrix = pca_result_R30$pca_result$x,
  mean_sd_data = mean_sd_R30,
  condition_labels = condition_labels_R30,
  pcs_to_plot = 2:10,
  color_palette = color_palette,
  file_name = "images/R30_combined_pairwise_mean_sd.png"
)
