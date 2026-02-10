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

load("~/mestrado-single-cell/mestrado-single-cell/t_matrix_batch_effect_correction_R0.RData")
load("~/mestrado-single-cell/mestrado-single-cell/t_matrix_batch_effect_correction_R30.RData")

df_t_matrix_batch_effect_correction_R0 <- t_matrix_batch_effect_correction_R0@transcriptogramS2
df_t_matrix_batch_effect_correction_R30 <- t_matrix_batch_effect_correction_R30@transcriptogramS2

# Verificar se todos os elementos numéricos são não negativos no primeiro dataframe
all_non_negative_R0 <- all(df_t_matrix_batch_effect_correction_R0 >= 0)
all_non_negative_R30 <- all(df_t_matrix_batch_effect_correction_R30 >= 0)

cat("Todos os valores em df_t_matrix_batch_effect_correction_R0 são não negativos?", all_non_negative_R0, "\n")
cat("Todos os valores em df_t_matrix_batch_effect_correction_R30 são não negativos?", all_non_negative_R30, "\n")

# Função para filtrar os barcodes
filter_barcodes <- function(data, conditions) {
  fixed_columns <- data[, 1:2, drop = FALSE]
  
  filtered_data <- lapply(conditions, function(cond) {
    subset <- data[, grepl(cond, colnames(data)), drop = FALSE]
    if (ncol(subset) > 20) {
      subset <- subset[, sample(ncol(subset), 20), drop = FALSE]
    }
    return(subset)
  })
  
  combined_data <- do.call(cbind, filtered_data)
  result <- cbind(fixed_columns, combined_data)
  return(result)
}

# Lista de condições para filtrar
conditions <- c("notreated-batch1", "notreated-batch2", "TGFbeta1-1day-batch2", 
                "TGFbeta1-2day-batch2", "TGFbeta1-3day-batch2", 
                "TGFbeta1-4day-batch1", "TGFbeta1-8day-batch1")

# Aplicar a função para os dois datasets
df_t_matrix_batch_effect_correction_R0_filtered <- filter_barcodes(df_t_matrix_batch_effect_correction_R0, conditions)
df_t_matrix_batch_effect_correction_R30_filtered <- filter_barcodes(df_t_matrix_batch_effect_correction_R30, conditions)


# Exportar os resultados
write.table(df_t_matrix_batch_effect_correction_R0_filtered, "df_t_matrix_batch_effect_correction_R0_filtered.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(df_t_matrix_batch_effect_correction_R30_filtered, "df_t_matrix_batch_effect_correction_R30_filtered.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Preparação para PCA
df_t_matrix_batch_effect_correction_R0_for_PCA <- df_t_matrix_batch_effect_correction_R0[,-2] 
df_t_matrix_batch_effect_correction_R30_for_PCA <- df_t_matrix_batch_effect_correction_R30[,-2]

# Definir nomes das linhas
rownames(df_t_matrix_batch_effect_correction_R0_for_PCA) <- df_t_matrix_batch_effect_correction_R0_for_PCA$Protein
df_t_matrix_batch_effect_correction_R0_for_PCA <- df_t_matrix_batch_effect_correction_R0_for_PCA[, -1]

rownames(df_t_matrix_batch_effect_correction_R30_for_PCA) <- df_t_matrix_batch_effect_correction_R30_for_PCA$Protein
df_t_matrix_batch_effect_correction_R30_for_PCA <- df_t_matrix_batch_effect_correction_R30_for_PCA[, -1]

# Função para realizar PCA
perform_pca <- function(data) {
  pca_result <- prcomp(t(data), scale. = FALSE)
  pca_df <- as.data.frame(pca_result$x)
  return(pca_df)
}

# Realizar PCA
pca_R0 <- perform_pca(df_t_matrix_batch_effect_correction_R0_for_PCA)
pca_R30 <- perform_pca(df_t_matrix_batch_effect_correction_R30_for_PCA)

# Adicionar a coluna "Day"
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
  return(pca_df)
}

pca_r0 <- add_day_column(pca_R0, colnames(df_t_matrix_batch_effect_correction_R0_for_PCA))
pca_r30 <- add_day_column(pca_R30, colnames(df_t_matrix_batch_effect_correction_R30_for_PCA))

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
