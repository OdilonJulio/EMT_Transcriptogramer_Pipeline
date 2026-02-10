# normalization.R
# Normaliza os dados e gera plots em qualidade de artigo.

library(Seurat)
library(ggplot2)
library(patchwork)
library(scales)

# Carregar objeto combinado filtrado
load("filtered_combined_seurat.RData")

# Extrair as camadas e criar matriz de expressão combinada
layers <- Layers(filtered_combined_seurat)
count_matrices <- lapply(layers, function(layer_name) {
  GetAssayData(filtered_combined_seurat, assay = "RNA", layer = layer_name)
})
combined_expression_matrix <- do.call(cbind, count_matrices)
combined_expression_matrix_full <- as.matrix(combined_expression_matrix)

# --- 1. Calcular Total Count por célula (pré-normalização) ---
column_sums <- colSums(combined_expression_matrix_full)

# --- 2. Normalização Total Count ---
normalized_matrix <- sweep(combined_expression_matrix_full, 2, column_sums, FUN = "/")

# --- 3. Calcular Total Count pós-normalização (deve ser 1 por célula) ---
column_sums_normalized <- colSums(normalized_matrix)

# Salvar matriz normalizada
save(normalized_matrix, file = "normalized_matrix.RData")

# --- 4. Preparar dados para visualização ---
meta_df <- data.frame(
  cell = colnames(combined_expression_matrix_full),
  total_counts_pre = column_sums,
  total_counts_post = column_sums_normalized
)

# Adicionar metadados do Seurat (batch, sample, etc.)
meta_df$sample <- filtered_combined_seurat$orig.ident
meta_df$batch <- filtered_combined_seurat$batch


# --- 5. Plots em qualidade de artigo ---
# 5.1. Gráfico PRÉ-normalização
p_counts <- ggplot(meta_df, aes(x = total_counts_pre, fill = batch)) +
  geom_density(alpha = 0.5) +
  scale_x_log10(labels = comma) +
  theme_classic(base_size = 14) +
  labs(
    x = "Total Counts per cell (Pre-normalization, log10)", 
    y = "Density", 
    fill = "Batch"
  ) +
  ggtitle("Distribution of Total Counts Before Normalization") +
  scale_fill_manual(values = c("Batch1" = "#E69F00", "Batch2" = "#56B4E9"))

# 5.2. Gráfico PÓS-normalização (corrigido)
p_counts_norm <- ggplot(meta_df, aes(x = total_counts_post)) +
  geom_histogram(
    aes(y = after_stat(density)),
    fill = "#4C72B0",
    bins = 30,
    alpha = 0.8,
    color = "white",
    linewidth = 0.3
  ) +
  geom_density(alpha = 0.3, color = "#2E5C8A", linewidth = 0.8) +
  scale_x_continuous(limits = c(0.85, 1.05)) +  # Fechamento corrigido aqui
  theme_classic(base_size = 14) +
  labs(
    x = "Total Counts per cell (Post-normalization)", 
    y = "Density"
  ) +
  ggtitle("Distribution of Total Counts After Normalization") +
  theme(legend.position = "none")

# 5.3. Combinar e salvar
p_combined <- p_counts + p_counts_norm + 
  plot_annotation(tag_levels = "A") +
  plot_layout(widths = c(1, 1))

ggsave(
  "images/Normalization_TotalCounts.pdf",
  plot = p_combined,
  width = 14,
  height = 6,
  dpi = 300
)


# --- 6. Boxplots por amostra para checar efeito da normalização ---

# Pré-normalização
p_box_pre <- ggplot(meta_df, aes(x = sample, y = total_counts_pre, fill = batch)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  scale_y_log10(labels = comma) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) +
  labs(
    x = "Sample",
    y = "Total Counts per cell (Pre-normalization, log10)",
    fill = "Batch"
  ) +
  ggtitle("Total Counts per Sample Before Normalization") +
  scale_fill_manual(values = c("Batch1" = "#E69F00", "Batch2" = "#56B4E9"))

# Pós-normalização
p_box_post <- ggplot(meta_df, aes(x = sample, y = total_counts_post)) +
  geom_boxplot(fill = "#4C72B0", outlier.size = 0.5, alpha = 0.8) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Sample",
    y = "Total Counts per cell (Post-normalization)"
  ) +
  ggtitle("Total Counts per Sample After Normalization")

# Combinar lado a lado
p_box <- p_box_pre | p_box_post

# Exportar em alta resolução
pdf("images/Normalization_Boxplots.pdf", width = 14, height = 6, useDingbats = FALSE)
print(p_box)
dev.off()
