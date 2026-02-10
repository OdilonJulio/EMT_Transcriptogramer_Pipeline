# batch_effect_correction.R
# Realiza a análise de transcriptograma.

# Carregar matriz normalizada
load("normalized_matrix.RData")

# Filtrar colunas para "notreated_batch1" e "notreated_batch2"
cols_batch1 <- grep("^notreated-batch1", colnames(normalized_matrix))
cols_batch2 <- grep("^notreated-batch2", colnames(normalized_matrix))

# Submatrizes para cada batch
sub_matrix_batch1 <- normalized_matrix[, cols_batch1]
sub_matrix_batch2 <- normalized_matrix[, cols_batch2]

# Calcular as médias por linha
mean_batch1 <- rowMeans(sub_matrix_batch1)
mean_batch2 <- rowMeans(sub_matrix_batch2)

# Criar um data frame com os resultados
result <- data.frame(
  Gene = rownames(normalized_matrix),
  Mean_Batch1 = mean_batch1,
  Mean_Batch2 = mean_batch2,
  Ratio_Batch2_Batch1 = ifelse(mean_batch1 == 0, 0, mean_batch2 / mean_batch1) # Divisão das médias
)

## Multiplicando todas as colunas de normalized_matrix que possuem "batch1" no colname.

colnames_batch1 <- grep("batch1", colnames(normalized_matrix))

# criando cópia de normalized_matrix

copy_matrix <- normalized_matrix

# Multiplicando a razão pelas colunas selecionadas.

copy_matrix[, colnames_batch1] <- sweep(normalized_matrix[, colnames_batch1], 1, result$Ratio_Batch2_Batch1, "*")


# Salvar matriz com efeito de lote corrigido.
save(copy_matrix, file = "copy_matrix.RData")


## VISUALIZAÇÃO


library(ggplot2)
library(ggpointdensity)
library(patchwork)
library(viridis)

# Load data (all untreated cells)
batch1 <- grep("notreated-batch1", colnames(copy_matrix), value = TRUE)
batch2 <- grep("notreated-batch2", colnames(copy_matrix), value = TRUE)

# 1. Density Plot (All Cells)
p_density <- ggplot(data.frame(
  Expression = c(copy_matrix[, c(batch1, batch2)]),
  Batch = ifelse(grepl("batch1", colnames(copy_matrix[, c(batch1, batch2)])), 
                 "Batch1", "Batch2")
), aes(x = Expression, fill = Batch)) +
  geom_density(alpha = 0.3, adjust = 0.8, linewidth = 0.7) +  # Aumente alpha e reduza adjust
  scale_x_log10(limits = c(1e-3, NA)) +
  scale_fill_manual(values = c("#E41A1C80", "#377EB880")) +   # Cores com transparência (hex +80)
  labs(title = "Expression Distribution (Untreated Cells)",
       x = expression(log[10]("Expression Level")), 
       y = "Density") +
  theme_classic(base_size = 12) +
  theme(legend.position = c(0.8, 0.8))
# 2. Scatter Plot (Gene Means)
mean_df <- data.frame(
  Mean_Batch1 = rowMeans(copy_matrix[, batch1]),
  Mean_Batch2 = rowMeans(copy_matrix[, batch2])
)

# 2. Scatter Plot (Gene Means) - VERSÃO CORRIGIDA
mean_df <- data.frame(
  Mean_Batch1 = rowMeans(copy_matrix[, batch1]),
  Mean_Batch2 = rowMeans(copy_matrix[, batch2])
)

p_scatter <- ggplot(mean_df, aes(x = Mean_Batch1, y = Mean_Batch2)) +
  {
    if(requireNamespace("ggpointdensity", quietly = TRUE)){
      ggpointdensity::geom_pointdensity(size = 0.3, alpha = 0.8)
    } else {
      geom_point(alpha = 0.3, color = "gray30")
    }
  } +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  scale_x_log10(limits = c(1e-3, NA)) +
  scale_y_log10(limits = c(1e-3, NA)) +
  {
    if(requireNamespace("viridis", quietly = TRUE)){
      scale_color_viridis(option = "plasma")
    } else {
      scale_color_gradient(low = "blue", high = "red")
    }
  } +
  labs(
    title = "Gene Expression Means Comparison",
    x = "Batch1 Mean Expression (log10)",
    y = "Batch2 Mean Expression (log10)",
    color = "Point Density"
  ) +
  theme_minimal(base_size = 12)

# 3. Combine plots side-by-side
combined_plot <- p_density + p_scatter + 
  plot_annotation(tag_levels = 'A') +
  plot_layout(widths = c(1, 1.2))

# Save
ggsave("images/batch_correction_validation.pdf", combined_plot, 
       width = 12, height = 5, dpi = 300)


library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(ggrepel)
library(stringr)  # For str_extract()

# Prepare metadata
metadata <- data.frame(
  Cell = colnames(copy_matrix),
  Batch = ifelse(grepl("batch1", colnames(copy_matrix)), "Batch1", "Batch2"),
  Day = str_extract(colnames(copy_matrix), "day[0-9]+"),
  TotalCounts = colSums(copy_matrix)
)

# 1. Enhanced Violin Plot by Day
p_violin <- ggplot(metadata, aes(x = Day, y = TotalCounts, fill = Batch)) +
  geom_violin(scale = "width", alpha = 0.7, width = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.5, outlier.size = 0.5) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
  scale_y_log10() +
  labs(title = "Total Counts Distribution by Treatment Day",
       y = "Total UMI Counts (log10)",
       x = "Treatment Day") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 2. Enhanced MDS Plot (with sample labeling)
set.seed(123)
sampled_cells <- sample(ncol(copy_matrix), min(1000, ncol(copy_matrix)))
mds_data <- cmdscale(dist(t(copy_matrix[, sampled_cells])))
mds_df <- data.frame(
  MDS1 = mds_data[,1],
  MDS2 = mds_data[,2],
  Batch = metadata$Batch[sampled_cells],
  Day = metadata$Day[sampled_cells],
  Sample = colnames(copy_matrix)[sampled_cells]
)

p_mds <- ggplot(mds_df, aes(x = MDS1, y = MDS2, color = Day, shape = Batch)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = Day), level = 0.8, linetype = 2) +
  scale_color_brewer(palette = "Set2") +
  scale_shape_manual(values = c(16, 17)) +
  labs(title = "Global Sample Similarity (MDS Plot)",
       x = "MDS Dimension 1",
       y = "MDS Dimension 2",
       color = "Treatment Day",
       shape = "Batch") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right") +
  guides(color = guide_legend(ncol = 2))

# 3. Combined Plot Layout
full_plot <- (p_violin + p_mds) +
  plot_layout(widths = c(1, 1.2)) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold"))

# Save high-quality output
ggsave("images/treatment_days_analysis.pdf", full_plot, 
       width = 14, height = 6, dpi = 300)
