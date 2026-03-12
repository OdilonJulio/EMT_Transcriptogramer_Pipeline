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





##########################
# ==============================================================================
# BATCH EFFECT VALIDATION USING THE TRANSCTOGRAMER PACKAGE
# ==============================================================================

library(transcriptogramer)
library(dplyr)
library(biomaRt)
library(vroom)
library(ggplot2)
library(patchwork)
library(tidyr)
library(Ropj)

# 1. Load matrices BEFORE and AFTER correction
# (Ensure both files are in the working directory)
load("normalized_matrix.RData") # Raw Matrix (Before Batch Correction)
load("copy_matrix.RData")       # Corrected Matrix (After Batch Correction)

# 2. Isolate Day 0 columns (Control conditions)
cols_batch1 <- grep("^notreated-batch1", colnames(normalized_matrix))
cols_batch2 <- grep("^notreated-batch2", colnames(normalized_matrix))

# 3. Calculate mean expression (To optimize Transcriptogramer processing)
# Matrix BEFORE correction
mat_before <- cbind(
  Batch1 = rowMeans(normalized_matrix[, cols_batch1]),
  Batch2 = rowMeans(normalized_matrix[, cols_batch2])
)

# Matrix AFTER correction
mat_after <- cbind(
  Batch1 = rowMeans(copy_matrix[, cols_batch1]),
  Batch2 = rowMeans(copy_matrix[, cols_batch2])
)

# 4. ENSEMBL Mapping
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

dictionary <- getBM(attributes = c("ensembl_peptide_id", "ensembl_gene_id"), mart = ensembl) %>%
  mutate(ensembl_peptide_id = ifelse(ensembl_peptide_id == "", NA, ensembl_peptide_id)) %>%
  na.omit()

gene_mapping <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                      filters = "external_gene_name",
                      values = rownames(mat_before), 
                      mart = ensembl) %>%
  mutate(ensembl_gene_id = ifelse(ensembl_gene_id == "", NA, ensembl_gene_id)) %>%
  na.omit()

mapped_genes <- gene_mapping$ensembl_gene_id[match(rownames(mat_before), gene_mapping$external_gene_name)]
valid_indices <- !is.na(mapped_genes)

# Filter and rename mean matrices with ENSEMBL IDs
mat_before <- mat_before[valid_indices, , drop = FALSE]
rownames(mat_before) <- mapped_genes[valid_indices]

mat_after <- mat_after[valid_indices, , drop = FALSE]
rownames(mat_after) <- mapped_genes[valid_indices]

# 5. Load Ordering and Association matrix
ord <- vroom("ordering_HomoSapiensScore800-2024-C.txt")
assoc <- read.opj("Associationmatrix.opj")
assoc <- assoc$associationMa

assoc <- inner_join(assoc, ord, by = c("A" = "dim1")) %>% 
  dplyr::rename(protein1 = Protein) %>% 
  dplyr::inner_join(ord, by = c("B" = "dim1")) %>% 
  dplyr::rename(protein2 = Protein) %>% 
  dplyr::select(protein1, protein2)

# ==============================================================================
# 6. TRANSCRIPTOGRAMER PROCESSING (Radius = 30)
# ==============================================================================

# Object BEFORE Batch Correction
t_before <- transcriptogramPreprocess(association = assoc, ordering = ord$Protein, radius = 30)
t_before <- transcriptogramStep1(object = t_before, expression = mat_before, dictionary = dictionary)
t_before <- transcriptogramStep2(object = t_before)

# Object AFTER Batch Correction
t_after <- transcriptogramPreprocess(association = assoc, ordering = ord$Protein, radius = 30)
t_after <- transcriptogramStep1(object = t_after, expression = mat_after, dictionary = dictionary)
t_after <- transcriptogramStep2(object = t_after)

# ==============================================================================
# 7. VISUALIZATION WITH GGPLOT2 AND PATCHWORK
# ==============================================================================

# Extract smoothed data (Step 2)
df_before <- as.data.frame(t_before@transcriptogramS2)
df_after  <- as.data.frame(t_after@transcriptogramS2)

# Format to Tidy Data (Long format) forcing dplyr usage
df_before_long <- df_before %>%
  dplyr::select(Position, Batch1, Batch2) %>%
  pivot_longer(cols = c(Batch1, Batch2), names_to = "Batch", values_to = "Expression")

df_after_long <- df_after %>%
  dplyr::select(Position, Batch1, Batch2) %>%
  pivot_longer(cols = c(Batch1, Batch2), names_to = "Batch", values_to = "Expression")

# Standardized colors
color_b1 <- "#E41A1C" # Red
color_b2 <- "#377EB8" # Blue

# Plot 1: Before Correction
p_before <- ggplot(df_before_long, aes(x = Position, y = Expression, color = Batch)) +
  geom_line(alpha = 0.8, linewidth = 0.7) +
  scale_color_manual(values = c("Batch1" = color_b1, "Batch2" = color_b2),
                     labels = c("Batch1" = "Batch 1", "Batch2" = "Batch 2")) +
  labs(title = "A. Before Batch Correction",
       subtitle = "Evident basal shift (Batch Effect)",
       x = "Protein Position in the Network",
       y = "Smoothed Expression (R30)") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"),
        axis.text = element_text(color = "black"))

# Plot 2: After Correction
p_after <- ggplot(df_after_long, aes(x = Position, y = Expression, color = Batch)) +
  geom_line(alpha = 0.8, linewidth = 0.7) +
  scale_color_manual(values = c("Batch1" = color_b1, "Batch2" = color_b2),
                     labels = c("Batch1" = "Batch 1", "Batch2" = "Batch 2")) +
  labs(title = "B. After Batch Correction",
       subtitle = "Perfect alignment of control conditions",
       x = "Protein Position in the Network",
       y = "") + 
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold"),
        axis.text = element_text(color = "black"))

# Combine and Save
combined_transcriptograms <- p_before + p_after + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold"))

# Create output directory if it doesn't exist and save
if (!dir.exists("images")) dir.create("images")

ggsave(filename = "images/batch_correction_transcriptograms_R30.pdf", 
       plot = combined_transcriptograms, 
       width = 12, height = 5, dpi = 300)