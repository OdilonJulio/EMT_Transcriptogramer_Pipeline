################################################################################
# TÍTULO: EMT_22_publication_figure.R
# DESCRIÇÃO: Gera a Figura Principal (4 painéis) em qualidade de publicação,
#            ilustrando o filtro topológico, a clusterização, 
#            os marcadores diferenciais (Volcano) e a onda do transcriptograma.
################################################################################

suppressMessages({
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(dplyr)
  library(transcriptogramer)
})

# Cria diretório para a figura final
dir.create("images/publication", showWarnings = FALSE)

# ==========================================
# 1. CARREGAR DADOS DOS SCRIPTS ANTERIORES
# ==========================================
# Carrega as PCAs (Bruta e Suavizada)
load("pca_result_R0.RData")
load("pca_result_R30.RData")

# Carrega os marcadores calculados no EMT_21
library(openxlsx)
markers_df <- read.xlsx("markers_analysis/EMT_Cell_and_Outlier_Markers.xlsx", sheet = "Marcadores_Estagios_EMT")

# Carrega as matrizes do Transcriptograma
load("t_matrix_R0.RData")
load("t_matrix_R30.RData")

# ==========================================
# 2. PAINEL A e B: ESPAÇO DE DIMENSÕES (PCA/t-SNE)
# ==========================================
# Extrair PC1 e PC2 do R0 (Controle / Ruído)
df_pca_r0 <- data.frame(
  Dim1 = pca_result_R0$pca_result$x[, 1],
  Dim2 = pca_result_R0$pca_result$x[, 2],
  Condition = "Raw (R0)"
)

# Extrair PC1 e PC2 do R30 (Suavizado / Topológico)
df_pca_r30 <- data.frame(
  Dim1 = pca_result_R30$pca_result$x[, 1],
  Dim2 = pca_result_R30$pca_result$x[, 2],
  Condition = "Smoothed (R30)"
)

# Vamos simular a cor pelo estágio usando o PC1 como fizemos no EMT_21
quantiles_r30 <- quantile(df_pca_r30$Dim1, probs = c(0, 0.33, 0.66, 1))
df_pca_r0$Stage <- case_when(
  df_pca_r30$Dim1 <= quantiles_r30[2] ~ "Epithelial",
  df_pca_r30$Dim1 > quantiles_r30[2] & df_pca_r30$Dim1 <= quantiles_r30[3] ~ "Transition",
  TRUE ~ "Mesenchymal"
)
df_pca_r30$Stage <- df_pca_r0$Stage

# Cores padronizadas para o artigo
stage_colors <- c("Epithelial" = "#4575b4", "Transition" = "#fdae61", "Mesenchymal" = "#d73027")

pA <- ggplot(df_pca_r0, aes(x = Dim1, y = Dim2, color = Stage)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = stage_colors) +
  labs(title = "A. Cell State (Raw R0)", x = "Component 1", y = "Component 2") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

pB <- ggplot(df_pca_r30, aes(x = Dim1, y = Dim2, color = Stage)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = stage_colors) +
  labs(title = "B. Topological Smoothing (R30)", x = "Component 1", y = "Component 2") +
  theme_classic(base_size = 14) +
  theme(legend.position = "right", plot.title = element_text(face = "bold"))

# ==========================================
# 3. PAINEL C: VOLCANO PLOT (Metabolic Switch)
# ==========================================
# Preparar os dados para o Volcano
volcano_data <- markers_df %>%
  mutate(
    Significance = -log10(p_val_adj + 1e-300), # Evita log(0)
    # Define categorias para colorir
    Status = case_when(
      avg_log2FC > 0.5 & p_val_adj < 0.05 ~ "Upregulated",
      avg_log2FC < -0.5 & p_val_adj < 0.05 ~ "Downregulated",
      TRUE ~ "Not Sig"
    )
  )

# Selecionar genes chave para rotular no gráfico (Ajuste com os seus melhores achados)
genes_to_label <- c("CDH1", "VIM", "FN1", "ATP5B", "NDUFB8", "MT1X", "MT2A", "SNAI1", "ZEB1")
labels_data <- volcano_data %>% filter(gene %in% genes_to_label)

pC <- ggplot(volcano_data, aes(x = avg_log2FC, y = Significance, color = Status)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "#d73027", "Downregulated" = "#4575b4", "Not Sig" = "grey80")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  # Adiciona os nomes dos genes sem sobrepor os pontos
  geom_text_repel(data = labels_data, aes(label = gene), 
                  size = 4, fontface = "italic", color = "black", 
                  box.padding = 0.5, max.overlaps = Inf) +
  labs(title = "C. Differential Expression", x = expression(log[2]("Fold Change")), y = expression(-log[10]("P-value"))) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

# ==========================================
# 4. PAINEL D: ONDA DO TRANSCRIPTOGRAMA
# ==========================================
# Extrair o perfil médio de uma célula mesenquimal extrema para mostrar a "onda"
mesen_cells <- rownames(df_pca_r30 %>% filter(Stage == "Mesenchymal"))
mat_R0_num <- as.matrix(t_matrix_R0@transcriptogramS2[, -(1:2)])
mat_R30_num <- as.matrix(t_matrix_R30@transcriptogramS2[, -(1:2)])

# Pega a média das células avançadas
wave_R0 <- rowMeans(mat_R0_num[, mesen_cells[1:50]], na.rm = TRUE)
wave_R30 <- rowMeans(mat_R30_num[, mesen_cells[1:50]], na.rm = TRUE)
pos <- t_matrix_R0@transcriptogramS2$Position

df_wave <- data.frame(Position = pos, R0 = wave_R0, R30 = wave_R30)

pD <- ggplot(df_wave, aes(x = Position)) +
  geom_line(aes(y = R0, color = "Control (R0)"), alpha = 0.5, linewidth = 0.5) +
  geom_line(aes(y = R30, color = "Metabolic Switch (R30)"), alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("Control (R0)" = "#4575b4", "Metabolic Switch (R30)" = "#d73027")) +
  # Destacar a região do cluster metabólico (exemplo: posições 2500 a 3500)
  annotate("rect", xmin = 2500, xmax = 3500, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "#fdae61") +
  annotate("text", x = 3000, y = max(wave_R30)*0.9, label = "OXPHOS / Metabolism", fontface = "bold", size = 4) +
  labs(title = "D. Transcriptogram: Metabolic Reprogramming", 
       x = "Protein-Protein Interaction (Network Position)", y = "Relative Expression", color = "") +
  theme_classic(base_size = 14) +
  theme(legend.position = c(0.8, 0.9), plot.title = element_text(face = "bold"))

# ==========================================
# 5. MONTAGEM FINAL COM PATCHWORK
# ==========================================
# O layout (pA ao lado do pB) em cima de (pC ao lado do pD)
final_plot <- (pA | pB) / (pC | pD)

# Salvar em alta resolução (TIFF e PDF) exigido pela IEEE
ggsave("images/publication/Figure_Main_Panel.pdf", final_plot, width = 16, height = 12, dpi = 300)
ggsave("images/publication/Figure_Main_Panel.tiff", final_plot, width = 16, height = 12, dpi = 300, compression = "lzw")

message("Figura gerada com sucesso! Verifique a pasta 'images/publication'.")