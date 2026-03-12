# pca_analysis_enhanced.R
# Versão final com ajustes para dissertação/artigo

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)
library(transcriptogramer)
# -----------------------------
# 1. Carregar objetos
# -----------------------------
load("t_matrix_R0.RData")
load("t_matrix_R30.RData")

# -----------------------------
# 2. Funções PCA
# -----------------------------



prepare_pca_data <- function(transcriptogram_obj) {
  df <- transcriptogram_obj@transcriptogramS2
  df <- df[, !colnames(df) %in% c("Protein", "Position")]
  if (any(is.na(df))) stop("O dataset contém valores NA.")
  return(t(df))
}

run_pca_enhanced <- function(pca_input) {
  pca_result <- prcomp(pca_input, center = TRUE, scale. = FALSE)
  eigenvalues <- pca_result$sdev^2
  explained_variance <- eigenvalues / sum(eigenvalues)
  cumulative_variance <- cumsum(explained_variance)
  return(list(
    pca_result = pca_result,
    eigenvalues = eigenvalues,
    explained_variance = explained_variance,
    cumulative_variance = cumulative_variance
  ))
}

extract_condition <- function(colnames) {
  conditions <- c("notreated-batch1", "notreated-batch2", "TGFbeta1-1day-batch2",
                  "TGFbeta1-2day-batch2", "TGFbeta1-3day-batch2",
                  "TGFbeta1-4day-batch1", "TGFbeta1-8day-batch1")
  sapply(colnames, function(col) {
    match <- grep(paste(conditions, collapse = "|"), col, value = TRUE)
    if (length(match) > 0) sub("_.*", "", match[1]) else "Unknown"
  })
}

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
  if (any(pca_df$Day == "Unknown")) warning(sum(pca_df$Day == "Unknown"), " amostras não mapeadas")
  return(pca_df)
}

# -----------------------------
# 3. Preparar dados
# -----------------------------
df_R0 <- prepare_pca_data(t_matrix_R0)
df_R30 <- prepare_pca_data(t_matrix_R30)

condition_labels_R0 <- extract_condition(rownames(df_R0))
condition_labels_R30 <- extract_condition(rownames(df_R30))

color_palette <- c(
  "notreated-batch1" = "#984EA3",
  "notreated-batch2" = "#377EB8",
  "TGFbeta1-1day-batch2" = "#4DAF4A",
  "TGFbeta1-2day-batch2" = "#3A5400",
  "TGFbeta1-3day-batch2" = "#A65628",
  "TGFbeta1-4day-batch1" = "#FF7F00",
  "TGFbeta1-8day-batch1" = "#E41A1C"
)

# -----------------------------
# 4. PCA + metadados
# -----------------------------
pca_result_R0 <- run_pca_enhanced(df_R0)
pca_result_R30 <- run_pca_enhanced(df_R30)

pca_r0_df <- add_day_column(as.data.frame(pca_result_R0$pca_result$x), rownames(df_R0))
pca_r30_df <- add_day_column(as.data.frame(pca_result_R30$pca_result$x), rownames(df_R30))

# -----------------------------
# 5. Calcular centroides
# -----------------------------
centroids_r0 <- pca_r0_df %>%
  group_by(Day) %>%
  summarize(across(starts_with("PC"), mean), .groups = "drop")

centroids_r30 <- pca_r30_df %>%
  group_by(Day) %>%
  summarize(across(starts_with("PC"), mean), .groups = "drop")

# -----------------------------
# 6. Plots com setas temporais
# -----------------------------
temporal_order <- c(
  "notreated-batch1", "notreated-batch2", "TGFbeta1-1day-batch2",
  "TGFbeta1-2day-batch2", "TGFbeta1-3day-batch2",
  "TGFbeta1-4day-batch1", "TGFbeta1-8day-batch1"
)

# -----------------------------
# 6a. Plots com setas temporais (PC1 vs PC2) com escalas iguais
# -----------------------------
generate_temporal_pca_plot <- function(pca_df, centroids, title, limits) {
  centroids_ordered <- centroids %>%
    mutate(Day = factor(Day, levels = temporal_order)) %>%
    arrange(Day)
  
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Day)) +
    geom_point(size = 0.5) +
    geom_point(data = centroids_ordered, size = 3, shape = 18) +
    geom_path(data = centroids_ordered, aes(group = 1),
              arrow = arrow(length = grid::unit(0.3, "cm")),
              color = "black", linewidth = 0.6) +
    scale_color_manual(values = color_palette) +
    coord_equal(xlim = limits$x, ylim = limits$y) +   # <-- mesma graduação
    ggtitle(title) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

# definir limites globais
xlim_all <- range(c(pca_r0_df$PC1, pca_r30_df$PC1))
ylim_all <- range(c(pca_r0_df$PC2, pca_r30_df$PC2))
limits <- list(x = xlim_all, y = ylim_all)

plot_R0 <- generate_temporal_pca_plot(pca_r0_df, centroids_r0, "PCA - Condition R0", limits)
plot_R30 <- generate_temporal_pca_plot(pca_r30_df, centroids_r30, "PCA - Condition R30", limits)

combined_plots <- plot_R0 + plot_R30 + plot_layout(guides = 'collect')
ggsave("images/pca_combined_plot_temporal.pdf", combined_plots, width = 12, height = 6, dpi = 300)


# -----------------------------
# 6a-extra. Versão única com os dois quadros (topo: mesma escala; baixo: escalas independentes)
# -----------------------------

# Função para a versão com escalas independentes (sem coord_equal)
generate_temporal_pca_plot_indep <- function(pca_df, centroids, title) {
  centroids_ordered <- centroids %>%
    mutate(Day = factor(Day, levels = temporal_order)) %>%
    arrange(Day)
  
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Day)) +
    geom_point(size = 0.5) +
    geom_point(data = centroids_ordered, size = 3, shape = 18) +
    geom_path(data = centroids_ordered, aes(group = 1),
              arrow = arrow(length = grid::unit(0.3, "cm")),
              color = "black", linewidth = 0.6) +
    scale_color_manual(values = color_palette) +
    ggtitle(title) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

# Linha de cima (mesma escala)
row_same_scale <- (plot_R0 + plot_R30) + 
  plot_annotation(title = "(A) Same scale")

# Linha de baixo (escalas independentes)
plot_R0_indep  <- generate_temporal_pca_plot_indep(pca_r0_df,  centroids_r0,  "Condition R0")
plot_R30_indep <- generate_temporal_pca_plot_indep(pca_r30_df, centroids_r30, "Condition R30")

row_indep_scale <- (plot_R0_indep + plot_R30_indep) + 
  plot_annotation(title = "(B) Independent scales")

# Junta tudo em um único arquivo, com título geral e legenda única
both_scales_plot <- (
  row_same_scale / row_indep_scale
) + plot_annotation(
  title = "Temporal PCA comparison"
) & theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Recolhe as legendas e coloca embaixo
both_scales_plot <- both_scales_plot + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("images/pca_combined_plot_temporal_both_scales.pdf",
       both_scales_plot, width = 12, height = 12, dpi = 300)

# Junta tudo em um único arquivo, com título geral e legenda única
indep_scales_plot <- row_indep_scale + plot_annotation(
  title = "Temporal PCA comparison"
) & theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Recolhe as legendas e coloca embaixo
indep_scales_plot <- indep_scales_plot + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("images/pca_combined_plot_temporal_indep_scales.pdf",
       both_scales_plot, width = 12, height = 12, dpi = 300)

# -----------------------------
# 6b. Novo grid 2x5 (PC1×PC2..6, R0 na linha de cima, R30 na linha de baixo)
# -----------------------------

# Definir limites globais por condição
limits_r0 <- list(
  x = range(pca_r0_df$PC1),
  y = range(pca_r0_df[, paste0("PC", 2:6)])
)

limits_r30 <- list(
  x = range(pca_r30_df$PC1),
  y = range(pca_r30_df[, paste0("PC", 2:6)])
)

# Alongar o eixo Y (PC1) multiplicando a amplitude
stretch_factor <- 0.8
mid_y_r0  <- mean(limits_r0$y)
range_y_r0 <- diff(limits_r0$y) * stretch_factor / 2
limits_r0$y <- c(mid_y_r0 - range_y_r0, mid_y_r0 + range_y_r0)

mid_y_r30  <- mean(limits_r30$y)
range_y_r30 <- diff(limits_r30$y) * stretch_factor / 2
limits_r30$y <- c(mid_y_r30 - range_y_r30, mid_y_r30 + range_y_r30)

# Função para gerar plots
generate_pair_plot <- function(pca_df, pc_y, title, limits) {
  ggplot(pca_df, aes(x = PC1, y = .data[[pc_y]], color = Day)) +
    geom_point(size = 0.2) +
    scale_color_manual(values = color_palette) +
    labs(title = title, x = "PC1", y = pc_y) +
    coord_cartesian(xlim = limits$x, ylim = limits$y) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      axis.text = element_text(size = 6),
      axis.title = element_text(size = 8)
    )
}

# Criar lista de plots
plots_list <- list()
for (pc in paste0("PC", 2:6)) {
  plots_list[[paste0("R0_", pc)]]  <- generate_pair_plot(pca_r0_df,  pc, paste("R0:",  "PC1 vs", pc), limits_r0)
  plots_list[[paste0("R30_", pc)]] <- generate_pair_plot(pca_r30_df, pc, paste("R30:", "PC1 vs", pc), limits_r30)
}


# Montar grid 2x5 (linha 1 = R0, linha 2 = R30)
grid_plot <- (
  plots_list$R0_PC2 + plots_list$R0_PC3 + plots_list$R0_PC4 + plots_list$R0_PC5 + plots_list$R0_PC6 +
    plots_list$R30_PC2 + plots_list$R30_PC3 + plots_list$R30_PC4 + plots_list$R30_PC5 + plots_list$R30_PC6
) +
  plot_layout(ncol = 5, guides = "collect") &
  theme(legend.position = "bottom")

# Salvar em PDF largo
ggsave("images/pca_grid_PC1_vs_PC2to6.pdf",
       grid_plot, width = 20, height = 10, dpi = 300)





# -----------------------------
# 7. Gráficos de variância PCA
# -----------------------------
find_elbow_point <- function(cum_var) {
  n <- length(cum_var)
  coords <- data.frame(x = 1:n, y = cum_var)
  line_start <- coords[1, ]
  line_end <- coords[n, ]
  
  distances <- sapply(1:n, function(i) {
    pt <- coords[i, ]
    abs((line_end$y - line_start$y) * pt$x - (line_end$x - line_start$x) * pt$y +
          line_end$x * line_start$y - line_end$y * line_start$x) /
      sqrt((line_end$y - line_start$y)^2 + (line_end$x - line_start$x)^2)
  })
  
  which.max(distances)
}

plot_combined_variance <- function(pca_result, n_pcs, condition_name = "", label_all = FALSE) {
  total_pcs <- length(pca_result$explained_variance)
  label <- if (label_all) "All PCs" else paste(n_pcs, "PCs")
  
  var_df <- data.frame(
    PC = 1:n_pcs,
    Explained = pca_result$explained_variance[1:n_pcs] * 100,
    Cumulative = pca_result$cumulative_variance[1:n_pcs] * 100
  )
  
  elbow_pc <- find_elbow_point(var_df$Cumulative)
  elbow_y <- var_df$Cumulative[elbow_pc]
  
  ggplot(var_df, aes(x = PC)) +
    geom_bar(aes(y = Explained, fill = "Explained Variance"),
             stat = "identity", width = 0.7, alpha = 0.8) +
    geom_line(aes(y = Cumulative, color = "Cumulative Variance"), linewidth = 1) +
    geom_point(aes(y = Cumulative, color = "Cumulative Variance"), size = 2) +
    geom_vline(xintercept = elbow_pc, linetype = "dotdash", color = "darkgreen") +
    annotate("text", x = elbow_pc, y = 105,
             label = sprintf("Elbow: PC%d (%.1f%%)", elbow_pc, elbow_y),
             color = "darkgreen", size = 4, vjust = 0) +
    geom_hline(yintercept = 80, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = 90, linetype = "dashed", color = "gray50") +
    scale_y_continuous(limits = c(0, 105), labels = percent_format(scale = 1)) +
    scale_fill_manual(values = c("Explained Variance" = "#1f78b4")) +
    scale_color_manual(values = c("Cumulative Variance" = "#e31a1c")) +
    labs(title = sprintf("PCA Variance - %s (%s)", condition_name, label),
         x = "Principal Component", y = "Percentage of Variance (%)",
         fill = "", color = "") +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.border = element_rect(fill = NA, color = "gray80")
    )
}

generate_all_plots <- function(pca_result, condition_name) {
  n_total <- length(pca_result$explained_variance)
  plot_100 <- plot_combined_variance(pca_result, 100, condition_name, label_all = FALSE)
  plot_all <- plot_combined_variance(pca_result, n_total, condition_name, label_all = TRUE)
  
  combined <- wrap_plots(plot_100, plot_all, ncol = 1) +
    plot_annotation(title = sprintf("PCA Variance Decomposition - %s", condition_name),
                    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))
  
  output_file <- sprintf("images/pca_variance_%s.pdf", gsub(" ", "_", tolower(condition_name)))
  ggsave(output_file, combined, width = 10, height = 12, dpi = 300)
  return(combined)
}

final_plot_R0 <- generate_all_plots(pca_result_R0, "Condition R0")
final_plot_R30 <- generate_all_plots(pca_result_R30, "Condition R30")

# Combinar ambas
all_conditions_plot <- wrap_plots(final_plot_R0, final_plot_R30, ncol = 2) +
  plot_annotation(title = "Comparative PCA Variance Analysis",
                  theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)))

ggsave("images/pca_variance_complete_analysis.pdf", all_conditions_plot, 
       width = 20, height = 12, dpi = 600)
# -----------------------------
# 8. Day 8 Heterogeneity Analysis (Endpoint) - BIG FONTS VERSION
# -----------------------------

generate_day8_heterogeneity_plot <- function(pca_df, title) {
  
  # Define subsets
  day8_data <- subset(pca_df, Day == "TGFbeta1-8day-batch1")
  background_data <- subset(pca_df, Day != "TGFbeta1-8day-batch1")
  
  day8_color <- "#E41A1C" 
  x_limits <- range(pca_df$PC1)
  y_limits <- range(pca_df$PC2)
  
  p <- ggplot() +
    # Background
    geom_point(data = background_data, aes(x = PC1, y = PC2), 
               color = "grey85", size = 1.0, alpha = 0.4) + # Aumentei levemente o ponto também (0.5 -> 1.0)
    
    # Foreground (Day 8)
    geom_point(data = day8_data, aes(x = PC1, y = PC2), 
               color = day8_color, size = 2.5, alpha = 0.8) + # Aumentei o ponto (1.2 -> 2.5)
    
    # Rug plot
    geom_rug(data = day8_data, aes(x = PC1), 
             sides = "b", color = day8_color, alpha = 0.5, length = unit(0.3, "cm")) +
    
    # Styling (BIG FONTS)
    labs(title = title,
         subtitle = "Highlight: Day 8 Cells (Red) over trajectory (Grey)",
         x = "PC1 (Pseudotime / EMT State)", 
         y = "PC2") +
    coord_cartesian(xlim = x_limits, ylim = y_limits) +
    
    # --- AQUI ESTÁ A MUDANÇA CRÍTICA ---
    theme_minimal(base_size = 24) + # Dobrei de 12 para 24
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5), # Herda o tamanho 24+
      plot.subtitle = element_text(size = 20, hjust = 0.5, color = "grey40"), # Dobrei de 10 para 20
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey80", linewidth = 1) # Borda mais grossa
    )
  
  return(p)
}

# Generate plots
plot_day8_R0 <- generate_day8_heterogeneity_plot(pca_r0_df, "Heterogeneity at Day 8 (R0)")
plot_day8_R30 <- generate_day8_heterogeneity_plot(pca_r30_df, "Heterogeneity at Day 8 (R30)")

# Combine
day8_combined_plot <- (plot_day8_R0 + plot_day8_R30) +
  plot_annotation(
    title = "Persistence of Hybrid States at Endpoint (Day 8)",
    theme = theme(plot.title = element_text(size = 32, face = "bold", hjust = 0.5)) # Dobrei de 16 para 32
  )

# Save (Aumentei a altura da imagem para caber a fonte grande)
ggsave("images/day8_heterogeneity_analysis.pdf", 
       day8_combined_plot, width = 16, height = 9, dpi = 300)

message("Gráfico com FONTES GIGANTES salvo com sucesso.")


# ==============================================================================
# 9. FINAL GRID 1x5: PC1 vs PC2...PC6 (R30 Only) - ARTICLE QUALITY JPG
# ==============================================================================

library(dplyr)
library(ggplot2)
library(patchwork) 

# 1. Prepare Data with Explicit Temporal Order
pca_r30_df_ordered <- pca_r30_df %>%
  mutate(Day = factor(Day, levels = temporal_order)) %>%
  arrange(Day)

centroids_ordered <- pca_r30_df_ordered %>%
  group_by(Day) %>%
  summarize(across(starts_with("PC"), mean), .groups = "drop")

# Global limits for the Y-axis (PC2 to PC6) to ensure comparability
y_cols <- c("PC2", "PC3", "PC4", "PC5", "PC6")
global_y_limits <- range(pca_r30_df_ordered[, y_cols])
y_margin <- diff(global_y_limits) * 0.05
global_y_limits[1] <- global_y_limits[1] - y_margin
global_y_limits[2] <- global_y_limits[2] + y_margin

# 2. Define Plotting Function
create_consistent_traj_plot_1x5 <- function(y_pc_name, y_pc_idx, shared_ylim) {
  var_x <- round(pca_result_R30$explained_variance[1] * 100, 1)
  var_y <- round(pca_result_R30$explained_variance[y_pc_idx] * 100, 1)
  
  p <- ggplot(pca_r30_df_ordered, aes(x = PC1, y = .data[[y_pc_name]], color = Day)) +
    
    # REDUCED POINT SIZE AND ALPHA FOR 13K+ CELLS
    # Creates a transparent "cloud" rather than a solid block of color
    geom_point(size = 0.6, alpha = 0.5) + 
    
    # Centroids and trajectory path (kept large and solid for visibility)
    geom_point(data = centroids_ordered, aes(x = PC1, y = .data[[y_pc_name]]), 
               size = 3, shape = 18, color = "black") + 
    geom_path(data = centroids_ordered, aes(x = PC1, y = .data[[y_pc_name]], group = 1),
              arrow = arrow(length = unit(0.4, "cm")),
              color = "black", linewidth = 0.8) +
    
    scale_color_manual(values = color_palette) +
    coord_cartesian(ylim = shared_ylim) + 
    
    labs(x = paste0("PC1 (", var_x, "%)"), 
         y = paste0(y_pc_name, " (", var_y, "%)")) +
    theme_minimal(base_size = 20) +
    theme(
      panel.border = element_rect(color = "grey60", fill = NA, linewidth = 1),
      axis.title = element_text(face = "bold", size = 22),
      axis.text = element_text(size = 16, color = "black"),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

# 3. Generate the 5 Plots
plot_pc2 <- create_consistent_traj_plot_1x5("PC2", 2, global_y_limits)
plot_pc3 <- create_consistent_traj_plot_1x5("PC3", 3, global_y_limits)
plot_pc4 <- create_consistent_traj_plot_1x5("PC4", 4, global_y_limits)
plot_pc5 <- create_consistent_traj_plot_1x5("PC5", 5, global_y_limits)
plot_pc6 <- create_consistent_traj_plot_1x5("PC6", 6, global_y_limits)

# 4. Assemble Grid 1x5 using Patchwork
final_1x5_grid <- (plot_pc2 | plot_pc3 | plot_pc4 | plot_pc5 | plot_pc6) +
  plot_layout(guides = "collect") & 
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_blank(), # Removes legend title to save space
    legend.text = element_text(size = 22, face = "bold"),
    legend.key.size = unit(1.5, "cm") 
  ) &
  # INCREASES LEGEND DOT SIZE AND OPACITY OVERRIDING THE GEOM_POINT SETTINGS
  guides(color = guide_legend(override.aes = list(size = 10, alpha = 1)))

# 5. Save as Publication Quality JPG (Banner format)
# Fixed the object name to final_1x5_grid
ggsave("images/pca_grid_PC1_vs_PC2to6_R30.jpg", plot = final_1x5_grid, 
       device = "jpeg", width = 30, height = 8, dpi = 600)

message("1x5 Grid JPG successfully saved: images/pca_grid_PC1_vs_PC2to6_R30.jpg")

# (Aqui entra a Seção 10 do seu script inalterada)

# ==============================================================================
# 10. R30 VARIANCE ANALYSIS (FOCUSED 35 PCs vs ALL PCs) - ARTICLE QUALITY
# ==============================================================================

library(ggplot2)
library(patchwork)
library(scales)

# 1. Configuração de Tema com Fontes Aumentadas (Article Style)
# ------------------------------------------------------------------------------
theme_variance_article <- theme_minimal(base_size = 20) + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 26),
    plot.subtitle = element_text(hjust = 0.5, size = 18, color = "grey40"),
    axis.title = element_text(face = "bold", size = 20),
    axis.text = element_text(size = 16, color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 16),
    legend.title = element_text(face = "bold", size = 18),
    panel.border = element_rect(color = "grey60", fill = NA, size = 1),
    panel.grid.minor = element_blank()
  )

# 2. Função Geradora Customizada (R30 Variance)
# ------------------------------------------------------------------------------
create_variance_plot_custom <- function(pca_res, n_limit, title_text, show_elbow = TRUE) {
  
  # Preparar dados
  total_vars <- length(pca_res$explained_variance)
  limit <- min(n_limit, total_vars)
  
  df <- data.frame(
    PC = 1:limit,
    Explained = pca_res$explained_variance[1:limit] * 100,
    Cumulative = pca_res$cumulative_variance[1:limit] * 100
  )
  
  # Identificar o "Elbow" (Cotovelo) matematicamente para anotação
  # (Reutilizando a lógica do seu script ou calculando simples)
  if (show_elbow) {
    coords <- data.frame(x = 1:limit, y = df$Cumulative)
    line_start <- coords[1, ]
    line_end <- coords[limit, ]
    distances <- sapply(1:limit, function(i) {
      pt <- coords[i, ]
      abs((line_end$y - line_start$y) * pt$x - (line_end$x - line_start$x) * pt$y + 
            line_end$x * line_start$y - line_end$y * line_start$x) /
        sqrt((line_end$y - line_start$y)^2 + (line_end$x - line_start$x)^2)
    })
    elbow_pc <- which.max(distances)
    elbow_y <- df$Cumulative[elbow_pc]
  }
  
  # Plot
  p <- ggplot(df, aes(x = PC)) +
    # Barra de Variância Explicada (Azul)
    geom_bar(aes(y = Explained, fill = "Explained Variance"), 
             stat = "identity", width = 0.8, alpha = 0.9) +
    
    # Linha de Variância Acumulada (Vermelha)
    geom_line(aes(y = Cumulative, group = 1, color = "Cumulative Variance"), size = 1.5) +
    geom_point(aes(y = Cumulative, color = "Cumulative Variance"), size = 3) +
    
    # Linhas de referência (80% e 90%)
    geom_hline(yintercept = 80, linetype = "dashed", color = "gray50", size = 0.8) +
    geom_hline(yintercept = 90, linetype = "dashed", color = "gray50", size = 0.8) +
    
    # Escalas e Cores (Mantendo o padrão do seu script)
    scale_y_continuous(limits = c(0, 105), expand = c(0, 0)) +
    scale_fill_manual(name = "", values = c("Explained Variance" = "#1f78b4")) +
    scale_color_manual(name = "", values = c("Cumulative Variance" = "#e31a1c")) +
    
    # Labels
    labs(
      title = title_text,
      x = "Principal Component (PC)",
      y = "Variance (%)"
    ) +
    theme_variance_article
  
  # Adicionar anotação do Cotovelo se solicitado
  if (show_elbow) {
    p <- p + 
      geom_vline(xintercept = elbow_pc, linetype = "dotdash", color = "darkgreen", size = 1) +
      annotate("text", x = elbow_pc + (limit*0.05), y = 50, 
               label = sprintf("Elbow:\nPC%d\n(%.1f%%)", elbow_pc, elbow_y), 
               color = "darkgreen", size = 6, hjust = 0, fontface = "bold")
  }
  
  return(p)
}

# 3. Gerar os dois gráficos para R30
# ------------------------------------------------------------------------------

# Gráfico da Esquerda: Zoom nas primeiras 35 PCs
plot_zoom <- create_variance_plot_custom(
  pca_result_R30, 
  n_limit = 25, 
  title_text = "R30: First 25 Components (Detail)", 
  show_elbow = TRUE
)

# Gráfico da Direita: Todas as PCs (Visão Global)
plot_all <- create_variance_plot_custom(
  pca_result_R30, 
  n_limit = length(pca_result_R30$explained_variance), 
  title_text = "R30: Global Spectrum (All PCs)", 
  show_elbow = FALSE # Elbow faz menos sentido visualmente no full spectrum
)

# 4. Combinar e Salvar
# ------------------------------------------------------------------------------
# Combinar lado a lado com patchwork
combined_r30_variance <- (plot_zoom | plot_all) +
  plot_layout(guides = "collect") & # Legenda única
  theme(legend.position = "bottom")

# Adicionar Título Principal
final_r30_variance_fig <- combined_r30_variance +
  plot_annotation(
    title = "Structural Decomposition of Variance: R30 Condition",
    theme = theme(plot.title = element_text(size = 30, face = "bold", hjust = 0.5))
  )

# Salvar em alta resolução
ggsave("images/pca_variance_R30_focused_25PC.pdf", 
       final_r30_variance_fig, width = 20, height = 10, dpi = 300)

message("Gráfico de variância R30 (25 PCs vs All) salvo com sucesso: images/pca_variance_R30_focused_25PC.pdf")

# ==============================================================================
# 11. TIME-COURSE ANIMATION (GIF) WITH GGANIMATE (ZOOMED)
# ==============================================================================

library(gganimate)
library(gifski) 

message("Generating Time-Course Animation (GIF)... This may take a few minutes.")

# 1. PREPARE TRAJECTORY SEGMENTS
centroids_segments <- centroids_ordered %>%
  mutate(
    xstart = PC1,
    ystart = PC2,
    xend = lead(PC1),
    yend = lead(PC2),
    Day = lead(Day) 
  ) %>%
  filter(!is.na(xend))

# 2. BASE PLOT SETUP
anim_plot <- ggplot(pca_r30_df_ordered, aes(x = PC1, y = PC2, color = Day)) +
  
  # A. Cells
  geom_point(size = 0.8, alpha = 0.6) +
  
  # B. Centroids
  geom_point(data = centroids_ordered, aes(x = PC1, y = PC2), 
             size = 6, shape = 18, color = "black") + 
  
  # C. Progressive Trajectory Path
  geom_segment(data = centroids_segments, 
               aes(x = xstart, y = ystart, xend = xend, yend = yend), 
               arrow = arrow(length = unit(0.5, "cm"), type = "closed"),
               color = "black", linewidth = 1.2, inherit.aes = FALSE) +
  
  # D. Scales and Aesthetics
  scale_color_manual(values = color_palette) +
  
  # --- O ZOOM APLICADO AQUI ---
  # coord_cartesian corta visualmente a tela sem remover dados estatísticos
  coord_cartesian(
    xlim = c(-0.008, max(pca_r30_df_ordered$PC1, na.rm = TRUE)), 
    ylim = c(-0.005, 0.01)
  ) + 
  
  # LEGEND
  guides(color = guide_legend(override.aes = list(size = 10, alpha = 1))) +
  
  # E. Dynamic Labels
  labs(
    title = "EMT Progression Sequence\nCurrent Stage: {closest_state}",
    x = "PC1 (Progression Proxy)", 
    y = "PC2 (Amplitude)"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "right", 
    legend.title = element_blank(),
    legend.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(face = "bold", size = 26, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 20),
    panel.border = element_rect(color = "grey50", fill = NA, linewidth = 1.5)
  ) +
  
  # 3. GGANIMATE COMMANDS
  transition_states(Day, 
                    transition_length = 1, 
                    state_length = 3,
                    wrap = FALSE) + 
  
  shadow_mark(past = TRUE, future = FALSE, alpha = 0.4) 

# 4. RENDER AND SAVE THE GIF
total_days <- length(unique(pca_r30_df_ordered$Day))
seconds_per_day <- 2
fps <- 10

anim_rendered <- animate(
  anim_plot, 
  # loop = TRUE faz o GIF reiniciar infinitamente
  renderer = gifski_renderer(loop = TRUE), 
  width = 1000, height = 700, res = 100, 
  fps = fps, 
  duration = total_days * seconds_per_day,
  # end_pause é medido em frames. 3 segundos * 10 fps = 30 frames
  end_pause = 30 
)

# Save the generated GIF
dir.create("images", showWarnings = FALSE)
anim_save("images/pca_trajectory_animation.gif", animation = anim_rendered)

message("Animation successfully saved (Loop enabled with 3s pause): images/pca_trajectory_animation.gif")