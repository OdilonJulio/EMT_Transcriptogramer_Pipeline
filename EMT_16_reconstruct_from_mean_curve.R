# ==============================================================================
# SCRIPT COMPLETO E FINAL: EIXO X LIMPO + ORDEM CORRETA + ALTA RESOLUÇÃO
# ==============================================================================

# --- 1. Requisitos e Carga de Dados ---
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(gifski)
library(magick)

# Carregar objetos (Ajuste os caminhos se necessário)
load("~/mestrado/pca_result_R30.RData")
load("~/mestrado/t_matrix_R30.RData")

# Carregar dados
mean_df <- read_csv("PC_means_by_interval_PC2_to_PC7_full.csv")

# Criar ID do intervalo
mean_df <- mean_df %>%
  mutate(interval = paste0("frame_", round(PC1_mean, 5))) 

# --- Configuração ---
rotation_matrix <- pca_result_R30$pca_result$rotation
gene_names <- t_matrix_R30@transcriptogramS2[, 1]

pcs_presentes <- unique(mean_df$PCn_index)
pcs_presentes <- pcs_presentes[order(as.integer(gsub("PC", "", pcs_presentes)))]

# --- Construir matriz larga ---
wide_matrix <- mean_df %>%
  dplyr::select(interval, PC1_mean, PCn_index, PCn_mean) %>%
  tidyr::pivot_wider(names_from = PCn_index, values_from = PCn_mean) %>%
  arrange(PC1_mean)

wide_matrix[is.na(wide_matrix)] <- 0

# --- Reconstrução ---
reconstructed_list <- lapply(seq_len(nrow(wide_matrix)), function(i) {
  row <- wide_matrix[i, ]
  pc1_coeff <- row$PC1_mean
  reconstructed_vector <- pc1_coeff * rotation_matrix[, 1]
  for (pc_index in pcs_presentes) {
    pc_number <- as.integer(gsub("PC", "", pc_index))
    pc_coeff <- row[[pc_index]]
    reconstructed_vector <- reconstructed_vector + (pc_coeff * rotation_matrix[, pc_number])
  }
  data.frame(
    Gene = gene_names,
    Expression = as.numeric(reconstructed_vector),
    PC1_value = pc1_coeff,
    IntervalID = row$interval
  )
})

reconstructed_from_means_df <- bind_rows(reconstructed_list)
write.csv(reconstructed_from_means_df, "transcriptograms_reconstructed_by_PC1_weighted.csv", row.names = FALSE)


# ==============================================================================
# --- 2. GERAÇÃO DOS FRAMES (SEM VALORES NO EIXO X) ---
# ==============================================================================

cat("Gerando frames limpos...\n")
dir.create("frames_pc1_real", showWarnings = FALSE, recursive = TRUE)

global_y_min <- min(reconstructed_from_means_df$Expression)
global_y_max <- max(reconstructed_from_means_df$Expression)

interval_df <- reconstructed_from_means_df %>%
  distinct(IntervalID, PC1_value) %>%
  arrange(PC1_value)

n_frames <- nrow(interval_df)
png_files_real <- character(n_frames)

# --- LAYOUT ---
y_range <- global_y_max - global_y_min
y_top_limit <- global_y_max + (0.45 * y_range) 
bar_ymin <- global_y_max + (0.38 * y_range)
bar_ymax <- global_y_max + (0.41 * y_range)
bar_text_y <- bar_ymin + (bar_ymax - bar_ymin)/2 
title_y <- global_y_max + (0.3 * y_range)

for (i in seq_len(n_frames)) {
  current_interval <- interval_df$IntervalID[i]
  current_pc1_val <- interval_df$PC1_value[i]
  progress <- i / n_frames
  percent_txt <- paste0(round(progress*100), "%")
  
  plot_df_real <- reconstructed_from_means_df %>%
    filter(IntervalID == current_interval) %>%
    mutate(ExpressionPlot = Expression)
  
  total_genes <- nrow(plot_df_real)
  center_x <- total_genes / 2
  
  p_real <- ggplot(plot_df_real, aes(x = seq_along(ExpressionPlot), y = ExpressionPlot)) +
    geom_line(color = "black", linewidth = 0.8) + 
    
    # Cabeçalho e Barra
    annotate("rect", xmin = 0, xmax = total_genes, ymin = bar_ymin, ymax = bar_ymax, color = "grey80", fill = "white") +
    annotate("rect", xmin = 0, xmax = total_genes * progress, ymin = bar_ymin, ymax = bar_ymax, fill = "black") +
    annotate("text", x = center_x, y = bar_text_y, label = percent_txt, size = 3, fontface = "bold", color = ifelse(progress > 0.5, "white", "black")) + 
    annotate("text", x = center_x, y = title_y, label = paste0("             Transcriptogram Real Values (PC1 = ", signif(current_pc1_val, 4), ")"), size = 5, fontface = "bold") +
    
    # Eixos (Note x = "" para remover o título "Gene Position")
    labs(x = "", y = "Expression") + 
    
    theme_minimal(base_size = 14) +
    
    # === AQUI ESTÁ A ALTERAÇÃO PARA REMOVER OS VALORES DO EIXO X ===
    theme(
      plot.title = element_blank(),
      axis.text.x = element_blank(),  # Remove os números (0, 2000, 4000...)
      axis.ticks.x = element_blank(), # Remove os tracinhos
      axis.title.x = element_blank()  # Garante que não tenha título no X
    ) +
    # ===============================================================
  
  scale_y_continuous(
    position = "right", 
    limits = c(global_y_min - 0.1*y_range, y_top_limit)
  )
  
  filename_real <- file.path("frames_pc1_real", paste0(current_interval, ".png"))
  ggsave(filename_real, plot = p_real, width = 10, height = 4)
  png_files_real[i] <- filename_real
  
  if (i %% 20 == 0) cat("\rFrame", i, "de", n_frames, "concluído.")
}
cat("\nTodos os frames gerados.\n")


# ==============================================================================
# --- 3. GERAÇÃO DAS IMAGENS DE ALTA RESOLUÇÃO ---
# ==============================================================================

# --- CONFIGURAÇÕES ---
input_frames_dir <- "frames_pc1_real"
bg_file          <- "Biologic.jpg"       
output_dir       <- "Figuras/HighRes_Frames" 

unlink(output_dir, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

TARGET_WIDTH <- 2500 

# Tamanho e Posição
fator <- 1.72
base_width  <- 1200
base_height <- 480
PLOT_WIDTH   <- base_width * fator
PLOT_HEIGHT  <- base_height * fator

POS_X        <- 340    
POS_Y        <- 120   

cat("\nPreparando background...\n")
if (!file.exists(bg_file)) stop("Erro: Biologic.jpg não encontrado.")

bg_base <- tryCatch({
  img <- image_read(bg_file)
  image_scale(img, geometry_size_pixels(width = TARGET_WIDTH))
}, error = function(e) stop("Erro ao ler Biologic.jpg."))
gc() 

# --- SELEÇÃO DE FRAMES (ORDEM CORRETA) ---
ordered_ids <- interval_df$IntervalID
indices <- round(seq(1, length(ordered_ids), length.out = 8))
selected_ids <- ordered_ids[indices]
selected_frames <- file.path(input_frames_dir, paste0(selected_ids, ".png"))

cat("Gerando 8 imagens finais de alta resolução...\n")

for (i in seq_along(selected_frames)) {
  cat("\rProcessando Imagem", i, "de 8...")
  
  if (!file.exists(selected_frames[i])) stop(paste("Arquivo faltando:", selected_frames[i]))
  
  # Processamento
  plot_img <- image_read(selected_frames[i])
  plot_img <- image_resize(plot_img, geometry_size_pixels(width = PLOT_WIDTH, height = PLOT_HEIGHT), filter = "Lanczos")
  
  offset_str <- paste0("+", POS_X, "+", POS_Y)
  final_comp <- image_composite(bg_base, plot_img, operator = "Over", offset = offset_str, gravity = "NorthWest")
  
  letra <- LETTERS[i]
  final_comp <- image_annotate(final_comp, letra, size = 100, color = "black", boxcolor = "white", location = "+50+50")
  
  filename <- file.path(output_dir, paste0("Fig_Part_", i, ".png"))
  image_write(final_comp, path = filename, format = "png", density = 300)
  
  rm(plot_img, final_comp)
  gc()
}

cat("\n\nSucesso! Imagens salvas em 'Figuras/HighRes_Frames'.\n")