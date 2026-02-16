# Sequential script unificando a etapa de cálculo por janelas sobrepostas (smoothing manual)
# com a etapa de plotagem/overlay e animação. Configurável no topo.
#
# Requer: pca_result_R30.RData no caminho especificado; objetos opcionais (pca_r30_df,
# color_palette, centroids_r30) serão usados se existirem — caso contrário são gerados fallback simples.

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(gganimate)
library(scales)

# -----------------------------
# 0. Parâmetros (edite aqui)
# -----------------------------
# -----------------------------
num_pcs        <- 6       # PCs beyond PC1 (PC2..PC7 when =6)
num_intervals  <- 50      # number of windows/intervals
overlap_prop   <- 0.3     # overlap (0.5 = 50%)
x_min_limit    <- -0.005   # X axis min
x_max_limit    <- 0.005    # X axis max
y_min_limit    <- -0.01   # Y axis min
y_max_limit    <- 0.01    # Y axis max

# Diretórios de saída

dir.create("frames", showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1. Carregar PCA
# -----------------------------
load("~/mestrado/pca_results_with_metadata.RData")
load("~/mestrado/pca_result_R0.RData")
load("~/mestrado/pca_result_R30.RData")   # cria pca_result_R30
pc_matrix <- pca_result_R30$pca_result$x  # matrix de scores (samples x PCs)

# -----------------------------
# 2. Function: overlapping windows
# -----------------------------
generate_pc1_vs_pcn_data <- function(n, num_intervals, overlap_prop, pc_matrix) {
  data <- data.frame(PC1 = pc_matrix[, 1], PCn = pc_matrix[, n])
  pc1_min <- min(data$PC1, na.rm = TRUE)
  pc1_max <- max(data$PC1, na.rm = TRUE)
  interval_width <- (pc1_max - pc1_min) / num_intervals
  step_size <- interval_width * (1 - overlap_prop)
  starts <- seq(pc1_min, pc1_max - interval_width, by = step_size)
  if (tail(starts, 1) + interval_width < pc1_max) {
    starts <- c(starts, pc1_max - interval_width)
    starts <- unique(starts)
  }
  res <- lapply(starts, function(start_val) {
    end_val <- start_val + interval_width
    window_data <- data[data$PC1 >= start_val & data$PC1 <= end_val, ]
    if (nrow(window_data) > 0) {
      data.frame(
        PC1_mean = mean(window_data$PC1, na.rm = TRUE),
        PCn_mean = mean(window_data$PCn, na.rm = TRUE)
      )
    } else NULL
  })
  res <- do.call(rbind, res)
  if (!is.null(res)) res$PCn_index <- paste0("PC", n)
  res
}
# -----------------------------
# 3. Generate data
# -----------------------------
pc_indices <- 2:(num_pcs + 1)
all_paths_data_list <- lapply(pc_indices,
                              function(n) generate_pc1_vs_pcn_data(n, num_intervals, overlap_prop, pc_matrix))

# CORREÇÃO: Renomeie para 'all_paths_df_full' e remova a filtragem de limites
all_paths_df_full <- bind_rows(all_paths_data_list) %>%
  filter(is.finite(PC1_mean))

# (A filtragem de x_min/x_max foi removida daqui para que o objeto _full seja salvo corretamente)


# -----------------------------
# 4. Aplicar corte X (manter só PC1_mean > x_min_limit)
# -----------------------------
# AGORA ISTO VAI FUNCIONAR, pois 'all_paths_df_full' existe:
all_paths_df <- all_paths_df_full %>%
  filter(is.finite(PC1_mean)) %>%
  filter(PC1_mean > x_min_limit)

# Salvar CSV (completo e filtrado)
write.csv(all_paths_df_full,
          file = paste0("PC_means_by_interval_PC2_to_PC", num_pcs + 1, "_full.csv"),
          row.names = FALSE)
write.csv(all_paths_df,
          file = paste0("PC_means_by_interval_PC2_to_PC", num_pcs + 1, "_filtered_xmin_", x_min_limit, ".csv"),
          row.names = FALSE)

# -----------------------------
# 5. Escalas globais Y (baseado nos dados filtrados)
# -----------------------------
amplitude_by_pc <- all_paths_df %>%
  group_by(PCn_index) %>%
  summarise(
    ymin = min(PCn_mean, na.rm = TRUE),
    ymax = max(PCn_mean, na.rm = TRUE),
    amplitude = ymax - ymin,
    .groups = "drop"
  ) %>%
  arrange(desc(amplitude))

if (nrow(amplitude_by_pc) == 0) {
  stop("Nenhum dado permaneceu após o filtro x_min_limit. Ajuste x_min_limit ou verifique os dados.")
}

ymin <- amplitude_by_pc$ymin[1]
ymax <- amplitude_by_pc$ymax[1]
max_pc <- amplitude_by_pc$PCn_index[1]

# -----------------------------
# 6. Plot facet (linhas calculadas, sem geom_smooth)
# -----------------------------
# -----------------------------
plot_paths <- ggplot(all_paths_df, aes(x = PC1_mean, y = PCn_mean)) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  facet_wrap(~PCn_index, scales = "fixed", ncol = 3) +
  labs(
    title = paste0("PC1 vs PC2..PC", num_pcs + 1,
                   " — ", num_intervals, " windows, overlap ", overlap_prop * 100, "%"),
    x = "PC1 (window mean)",
    y = "PCn (window mean)"
  ) +
  coord_cartesian(xlim = c(x_min_limit, x_max_limit),
                  ylim = c(y_min_limit, y_max_limit)) +
  theme_minimal(base_size = 13)

print(plot_paths)

ggsave("images/PC1_vs_PC2_to_PC_windowed.png", plot_paths,
       width = 16, height = ceiling(num_pcs/3)*3.5, dpi = 300)


# -----------------------------
# 7. Preparar dados de pontos brutos (pca_r30_df) e paleta
#    - Se os objetos já existem no ambiente, usamos; senão criamos fallback.
# -----------------------------
if (!exists("pca_r30_df")) {
  # cria tabela mínima com PCs 1..(num_pcs+1) e Day=NA
  pca_r30_df <- as.data.frame(pc_matrix[, 1:(num_pcs + 1)])
  colnames(pca_r30_df) <- paste0("PC", seq_len(ncol(pca_r30_df)))
  pca_r30_df$Day <- factor(rep("1", nrow(pca_r30_df)))
  message("Aviso: 'pca_r30_df' não existia — criado fallback contendo colunas PC1..")
}

# garantir colunas PC1..PC(num_pcs+1)
needed_cols <- paste0("PC", 1:(num_pcs + 1))
missing_cols <- setdiff(needed_cols, colnames(pca_r30_df))
if (length(missing_cols) > 0) {
  stop("pca_r30_df não contém colunas necessárias: ", paste(missing_cols, collapse = ", "))
}


color_palette <- c(
  "notreated-batch1" = "#984EA3",
  "notreated-batch2" = "#377EB8",
  "TGFbeta1-1day-batch2" = "#4DAF4A",
  "TGFbeta1-2day-batch2" = "#3A5400",
  "TGFbeta1-3day-batch2" = "#A65628",
  "TGFbeta1-4day-batch1" = "#FF7F00",
  "TGFbeta1-8day-batch1" = "#E41A1C"
)

# paleta: se não existir, criar automática com níveis de Day
if (!exists("color_palette")) {
  days_levels <- unique(as.character(pca_r30_df$Day))
  palette_colors <- hue_pal()(max(1, length(days_levels)))
  names(palette_colors) <- days_levels
  color_palette <- palette_colors
  message("Aviso: 'color_palette' não existia — criada paleta automática baseada em Day.")
}

# -----------------------------
# Calcular centroides
# -----------------------------
centroids_r0 <- pca_r0_df %>%
  group_by(Day) %>%
  summarize(across(starts_with("PC"), mean), .groups = "drop")

centroids_r30 <- pca_r30_df %>%
  group_by(Day) %>%
  summarize(across(starts_with("PC"), mean), .groups = "drop")

# centroids_r30: se não existir, calcular médias por Day (fallback)
if (!exists("centroids_r30")) {
  centroids_r30 <- pca_r30_df %>%
    group_by(Day) %>%
    summarise(across(starts_with("PC"), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  message("Aviso: 'centroids_r30' não existia — criado a partir de pca_r30_df (médias por Day).")
}

# -----------------------------
# 8. Criar overlay plots (PC2..PC(num_pcs+1)): pontos brutos + linha suavizada (nossa)
# -----------------------------
overlay_plots <- list()
for (n in pc_indices) {
  pc_name <- paste0("PC", n)
  point_data <- pca_r30_df %>%
    transmute(PC1 = .data$PC1, PCn = .data[[pc_name]], Day = as.factor(Day)) %>%
    filter(PC1 >= x_min_limit & PC1 <= x_max_limit)
  line_data <- all_paths_df %>% filter(PCn_index == pc_name)
  
  p <- ggplot() +
    geom_point(data = point_data,
               aes(x = PC1, y = PCn, color = Day),
               size = 0.3, alpha = 0.5) +  # smaller points
    geom_line(data = line_data,
              aes(x = PC1_mean, y = PCn_mean),
              color = "black", linewidth = 1) +
    # geom_point(data = centroids_r30,
    #            aes(x = PC1, y = .data[[pc_name]]),
    #            size = 1.5, shape = 18, color = "black") +
    scale_color_manual(values = color_palette) +
    labs(title = paste("PC1 vs", pc_name),
         x = "PC1", y = pc_name) +
    coord_fixed(
      ratio = 1,
      xlim = c(x_min_limit, x_max_limit),
      ylim = c(y_min_limit, y_max_limit)
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  overlay_plots[[pc_name]] <- p
}


# -----------------------------
# 9. Combinar grids e salvar (1x6 e 2x3 adaptáveis)
# -----------------------------
combined_grid_1row <- wrap_plots(overlay_plots, nrow = 1, guides = "collect") +
  plot_annotation(title = "PCA R30: Raw points + Smoothed trajectory (PC1 vs PC2-PCn)")

ggsave("images/overlay_pca_R30_PC1_vs_PC2to7_1row.pdf",
       combined_grid_1row, width = 18, height = 6, dpi = 300)

combined_grid_2x3 <- wrap_plots(overlay_plots, nrow = 2, ncol = 3, guides = "collect") +
  plot_annotation(title = "PCA R30: Raw points + Smoothed trajectory")

ggsave("images/overlay_pca_R30_PC1_vs_PC2to7_2x3.pdf",
       combined_grid_2x3, width = 12, height = 8, dpi = 300)

message("Plots saved with fixed axis limits X=[-0.01,0.01], Y=[-0.01,0.01]")

# -----------------------------
# 10. Animação: transição por PCn_index (usa as linhas calculadas)
# -----------------------------
# Dados para animação (usamos all_paths_df)
anim_df <- all_paths_df %>% arrange(PCn_index, PC1_mean)

# Plano de fundo: pontos brutos (filtrados)
bg_points <- pca_r30_df %>%
  transmute(PC1 = .data$PC1, PCn = .data$PC2, Day = as.factor(Day)) %>% # PCn aqui só para estrutura; pontos usados apenas como fundo
  filter(PC1 > x_min_limit)

anim <- ggplot() +
  # fundo com pontos (pontos brutos filtrados)
  # geom_point(data = pca_r30_df %>% filter(PC1 > x_min_limit),
  #            aes(x = PC1, y = .data[[paste0("PC", 2)]]), # somente para preencher fundo; visível em cada estado
  #            size = 0.3, alpha = 0.15) +
  # linha atual (cada estado traz sua própria PCn)
  geom_line(
    data = anim_df,
    aes(x = PC1_mean, y = PCn_mean, group = 1),
    color = "steelblue", size = 1
  ) +
  labs(
    title = 'PC1 vs {closest_state} — overlap {round(overlap_prop*100)}%',
    x = "PC1 (média janela)",
    y = "PCn (média janela)"
  ) +
  coord_fixed(
    ratio = 1,
    xlim = c(x_min_limit, x_max_limit),
    ylim = c(ymin, ymax + 0.005)
  ) +
  theme_minimal(base_size = 14) +
  transition_states(PCn_index, transition_length = 2, state_length = 1) +
  ease_aes('cubic-in-out')

# Renderizar e salvar GIF (gifski)
anim_rendered <- animate(anim,
                         renderer = gifski_renderer(),
                         width = 1000, height = 700,
                         duration = max(6, length(unique(anim_df$PCn_index)) * 1.5),
                         fps = 10)
anim_save("images/PC1_vs_PCn_animation_windowed.gif", animation = anim_rendered)

message("Script finalizado. Arquivos gerados em 'images/':")
message(" - plot de janelas (PNG)")
message(" - overlay PDFs (1x6 e 2x3)")
message(" - CSVs com médias (full e filtered)")
message(" - animação GIF: images/PC1_vs_PCn_animation_windowed.gif")


### ETAPA PEDIDA PELA PROF RITA EM 06/11/2025.

# gerar_tabela_media_pc_por_dia.R
#
# Carrega os resultados do transcriptograma (t_matrix) e roda a análise de PCA
# para calcular o valor médio (centroide) de cada Componente Principal (PC)
# para cada dia de tratamento.
# Salva os resultados em arquivos CSV.

library(dplyr)
library(transcriptogramer) # Necessário para carregar o objeto t_matrix

# --- 1. Carregar dados de entrada (gerados por transcriptogram_analysis.R) ---
message("Carregando dados de entrada (t_matrix R0 e R30)...")

# Verifica se os arquivos de entrada existem
if (!file.exists("t_matrix_R0.RData") || !file.exists("t_matrix_R30.RData")) {
  stop("Arquivos 't_matrix_R0.RData' e 't_matrix_R30.RData' não encontrados. 
       Rode o script 'transcriptogram_analysis.R' primeiro.")
}

load("t_matrix_R0.RData")
load("t_matrix_R30.RData")

# --- 2. Funções copiadas de 'pca_analysis_enhanced.R' ---

# Função para preparar dados (amostras x posições)
prepare_pca_data <- function(transcriptogram_obj) {
  df <- transcriptogram_obj@transcriptogramS2
  df <- df[, !colnames(df) %in% c("Protein", "Position")]
  if (any(is.na(df))) stop("O dataset contém valores NA.")
  return(t(df)) # Retorna Amostras (linhas) x Posições (colunas)
}

# Função para rodar PCA
run_pca_enhanced <- function(pca_input) {
  pca_result <- prcomp(pca_input, center = TRUE, scale. = FALSE)
  # Retorna apenas o objeto pca (simplificado)
  return(list(pca_result = pca_result)) 
}

# Função para adicionar a coluna 'Day' (ajustada para ordenação correta)
add_day_column <- function(pca_df, sample_names) {
  # Mapeia nomes das amostras para os dias
  pca_df$Day <- case_when(
    grepl("notreated-batch1", sample_names) ~ "Dia 0 (B1)",
    grepl("notreated-batch2", sample_names) ~ "Dia 0 (B2)",
    grepl("TGFbeta1-1day-batch2", sample_names) ~ "Dia 1",
    grepl("TGFbeta1-2day-batch2", sample_names) ~ "Dia 2",
    grepl("TGFbeta1-3day-batch2", sample_names) ~ "Dia 3",
    grepl("TGFbeta1-4day-batch1", sample_names) ~ "Dia 4",
    grepl("TGFbeta1-8day-batch1", sample_names) ~ "Dia 8",
    TRUE ~ "Unknown"
  )
  
  # Define a ordem temporal correta (essencial para a tabela final)
  day_order <- c("Dia 0 (B1)", "Dia 0 (B2)", "Dia 1", "Dia 2", "Dia 3", "Dia 4", "Dia 8")
  pca_df$Day <- factor(pca_df$Day, levels = day_order)
  
  return(pca_df)
}

# --- 3. Executar PCA e preparação de dados ---

message("Preparando dados R0...")
df_R0 <- prepare_pca_data(t_matrix_R0)
message("Rodando PCA R0...")
pca_result_R0 <- run_pca_enhanced(df_R0)
# Adiciona metadados de 'Day' aos scores de PCA
pca_r0_df <- add_day_column(as.data.frame(pca_result_R0$pca_result$x), rownames(df_R0))

message("Preparando dados R30...")
df_R30 <- prepare_pca_data(t_matrix_R30)
message("Rodando PCA R30...")
pca_result_R30 <- run_pca_enhanced(df_R30)
# Adiciona metadados de 'Day' aos scores de PCA
pca_r30_df <- add_day_column(as.data.frame(pca_result_R30$pca_result$x), rownames(df_R30))

# --- 4. Calcular Centroides (A tabela que você quer) ---

message("Calculando médias (centroides) por dia...")

# R0
centroids_r0 <- pca_r0_df %>%
  group_by(Day) %>%
  summarize(across(starts_with("PC"), mean, na.rm = TRUE), .groups = "drop") %>%
  arrange(Day) # Ordenar pela ordem temporal definida no fator

# R30
centroids_r30 <- pca_r30_df %>%
  group_by(Day) %>%
  summarize(across(starts_with("PC"), mean, na.rm = TRUE), .groups = "drop") %>%
  arrange(Day) # Ordenar pela ordem temporal definida no fator

# --- 5. Salvar as tabelas em arquivos CSV ---

output_file_r0 <- "media_coeficientes_pca_por_dia_R0.csv"
output_file_r30 <- "media_coeficientes_pca_por_dia_R30.csv"

write.csv(centroids_r0, file = output_file_r0, row.names = FALSE)
write.csv(centroids_r30, file = output_file_r30, row.names = FALSE)

message(paste("\nSucesso! Tabelas salvas em:"))
message(paste("-", output_file_r0))
message(paste("-", output_file_r30))

# Visualizar a tabela R30 no console
print("Tabela R30 (Médias por Dia):")
print(centroids_r30)



# -------------------------------------------------------------------------
# ETAPA EXTRA: Plot Estático PC2 a PC13 com Zoom no Eixo Y
# Solicitado: Intervalo PC2-PC13, Y fixo em [-0.004, 0.004]
# -------------------------------------------------------------------------

message("Gerando plot extra: PC1 vs PC2..PC13 com Y limitado...")

# 1. Gerar dados especificamente para PC2 até PC13
# (Ignoramos o 'num_pcs' do topo para garantir que vá até o 13)
target_pcs <- 2:13

# Verifica se a matriz PCA tem colunas suficientes
if (ncol(pc_matrix) < 13) {
  stop("A matriz PCA (pc_matrix) tem menos de 13 colunas. Não é possível plotar até PC13.")
}

# Reutiliza a função 'generate_pc1_vs_pcn_data' já definida no script
data_13_list <- lapply(target_pcs, function(n) {
  generate_pc1_vs_pcn_data(n, num_intervals, overlap_prop, pc_matrix)
})

# Unir e filtrar (aplica o mesmo filtro de X do restante do script)
df_13 <- bind_rows(data_13_list) %>%
  filter(is.finite(PC1_mean)) %>%
  filter(PC1_mean > x_min_limit)

# Ajustar ordem dos fatores para o Facet não bagunçar (ex: PC10 vir antes de PC2)
df_13$PCn_index <- factor(df_13$PCn_index, levels = paste0("PC", target_pcs))

# 2. Criar o gráfico com as escalas solicitadas
plot_extra <- ggplot(df_13, aes(x = PC1_mean, y = PCn_mean)) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  # Linha tracejada no zero para referência
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", size = 0.3) +
  facet_wrap(~PCn_index, ncol = 4) + # Grid 3x4 ou 4x3 fica bom
  labs(
    title = paste0("PC1 vs PC2 to PC13 - intervals: ", num_intervals),
    x = "PC1 (mean per interval)",
    y = "PCn (mean per interval)"
  ) +
  # Força os limites exatos solicitados
  coord_cartesian(
    xlim = c(x_min_limit, x_max_limit),
    ylim = c(-0.004, 0.004)
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 6)
  )

# 3. Salvar
ggsave("images/PC1_vs_PC2_to_PC13.png", 
       plot_extra, 
       width = 12, height = 8, dpi = 300)

message("Imagem extra salva: images/PC1_vs_PC2_to_PC13.png")
