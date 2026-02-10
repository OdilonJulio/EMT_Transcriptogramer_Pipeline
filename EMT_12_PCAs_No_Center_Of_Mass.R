# PCAs_No_Center_Of_Mass.R

# Carregar os dados
load("~/mestrado/pca_result_R30.RData")
load("~/mestrado/t_matrix_R30.RData")

# Carregar bibliotecas necessárias
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(RColorBrewer)
library(dplyr)
library(patchwork)

## 1. Preparação dos dados ----
gene_names <- t_matrix_R30@transcriptogramS2[, 1]
gene_positions <- t_matrix_R30@transcriptogramS2[, 2]
pca_rotation <- pca_result_R30[["pca_result"]][["rotation"]]
rownames(pca_rotation) <- gene_names

n_pcs <- 6
top_percent <- 0.1  # 10% superior e inferior
threshold_abs <- 0.005  # Novo parâmetro para o limite em módulo

## 2. Plot das 10 primeiras PCs (centradas) ----
pc_means <- colMeans(pca_rotation[, 1:n_pcs])

pc_df <- data.frame(
  Gene = rep(gene_names, n_pcs),
  Position = rep(gene_positions, n_pcs),
  PC = rep(paste0("PC", 1:n_pcs), each = length(gene_names)),
  Loading = as.vector(pca_rotation[, 1:n_pcs]) - rep(pc_means, each = length(gene_names))
)

# Padronizar eixo Y
y_range <- range(pc_df$Loading)
y_padding <- diff(y_range) * 0.05
y_limits <- c(y_range[1] - y_padding, y_range[2] + y_padding)

pc_plots <- lapply(1:n_pcs, function(i) {
  current_pc <- paste0("PC", i)
  subset_data <- pc_df[pc_df$PC == current_pc, ]
  
  ggplot(subset_data, aes(x = Position, y = Loading)) +
    geom_line(color = "steelblue", alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    coord_cartesian(ylim = y_limits) +
    labs(title = current_pc, x = "Gene Position", y = "Rotation PCA (centered)") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.margin = unit(c(5, 5, 5, 5), "mm"))
})

combined_pc_plot <- wrap_plots(pc_plots, ncol = 2) + 
  plot_annotation(title = "PCs Centralized")

ggsave("images/PCs_centradas_alta_resolucao.png", combined_pc_plot, 
       width = 16, height = 20, dpi = 600, units = "in")

## 3. Identificar genes acima do limite absoluto ----
get_genes_by_threshold <- function(pc_scores, threshold) {
  peaks <- which(pc_scores >= threshold)
  valleys <- which(pc_scores <= -threshold)
  list(peaks = peaks, valleys = valleys)
}

extreme_genes <- lapply(1:n_pcs, function(i) {
  get_genes_by_threshold(pca_rotation[, i], threshold_abs)
})


extreme_df <- data.frame()
for (i in 1:n_pcs) {
  pc_name <- paste0("PC", i)
  peaks <- extreme_genes[[i]]$peaks
  valleys <- extreme_genes[[i]]$valleys
  
  extreme_df <- rbind(
    extreme_df,
    data.frame(PC = pc_name, Type = "Peak", Gene = gene_names[peaks],
               Position = gene_positions[peaks], Loading = pca_rotation[peaks, i]),
    data.frame(PC = pc_name, Type = "Valley", Gene = gene_names[valleys],
               Position = gene_positions[valleys], Loading = pca_rotation[valleys, i])
  )
}

## 4. Converter ENSP para Gene Names ----
# Se você já tiver um vetor com os nomes dos genes, substitua aqui:
# extreme_df$GeneName <- seus_gene_names[match(extreme_df$Gene, gene_names)]

# Se não tiver, use biomaRt para converter (requer internet):
if (!requireNamespace("biomaRt", quietly = TRUE)) install.packages("biomaRt")
library(biomaRt)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_info <- getBM(
  attributes = c("ensembl_peptide_id", "external_gene_name"),
  filters = "ensembl_peptide_id",
  values = unique(extreme_df$Gene),
  mart = mart
)

extreme_df <- left_join(extreme_df, gene_info, by = c("Gene" = "ensembl_peptide_id"))

## 5. Gráfico dos genes extremos com nomes ----
y_range_extremos <- range(extreme_df$Loading)
y_padding_extremos <- diff(y_range_extremos) * 0.1
y_limits_extremos <- c(y_range_extremos[1] - y_padding_extremos, 
                       y_range_extremos[2] + y_padding_extremos)
# Função para agrupar genes por vizinhança (clusters)
assign_groups <- function(df) {
  if (nrow(df) == 0) return(df)
  df <- df %>% arrange(Position)
  group_id <- c(1, 1 + cumsum(diff(df$Position) > 1))
  df$Cluster <- group_id
  return(df)
}

# Nova função de plotagem com clusters (cores visuais, legenda resumida)
plot_pc_with_clusters <- function(pc_num) {
  pc_name <- paste0("PC", pc_num)
  pc_data <- pc_df[pc_df$PC == pc_name, ]
  extremes <- extreme_df[extreme_df$PC == pc_name, ]
  
  # Separar picos e vales, agrupar em clusters
  peaks <- extremes %>% filter(Type == "Peak") %>% assign_groups()
  valleys <- extremes %>% filter(Type == "Valley") %>% assign_groups()
  
  # Prefixos para distinguir up/down
  if (nrow(peaks) > 0) peaks <- peaks %>% mutate(Cluster = paste0("Up_", Cluster))
  if (nrow(valleys) > 0) valleys <- valleys %>% mutate(Cluster = paste0("Down_", Cluster))
  
  extremes_clustered <- bind_rows(peaks, valleys)
  
  # Contagem de clusters
  n_up <- ifelse(nrow(peaks) > 0, length(unique(peaks$Cluster)), 0)
  n_down <- ifelse(nrow(valleys) > 0, length(unique(valleys$Cluster)), 0)
  
  legend_text <- paste0("Up clusters: ", n_up, " | Down clusters: ", n_down)
  
  # Plot
  pc_plot <- ggplot() +
    geom_line(data = pc_data, aes(x = Position, y = Loading), color = "gray70") +
    geom_hline(yintercept = c(threshold_abs, -threshold_abs), 
               linetype = "dashed", color = "black") +
    geom_point(data = extremes_clustered, 
               aes(x = Position, y = Loading, color = Cluster),
               size = 2, alpha = 0.8) +
    coord_cartesian(ylim = y_limits_extremos) +
    labs(title = pc_name, x = "Gene Position", y = "Rotation PCA") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom",
      plot.margin = unit(c(5, 5, 5, 5), "mm")
    ) +
    guides(color = "none") + # remove legenda dos clusters
    annotate("text", x = Inf, y = Inf, label = legend_text,
             hjust = 1.1, vjust = 1.5, size = 4, fontface = "bold")
  
  return(pc_plot)
}


extreme_plots <- lapply(1:n_pcs, plot_pc_with_clusters)


combined_extreme_plot <- wrap_plots(extreme_plots, ncol = 2) + 
  plot_annotation(title = "Extreme Genes per PC (|loading| > 0.005)")

ggsave("images/extreme_genes_per_PC_loading005.png", combined_extreme_plot, 
       width = 16, height = 20, dpi = 600, units = "in")


# 6 Exportar todos os genes acima do threshold
all_extreme_table <- extreme_df %>%
  mutate(Direction = ifelse(Type == "Peak", "Up", "Down")) %>%
  arrange(PC, Direction, desc(abs(Loading))) %>%
  dplyr::select(PC, Direction, GeneID = Gene, GeneName = external_gene_name, Position, Loading)

write.csv(all_extreme_table, "genes_destaque_threshold_por_PC.csv", row.names = FALSE)

