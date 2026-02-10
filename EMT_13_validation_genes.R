############################################
### 1. CARREGAMENTO DE DADOS E PACOTES ####
############################################

load("~/mestrado/pca_result_R30.RData")
load("~/mestrado/t_matrix_R30.RData")
load("~/mestrado/t_matrix_R0.RData") # Para cálculo do SD

required_packages <- c("transcriptogramer", "ggplot2", "ggrepel", "dplyr", 
                       "RColorBrewer", "patchwork", "biomaRt", 
                       "matrixStats", "gridExtra", "cowplot")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

############################################
### 2. PREPARAÇÃO DOS DADOS DE PCA #########
############################################

gene_names <- t_matrix_R30@transcriptogramS2[, 1]
gene_positions <- t_matrix_R30@transcriptogramS2[, 2]
pca_rotation <- pca_result_R30[["pca_result"]][["rotation"]]
rownames(pca_rotation) <- gene_names

n_pcs <- 6
top_percent <- 0.1

top_extreme_table <- read.csv("genes_destaque_threshold_por_PC.csv", stringsAsFactors = FALSE) %>%
  filter(PC %in% paste0("PC", 1:n_pcs)) %>%
  mutate(Type = ifelse(Direction == "Up", "Peak", "Valley"))

############################################
### 3. CÁLCULO DO DESVIO PADRÃO (SD) ######
############################################

calculate_gene_sd <- function() {
  exp_data <- as.matrix(t_matrix_R0@transcriptogramS2[, -c(1:2)])
  data.frame(
    Gene = t_matrix_R0@transcriptogramS2[, 1],
    Position = t_matrix_R0@transcriptogramS2[, 2],
    SD = rowSds(exp_data, na.rm = TRUE)
  )
}

gene_sd <- calculate_gene_sd()
high_sd_genes <- gene_sd %>%
  arrange(desc(SD)) %>%
  slice_head(prop = top_percent) %>%
  mutate(Type = "High_SD")

############################################
### 4. ANOTAÇÃO DOS GENES ##################
############################################

annotate_genes <- function(gene_ids) {
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  getBM(
    attributes = c("ensembl_peptide_id", "external_gene_name"),
    filters = "ensembl_peptide_id",
    values = gene_ids,
    mart = mart
  )
}

gene_ids <- unique(c(top_extreme_table$GeneID, high_sd_genes$Gene))
gene_annotations <- annotate_genes(gene_ids)

top_extreme_table <- left_join(top_extreme_table, gene_annotations, 
                               by = c("GeneID" = "ensembl_peptide_id"))
high_sd_genes <- left_join(high_sd_genes, gene_annotations,
                           by = c("Gene" = "ensembl_peptide_id"))

write.csv(high_sd_genes %>% dplyr::select(GeneID = Gene, GeneName = external_gene_name, Position, SD),
          "top_high_SD_genes.csv", row.names = FALSE)

############################################
### 5. IDENTIFICAÇÃO DE GENES COMUNS #######
############################################

common_genes <- inner_join(
  top_extreme_table %>% dplyr::select(GeneID, GeneName = external_gene_name, PC, Type, Loading),
  high_sd_genes %>% dplyr::select(Gene, GeneName = external_gene_name),
  by = c("GeneID" = "Gene", "GeneName")
) %>%
  left_join(gene_sd, by = c("GeneID" = "Gene")) %>%
  mutate(Common = TRUE)

write.csv(common_genes %>% dplyr::select(GeneID, GeneName, PC, Loading, SD),
          "genes_in_both_groups.csv", row.names = FALSE)

############################################
### 6. PREPARAÇÃO DO GRÁFICO PRINCIPAL ####
############################################

prepare_pc_data <- function() {
  pc_means <- colMeans(pca_rotation[, 1:n_pcs])
  data.frame(
    Gene = rep(gene_names, n_pcs),
    Position = rep(gene_positions, n_pcs),
    PC = factor(rep(paste0("PC", 1:n_pcs), each = length(gene_names))),
    Loading = as.vector(pca_rotation[, 1:n_pcs]) - rep(pc_means, each = length(gene_names))
  )
}

pc_df <- prepare_pc_data()
pc_colors <- brewer.pal(n_pcs, "Set1")
names(pc_colors) <- paste0("PC", 1:n_pcs)

############################################
### 7. FUNÇÃO PARA TABELAS AUTOMÁTICAS ####
############################################

split_genes_table <- function(df, type_label="Up") {
  if (nrow(df) == 0) {
    return(matrix("-", nrow = 1, ncol = 1))
  }
  
  n_genes <- nrow(df)
  
  # Encontra o melhor layout para a tabela (mais colunas possível, máximo 40)
  max_cols <- min(40, n_genes)
  best_cols <- max_cols
  
  # Tenta encontrar um número de colunas que minimize células vazias
  for (cols in max_cols:1) {
    remainder <- n_genes %% cols
    if (remainder == 0 || (n_genes / cols) <= 1.5) {
      best_cols <- cols
      break
    }
  }
  
  n_col <- best_cols
  n_row <- ceiling(n_genes / n_col)
  
  # Preenche com "-" os espaços vazios
  genes_filled <- c(df$GeneName, rep("-", n_row * n_col - n_genes))
  pos_filled <- c(df$Position, rep("-", n_row * n_col - n_genes))
  
  # Monta matriz
  mat <- matrix(genes_filled, nrow = n_row, ncol = n_col, byrow = TRUE)
  pos_mat <- matrix(pos_filled, nrow = n_row, ncol = n_col, byrow = TRUE)
  
  # Junta nome e posição
  mat_final <- matrix(paste0(mat, " (", pos_mat, ")"), nrow = n_row, ncol = n_col)
  colnames(mat_final) <- paste0(type_label, "_", 1:n_col)
  
  return(mat_final)
}

############################################
### 8. CRIAÇÃO DO GRÁFICO + TABELAS   #####
############################################

# Helper: constrói uma "tabela" em grid usando ggplot (melhor controle de colunas)
.make_table_plot <- function(df, title, group_colors,
                             target_rows = 6, max_cols = 12, min_cols = 3) {
  if (nrow(df) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = paste0(title, " (0 clusters, 0 genes)")) +
        theme_void() +
        theme(plot.margin = margin(10, 15, 10, 15))
    )
  }
  
  df <- df %>%
    dplyr::arrange(Group, Position) %>%
    mutate(Cell = paste0(GeneName, " (", Position, ")"))
  
  n <- nrow(df)
  # Escolhe nº de colunas para ~target_rows linhas, com limites
  n_col <- max(min_cols, min(max_cols, ceiling(n / target_rows)))
  n_row <- ceiling(n / n_col)
  
  # Posiciona células coluna-a-coluna (equilibra altura)
  Row <- rep(1:n_row, times = n_col)[1:n]
  Col <- rep(1:n_col, each = n_row)[1:n]
  
  cells <- tibble::tibble(Row = Row, Col = Col,
                          Cell = df$Cell, Group = df$Group)
  
  # Preenche células vazias com "-"
  n_fill <- n_row * n_col - n
  if (n_fill > 0) {
    Row_fill <- rep(1:n_row, times = n_col)[(n + 1):(n_row * n_col)]
    Col_fill <- rep(1:n_col, each = n_row)[(n + 1):(n_row * n_col)]
    cells <- dplyr::bind_rows(
      cells,
      tibble::tibble(Row = Row_fill, Col = Col_fill,
                     Cell = "-", Group = NA_character_)
    )
  }
  
  p <- ggplot(cells, aes(x = Col, y = -Row, fill = Group, label = Cell)) +
    geom_tile(color = "white", linewidth = 0.2) +
    geom_text(size = 2.4, lineheight = 0.9) +
    scale_x_continuous(expand = c(0, 0), breaks = NULL) +
    scale_y_continuous(expand = c(0, 0), breaks = NULL) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
          plot.margin = margin(10, 15, 10, 15),
          legend.position = "none") +
    labs(title = sprintf("%s (%d clusters, %d genes)",
                         title,
                         dplyr::n_distinct(na.omit(df$Group)), n))
  
  if (length(group_colors)) {
    p <- p + scale_fill_manual(values = group_colors, na.value = "grey95")
  } else {
    p <- p + scale_fill_manual(values = c("grey80"), na.value = "grey95")
  }
  p
}

create_main_plot_with_tables <- function(show_labels = FALSE, cutoff = 0.005) {
  # 1) Filtra extremos pelo cutoff
  filtered_extremes <- top_extreme_table %>%
    dplyr::filter(abs(Loading) > cutoff)
  
  # 2) Interseções com high SD
  common_genes_cut <- dplyr::inner_join(
    filtered_extremes %>%
      dplyr::select(GeneID, GeneName = external_gene_name, PC, Type, Loading, Position),
    high_sd_genes %>%
      dplyr::select(Gene, GeneName = external_gene_name, Position),
    by = c("GeneID" = "Gene", "GeneName", "Position")
  ) %>%
    dplyr::left_join(gene_sd %>% dplyr::select(Gene, SD),
                     by = c("GeneID" = "Gene")) %>%
    dplyr::mutate(Common = TRUE)
  
  # 3) Agrupamento por vizinhança (clusters) — separado para up/down
  assign_groups <- function(df) {
    if (nrow(df) == 0) return(df)
    df <- df %>% dplyr::arrange(.data$Position)
    group_id <- c(1, 1 + cumsum(diff(df$Position) > 2))
    df$Group <- paste0("G", group_id)
    df
  }
  
  down_table <- common_genes_cut %>%
    dplyr::filter(Loading < 0) %>%
    dplyr::distinct(GeneID, GeneName, Position, .keep_all = TRUE) %>%
    assign_groups()
  
  up_table <- common_genes_cut %>%
    dplyr::filter(Loading > 0) %>%
    dplyr::distinct(GeneID, GeneName, Position, .keep_all = TRUE) %>%
    assign_groups()
  
  # 4) Paleta de cores por cluster (compartilhada entre up/down)
  all_groups <- unique(c(down_table$Group, up_table$Group))
  group_colors <- if (length(all_groups)) {
    stats::setNames(colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(all_groups)),
                    all_groups)
  } else {
    character(0)
  }
  
  # 5) Contagem de clusters para a legenda/descrição
  n_down_clusters <- dplyr::n_distinct(down_table$Group)
  n_up_clusters   <- dplyr::n_distinct(up_table$Group)
  
  # 6) Dados p/ plotar interseções com as MESMAS CORES dos clusters
  intersections <- dplyr::bind_rows(
    down_table %>% dplyr::mutate(Direction = "Down"),
    up_table %>% dplyr::mutate(Direction = "Up")
  )
  
  # 7) Gráfico principal
  y_range  <- range(pc_df$Loading)
  y_limits <- y_range + c(-1, 1) * diff(y_range) * 0.05
  
  p_main <- ggplot() +
    # PC lines (translúcidas)
    geom_line(data = pc_df, aes(x = Position, y = Loading, color = PC),
              alpha = 0.3, linewidth = 0.5) +
    # linhas de referência
    geom_hline(yintercept = c(-cutoff, cutoff), linetype = "dashed", color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    # extremos não-interseção (translúcidos)
    geom_point(data = filtered_extremes,
               aes(x = Position, y = Loading),
               shape = ifelse(filtered_extremes$Type == "Peak", 24, 25),
               color = "gray60", fill = "gray85", alpha = 0.45, size = 1.3) +
    # interseções por cluster (cores vivas e consistentes com as tabelas)
    geom_point(data = intersections %>% dplyr::filter(Direction == "Down"),
               aes(x = Position, y = Loading, fill = Group),
               shape = 25, color = "black", size = 2.2) +
    geom_point(data = intersections %>% dplyr::filter(Direction == "Up"),
               aes(x = Position, y = Loading, fill = Group),
               shape = 24, color = "black", size = 2.2) +
    # escalas
    scale_color_manual(values = pc_colors) +
    { if (length(group_colors)) scale_fill_manual(values = group_colors, guide = "none")
      else scale_fill_manual(values = c("grey80"), guide = "none") } +
    scale_y_continuous(limits = y_limits,
                       sec.axis = sec_axis(~ . / max(y_limits) * max(high_sd_genes$SD),
                                           name = "Standard Deviation")) +
    labs(
      title = "Principal Components 1–6 with Intersection Clusters",
      subtitle = paste0(
        "Cutoff: |loading| > ", cutoff,
        "   —   Intersection clusters: Down = ", n_down_clusters,
        " | Up = ", n_up_clusters,
        "\nGray triangles: non-intersection extremes; colored triangles: intersection clusters"
      ),
      x = "Gene Position", y = "PCA Rotation (centered)",
      color = "PC"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom",
          plot.margin = margin(20, 20, 12, 20))
  
  if (show_labels) {
    p_main <- p_main +
      ggrepel::geom_text_repel(
        data = intersections,
        aes(x = Position, y = Loading, label = GeneName),
        size = 2.3, fontface = "bold", color = "black",
        max.overlaps = Inf, box.padding = 0.6
      )
  }
  
  # 8) Tabelas (com colunas calculadas, margens e cores por cluster)
  down_plot <- .make_table_plot(down_table, "Down-regulated", group_colors,
                                target_rows = 6, max_cols = 12, min_cols = 3)
  up_plot   <- .make_table_plot(up_table,   "Up-regulated",   group_colors,
                                target_rows = 6, max_cols = 12, min_cols = 3)
  
  # 9) Alturas relativas (cresce um pouco com o nº de linhas estimado)
  estimate_rows <- function(n) {
    if (n == 0) return(1)
    n_col <- max(3, min(12, ceiling(n / 6)))
    ceiling(n / n_col)
  }
  h_down <- 1 + 0.15 * estimate_rows(nrow(down_table))
  h_up   <- 1 + 0.15 * estimate_rows(nrow(up_table))
  
  # primeiro empilha as tabelas
  tables_plot <- cowplot::plot_grid(
    down_plot, up_plot,
    ncol = 1, rel_heights = c(h_down, h_up), align = "v"
  )
  
  # depois coloca lado a lado: gráfico à esquerda, tabelas à direita
  final_plot <- cowplot::plot_grid(
    p_main, tables_plot,
    nrow = 1,               # um ao lado do outro
    rel_widths = c(2, 1)    # mais espaço para o gráfico (2x) do que para as tabelas
  )
  
  
  return(final_plot)
}

# Exemplo de uso e salvamento
final_plot <- create_main_plot_with_tables(show_labels = FALSE, cutoff = 0.001)
ggsave("images/PCs_linhas_intersections_clusters.png", final_plot,
       width = 32, height = 18, dpi = 600, units = "in")

#### PER PC

create_plots_per_pc <- function(show_labels = FALSE, cutoff = 0.005) {
  
  # Filtra extremos acima do cutoff
  filtered_extremes <- top_extreme_table %>%
    dplyr::filter(abs(Loading) > cutoff)
  
  # lista PCs
  pcs <- paste0("PC", 1:n_pcs)
  
  plots_list <- lapply(pcs, function(pc) {
    # dados apenas da PC específica
    pc_data <- pc_df %>% dplyr::filter(PC == pc)
    extremes_pc <- filtered_extremes %>% dplyr::filter(PC == pc)
    
    # interseções high SD só para esta PC
    common_genes_cut <- dplyr::inner_join(
      extremes_pc %>%
        dplyr::select(GeneID, GeneName = external_gene_name, PC, Type, Loading, Position),
      high_sd_genes %>%
        dplyr::select(Gene, GeneName = external_gene_name, Position),
      by = c("GeneID" = "Gene", "GeneName", "Position")
    ) %>%
      dplyr::left_join(gene_sd %>% dplyr::select(Gene, SD),
                       by = c("GeneID" = "Gene")) %>%
      dplyr::mutate(Common = TRUE)
    
    assign_groups <- function(df) {
      if (nrow(df) == 0) return(df)
      df <- df %>% dplyr::arrange(Position)
      group_id <- c(1, 1 + cumsum(diff(df$Position) > 1))
      df$Group <- paste0("G", group_id)
      df
    }
    
    down_table <- common_genes_cut %>%
      dplyr::filter(Loading < 0) %>%
      dplyr::distinct(GeneID, GeneName, Position, .keep_all = TRUE) %>%
      assign_groups()
    
    up_table <- common_genes_cut %>%
      dplyr::filter(Loading > 0) %>%
      dplyr::distinct(GeneID, GeneName, Position, .keep_all = TRUE) %>%
      assign_groups()
    
    # cores por cluster desta PC
    all_groups <- unique(c(down_table$Group, up_table$Group))
    group_colors <- if (length(all_groups)) {
      stats::setNames(colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(all_groups)),
                      all_groups)
    } else {
      character(0)
    }
    
    intersections <- dplyr::bind_rows(
      down_table %>% dplyr::mutate(Direction = "Down"),
      up_table %>% dplyr::mutate(Direction = "Up")
    )
    
    y_range  <- range(pc_data$Loading)
    y_limits <- y_range + c(-1, 1) * diff(y_range) * 0.05
    
    # gráfico principal desta PC
    p_main <- ggplot() +
      geom_line(data = pc_data, aes(x = Position, y = Loading),
                color = pc_colors[pc], alpha = 0.6, linewidth = 0.6) +
      geom_hline(yintercept = c(-cutoff, cutoff), linetype = "dashed", color = "black") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_point(data = extremes_pc,
                 aes(x = Position, y = Loading),
                 shape = ifelse(extremes_pc$Type == "Peak", 24, 25),
                 color = "gray60", fill = "gray85", alpha = 0.45, size = 1.3) +
      geom_point(data = intersections %>% dplyr::filter(Direction == "Down"),
                 aes(x = Position, y = Loading, fill = Group),
                 shape = 25, color = "black", size = 2.2) +
      geom_point(data = intersections %>% dplyr::filter(Direction == "Up"),
                 aes(x = Position, y = Loading, fill = Group),
                 shape = 24, color = "black", size = 2.2) +
      { if (length(group_colors)) scale_fill_manual(values = group_colors, guide = "none")
        else scale_fill_manual(values = c("grey80"), guide = "none") } +
      scale_y_continuous(limits = y_limits) +
      labs(
        title = paste0(pc, " – Extreme Loadings and Intersection Clusters"),
        x = "Gene Position", y = "PCA Rotation (centered)"
      ) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "none",
            plot.margin = margin(20, 20, 12, 20))
    
    if (show_labels) {
      p_main <- p_main +
        ggrepel::geom_text_repel(
          data = intersections,
          aes(x = Position, y = Loading, label = GeneName),
          size = 2.3, fontface = "bold", color = "black",
          max.overlaps = Inf, box.padding = 0.6
        )
    }
    
    down_plot <- .make_table_plot(down_table, "Down-regulated", group_colors,
                                  target_rows = 6, max_cols = 12, min_cols = 3)
    up_plot   <- .make_table_plot(up_table,   "Up-regulated",   group_colors,
                                  target_rows = 6, max_cols = 12, min_cols = 3)
    
    # altura relativa das tabelas
    estimate_rows <- function(n) {
      if (n == 0) return(1)
      n_col <- max(3, min(12, ceiling(n / 6)))
      ceiling(n / n_col)
    }
    h_down <- 1 + 0.15 * estimate_rows(nrow(down_table))
    h_up   <- 1 + 0.15 * estimate_rows(nrow(up_table))
    
    cowplot::plot_grid(
      p_main, down_plot, up_plot,
      ncol = 1, rel_heights = c(6, h_down, h_up), align = "v", axis = "lr"
    )
  })
  
  plots_list
}

# Uso:
plots_per_pc <- create_plots_per_pc(show_labels = FALSE, cutoff = 0.001)

# Salvar cada gráfico em arquivo separado
for (i in seq_along(plots_per_pc)) {
  pc_name <- paste0("PC", i)
  ggsave(paste0("images/", pc_name, "_linhas_intersections_clusters.png"),
         plots_per_pc[[i]], width = 16, height = 16, dpi = 600, units = "in")
}


############################################
### 9. EXPORTAÇÃO DAS TABELAS P/ CSV ######
############################################

export_tables_to_csv <- function(cutoff = 0.005, file = "genes_clusters.csv") {
  filtered_extremes <- top_extreme_table %>%
    dplyr::filter(abs(Loading) > cutoff)
  
  common_genes_cut <- dplyr::inner_join(
    filtered_extremes %>%
      dplyr::select(GeneID, GeneName = external_gene_name, PC, Type, Loading, Position),
    high_sd_genes %>%
      dplyr::select(Gene, GeneName = external_gene_name, Position),
    by = c("GeneID" = "Gene", "GeneName", "Position")
  ) %>%
    dplyr::left_join(gene_sd %>% dplyr::select(Gene, SD),
                     by = c("GeneID" = "Gene")) %>%
    dplyr::mutate(Common = TRUE)
  
  assign_groups <- function(df) {
    if (nrow(df) == 0) return(df)
    df <- df %>% dplyr::arrange(.data$Position)
    group_id <- c(1, 1 + cumsum(diff(df$Position) > 1))
    df$Group <- paste0("G", group_id)
    df
  }
  
  down_table <- common_genes_cut %>%
    dplyr::filter(Loading < 0) %>%
    dplyr::distinct(GeneID, GeneName, Position, .keep_all = TRUE) %>%
    assign_groups() %>%
    mutate(Direction = "Down")
  
  up_table <- common_genes_cut %>%
    dplyr::filter(Loading > 0) %>%
    dplyr::distinct(GeneID, GeneName, Position, .keep_all = TRUE) %>%
    assign_groups() %>%
    mutate(Direction = "Up")
  
  all_clusters <- bind_rows(down_table, up_table) %>%
    arrange(Direction, Group, Position) %>%
    dplyr::select(Direction, Group, GeneID, GeneName, Position, PC, Loading, SD)
  
  write.csv(all_clusters, file, row.names = FALSE)
  
  message("Arquivo salvo em: ", file)
  return(all_clusters)
}

# Exemplo de uso:
clusters_df <- export_tables_to_csv(cutoff = 0.005,
                                    file = "genes_clusters.csv")


