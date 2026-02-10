library(transcriptogramer)
library(ggplot2)
library(gridExtra)

## Funções Compartilhadas -----

# Função para plotagem (MANTIDA EXATAMENTE COMO ESTAVA)
create_comparison_plot <- function(original, reconstructed, title, y_lim = NULL) {
  df <- data.frame(
    Position = original$Position,
    Original = 0,
    Reconstructed = reconstructed - original$Expression
  )
  
  p <- ggplot(df, aes(x = Position)) +
    geom_line(aes(y = Original, color = "Original"), linewidth = 0.2, alpha = 0.8) +
    geom_line(aes(y = Reconstructed, color = "Reconstructed"), linewidth = 0.15, alpha = 0.6) +
    scale_color_manual(values = c("Original" = "black", "Reconstructed" = "red")) +
    labs(title = title, x = "Gene", y = "Δ Expression", color = "") +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(size = 9, face = "bold"),
      axis.text.x = element_text(size = 6, angle = 0),
      axis.text.y = element_text(size = 7),
      legend.position = "top",
      panel.grid.major = element_line(linewidth = 0.05),
      panel.grid.minor = element_blank()
    )
  
  if (!is.null(y_lim)) {
    p <- p + ylim(-y_lim, y_lim)
  }
  return(p)
}

## Funções R0 (MANTIDAS AS FUNCIONALIDADES, APENAS OTIMIZADAS) -----

load_required_data_R0 <- function() {
  if (!file.exists("t_matrix_R0.RData") || !file.exists("pca_result_R0.RData")) {
    stop("Arquivos R0 faltando: t_matrix_R0.RData e/ou pca_result_R0.RData")
  }
  load("t_matrix_R0.RData", envir = .GlobalEnv)
  load("pca_result_R0.RData", envir = .GlobalEnv)
  message("Dados R0 carregados com sucesso")
}

get_cell_data_R0 <- function(cell_name) {
  cell_columns <- setdiff(colnames(t_matrix_R0@transcriptogramS2), c("Protein", "Position"))
  
  if (!cell_name %in% cell_columns) {
    stop("Célula não encontrada em R0. Exemplos: ", paste(head(cell_columns), collapse = ", "))
  }
  
  data.frame(
    Position = t_matrix_R0@transcriptogramS2$Position,
    Expression = t_matrix_R0@transcriptogramS2[[cell_name]],
    Gene = rownames(t_matrix_R0@transcriptogramS2)
  )
}

reconstruct_expression_R0 <- function(cell_name, n_pcs) {
  if (!cell_name %in% rownames(pca_result_R0$pca_result$x)) {
    stop("Célula não encontrada no PCA de R0")
  }
  
  rotation <- pca_result_R0$pca_result$rotation[, 1:n_pcs, drop = FALSE]
  scores <- pca_result_R0$pca_result$x[cell_name, 1:n_pcs, drop = FALSE]
  reconstructed <- scores %*% t(rotation)
  
  # Centraliza a reconstrução para média zero
  # reconstructed <- reconstructed - mean(reconstructed)
  
  if (!is.null(pca_result_R0$pca_result$center)) {
    reconstructed <- reconstructed + pca_result_R0$pca_result$center
  }
  
  setNames(as.numeric(reconstructed), colnames(reconstructed))
}

analyze_cell_R0 <- function(cell_name, max_pcs = 5) {
  load_required_data_R0()
  original_data <- get_cell_data_R0(cell_name)
  
  results <- list()
  diffs <- list()
  
  # Primeiro, calcular todas as diferenças para achar o maior valor em módulo
  for (n in 1:max_pcs) {
    reconstructed <- reconstruct_expression_R0(cell_name, n)
    diffs[[n]] <- reconstructed - original_data$Expression
  }
  y_lim <- max(abs(unlist(diffs)), na.rm = TRUE)
  
  plots <- list()
  for (n in 1:max_pcs) {
    reconstructed <- reconstruct_expression_R0(cell_name, n)
    
    plots[[n]] <- create_comparison_plot(
      original_data,
      reconstructed,
      paste("Using", n, ifelse(n == 1, "PC", "PCs")),
      y_lim = y_lim
    )
    
    results[[paste0("PC", n)]] <- reconstructed
  }
  
  combined_plot <- grid.arrange(
    grobs = plots, 
    ncol = 2, 
    top = paste("R0 - Analysis for cell:", cell_name)
  )
  
  output_dir <- "images"
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  filename <- paste0(output_dir, "/R0_",
                     gsub("[^[:alnum:]]", "_", cell_name), "_", max_pcs, "PCs.png")
  
  ggsave(
    filename = filename,
    plot = combined_plot,
    width = 12, 
    height = 6,
    dpi = 600,
    bg = "white"
  )
  
  print(combined_plot)
  
  list(
    original = original_data,
    reconstructions = results,
    combined_plot = combined_plot
  )
}


## Funções R30 (MANTIDAS AS FUNCIONALIDADES, APENAS OTIMIZADAS) -----

load_required_data_R30 <- function() {
  if (!file.exists("t_matrix_R30.RData") || !file.exists("pca_result_R30.RData")) {
    stop("Arquivos R30 faltando: t_matrix_R30.RData e/ou pca_result_R30.RData")
  }
  load("t_matrix_R30.RData", envir = .GlobalEnv)
  load("pca_result_R30.RData", envir = .GlobalEnv)
  message("Dados R30 carregados com sucesso")
}

get_cell_data_R30 <- function(cell_name) {
  cell_columns <- setdiff(colnames(t_matrix_R30@transcriptogramS2), c("Protein", "Position"))
  
  if (!cell_name %in% cell_columns) {
    stop("Célula não encontrada em R30. Exemplos: ", paste(head(cell_columns), collapse = ", "))
  }
  
  data.frame(
    Position = t_matrix_R30@transcriptogramS2$Position,
    Expression = t_matrix_R30@transcriptogramS2[[cell_name]],
    Gene = rownames(t_matrix_R30@transcriptogramS2)
  )
}

reconstruct_expression_R30 <- function(cell_name, n_pcs) {
  if (!cell_name %in% rownames(pca_result_R30$pca_result$x)) {
    stop("Célula não encontrada no PCA de R30")
  }
  
  rotation <- pca_result_R30$pca_result$rotation[, 1:n_pcs, drop = FALSE]
  scores <- pca_result_R30$pca_result$x[cell_name, 1:n_pcs, drop = FALSE]
  reconstructed <- scores %*% t(rotation)
  
  # Centraliza a reconstrução para média zero
  # reconstructed <- reconstructed - mean(reconstructed)
  
  if (!is.null(pca_result_R30$pca_result$center)) {
    reconstructed <- reconstructed + pca_result_R30$pca_result$center
  }
  
  setNames(as.numeric(reconstructed), colnames(reconstructed))
}

analyze_cell_R30 <- function(cell_name, max_pcs = 5) {
  load_required_data_R30()
  original_data <- get_cell_data_R30(cell_name)
  
  results <- list()
  diffs <- list()
  
  # Primeiro, calcular todas as diferenças para achar o maior valor em módulo
  for (n in 1:max_pcs) {
    reconstructed <- reconstruct_expression_R30(cell_name, n)
    diffs[[n]] <- reconstructed - original_data$Expression
  }
  y_lim <- max(abs(unlist(diffs)), na.rm = TRUE)
  
  plots <- list()
  for (n in 1:max_pcs) {
    reconstructed <- reconstruct_expression_R30(cell_name, n)
    
    plots[[n]] <- create_comparison_plot(
      original_data,
      reconstructed,
      paste("Using", n, ifelse(n == 1, "PC", "PCs")),
      y_lim = y_lim
    )
    
    results[[paste0("PC", n)]] <- reconstructed
  }
  
  combined_plot <- grid.arrange(
    grobs = plots, 
    ncol = 2, 
    top = paste("R30 - Analysis for cell:", cell_name)
  )
  
  output_dir <- "images"
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  filename <- paste0(output_dir, "/R30_",
                     gsub("[^[:alnum:]]", "_", cell_name), "_", max_pcs, "PCs.png")
  
  ggsave(
    filename = filename,
    plot = combined_plot,
    width = 12, 
    height = 6,
    dpi = 600,
    bg = "white"
  )
  
  print(combined_plot)
  
  list(
    original = original_data,
    reconstructions = results,
    combined_plot = combined_plot
  )
}


## Execução das análises -----

# 1. Análise para células R0
cells_R0 <- c(
  "TGFbeta1-8day-batch1_AACACGTCAGATCGGA-1",
  "TGFbeta1-4day-batch1_CACACAACAGGAACGT-1",
  "TGFbeta1-3day-batch2_ACGTCCTCAAAGGCGT-1",
  "notreated-batch2_ATAGGCTGTACGTGAG-1",
  "TGFbeta1-8day-batch1_GTAGGCCTCAACACTG-1",
  "notreated-batch1_GTGCAGCTCAACGGCC-1"
)

for (cell in cells_R0) {
  tryCatch({
    message("\nAnalisando célula R0: ", cell)
    analyze_cell_R0(cell, max_pcs = 6)
  }, error = function(e) {
    message("Erro ao analisar ", cell, " em R0: ", e$message)
  })
}

# 2. Análise para células R30
cells_R30 <- c(
  "TGFbeta1-4day-batch1_CACACAACAGGAACGT-1",
  "TGFbeta1-8day-batch1_CACTCCAAGCTGCAAG-1",
  "TGFbeta1-1day-batch2_GCCCGAATCTACTTCA-1",
  "TGFbeta1-1day-batch2_AGAAGCGGTCATCAGT-1",
  "TGFbeta1-2day-batch2_CGGAATTTCGATACTG-1",
  "notreated-batch1_GATCTAGGTGTAATGA-1"
)

for (cell in cells_R30) {
  tryCatch({
    message("\nAnalisando célula R30: ", cell)
    analyze_cell_R30(cell, max_pcs = 6)
  }, error = function(e) {
    message("Erro ao analisar ", cell, " em R30: ", e$message)
  })
}

# 3. Análise para célula compartilhada (se aplicável)
shared_cell <- "TGFbeta1-4day-batch1_CACACAACAGGAACGT-1"

# Verificar se a célula existe em ambos os conjuntos
load_required_data_R0()
load_required_data_R30()

if (shared_cell %in% colnames(t_matrix_R0@transcriptogramS2) && 
    shared_cell %in% colnames(t_matrix_R30@transcriptogramS2)) {
  tryCatch({
    message("\nAnalisando célula compartilhada em R0: ", shared_cell)
    analyze_cell_R0(shared_cell, max_pcs = 6)
    
    message("\nAnalisando célula compartilhada em R30: ", shared_cell)
    analyze_cell_R30(shared_cell, max_pcs = 6)
  }, error = function(e) {
    message("Erro ao analisar célula compartilhada ", shared_cell, ": ", e$message)
  })
} else {
  message("\nA célula ", shared_cell, " não está presente em ambos R0 e R30")
}