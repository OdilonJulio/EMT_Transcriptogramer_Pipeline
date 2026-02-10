# transcriptogram_reconstruction_validation_parallel.R

library(transcriptogramer, lib.loc="/home/ojdsantos/R/rocker-rstudio")
library(dplyr, lib.loc="/home/ojdsantos/R/rocker-rstudio")
library(ggplot2, lib.loc="/home/ojdsantos/R/rocker-rstudio")
library(parallel, lib.loc="/home/ojdsantos/R/rocker-rstudio")
library(progress, lib.loc="/home/ojdsantos/R/rocker-rstudio")

## ===============================
## 1. Função para carregar dados
## ===============================
load_data <- function() {
  files <- c(
    "reconstructed_matrix_R0.RData", "reconstructed_matrix_R30.RData",
    "t_matrix_R0.RData", "t_matrix_R30.RData",
    "pca_result_R0.RData", "pca_result_R30.RData"
  )
  
  for (file in files) {
    if (!file.exists(file)) stop("File not found: ", file)
    load(file, envir = .GlobalEnv)
    gc()
  }
}

## ==========================================
## 2. Preparar matriz original
## ==========================================
prepare_original_matrix <- function(obj) {
  mat <- as.matrix(obj@transcriptogramS2[, !names(obj@transcriptogramS2) %in% c("Protein", "Position")])
  rownames(mat) <- obj@transcriptogramS2$Protein
  gc()
  t(mat)
}

## ==========================================
## 3. Função para calcular R²
## ==========================================
calculate_r_squared <- function(orig, recon) {
  gene_means <- rowMeans(orig, na.rm = TRUE)
  mode_val <- as.numeric(names(which.max(table(orig)))) / 1000
  zero_means <- gene_means == 0
  gene_means[zero_means] <- mode_val
  
  chunk_size <- 1000
  n_genes <- nrow(orig)
  result <- numeric(n_genes)
  
  for (i in seq(1, n_genes, chunk_size)) {
    end <- min(i + chunk_size - 1, n_genes)
    chunk <- (orig[i:end, ] - recon[i:end, ])^2 / (gene_means[i:end]^2)
    result[i:end] <- rowMeans(chunk)
  }
  
  list(
    global = mean(result),
    by_gene = result,
    by_cell = colMeans((orig - recon)^2 / (gene_means^2)),
    adjusted_zero_genes = sum(zero_means),
    adjustment_value = mode_val
  )
}

## ==========================================
## 4. Função para reconstruir PCA
## ==========================================
reconstruct_with_n_pcs <- function(pca_res, n_pcs) {
  rot <- pca_res$pca_result$rotation[, 1:n_pcs, drop = FALSE]
  comp <- pca_res$pca_result$x[, 1:n_pcs, drop = FALSE]
  ctr <- pca_res$pca_result$center
  
  recon <- tcrossprod(comp, rot)
  if (!is.null(ctr)) recon <- sweep(recon, 2, ctr, "+")
  recon
}

## ==========================================
## 5. Paralelização com mclapply
## ==========================================
analyze_pc_components_mclapply <- function(pca_res, orig_mat, label, n_cores = 50) {
  max_pcs <- ncol(pca_res$pca_result$x)
  message(paste0("Starting parallel analysis for ", label, " with ", max_pcs, " PCs..."))
  
  results <- mclapply(1:max_pcs, function(n) {
    recon <- reconstruct_with_n_pcs(pca_res, n)
    r2 <- calculate_r_squared(orig_mat, recon)$global
    r2
  }, mc.cores = n_cores)
  
  results <- unlist(results)
  
  plot <- ggplot(data.frame(n_pcs = 1:max_pcs, r_squared = results),
                 aes(n_pcs, r_squared)) +
    geom_line(color = "steelblue") +
    geom_point(color = "steelblue", size = 1) +
    labs(title = paste("Reconstruction Error by PCs -", label),
         x = "Principal Components", y = "R² Error") +
    theme_minimal()
  
  list(results = data.frame(n_pcs = 1:max_pcs, r_squared = results), plot = plot)
}

## ==========================================
## 5. Paralelização com mclapply e barra de progresso
## ==========================================
analyze_pc_components_mclapply <- function(pca_res, orig_mat, label, n_cores = 50) {
  max_pcs <- ncol(pca_res$pca_result$x)
  message(paste0("Starting parallel analysis for ", label, " with ", max_pcs, " PCs..."))
  
  # Cria uma barra de progresso
  pb <- progress_bar$new(
    format = paste(label, "[:bar] :percent (:eta)"),
    total = max_pcs, 
    clear = FALSE,
    width = 60
  )
  
  # Função wrapper para atualizar o progresso
  process_pc <- function(n) {
    recon <- reconstruct_with_n_pcs(pca_res, n)
    r2 <- calculate_r_squared(orig_mat, recon)$global
    pb$tick()  # Atualiza a barra de progresso
    return(r2)
  }
  
  # Execução paralela
  results <- mclapply(1:max_pcs, process_pc, mc.cores = n_cores)
  
  results <- unlist(results)
  
  plot <- ggplot(data.frame(n_pcs = 1:max_pcs, r_squared = results),
                 aes(n_pcs, r_squared)) +
    geom_line(color = "steelblue") +
    geom_point(color = "steelblue", size = 1) +
    labs(title = paste("Reconstruction Error by PCs -", label),
         x = "Principal Components", y = "R² Error") +
    theme_minimal()
  
  list(results = data.frame(n_pcs = 1:max_pcs, r_squared = results), plot = plot)
}

## ==========================================
## 6. Função para gerar Elbow Plot (atualizada)
## ==========================================
generate_elbow_plots <- function(pca_res, r2_res, label) {
  sdev <- pca_res$pca_result$sdev[1:nrow(r2_res)]
  df <- data.frame(
    PC = r2_res$n_pcs,
    R2 = r2_res$r_squared / max(r2_res$r_squared),
    CumVar = cumsum(sdev^2) / sum(pca_res$pca_result$sdev^2)
  )
  
  # Calculate elbow point
  fit1 <- lm(CumVar ~ PC, df[1:2, ])
  fit2 <- lm(CumVar ~ PC, tail(df, 2))
  elbow_x <- (coef(fit2)[1] - coef(fit1)[1]) / (coef(fit1)[2] - coef(fit2)[2])
  elbow_y <- coef(fit1)[2] * elbow_x + coef(fit1)[1]
  
  # Define os intervalos de PCs para plotar
  pc_breaks <- list(
    "10 PCs" = 10,
    "50 PCs" = 50,
    "100 PCs" = 100,
    "All PCs" = nrow(df)
  )
  
  # Função para criar cada plot individual
  create_single_plot <- function(max_pc, label_suffix) {
    plot_data <- df[1:max_pc, ]
    
    ggplot(plot_data, aes(factor(PC))) +
      # Main components
      geom_col(aes(y = R2 * CumVar[1]), fill = "#1f77b4", alpha = 0.6, width = 0.7) +
      geom_line(aes(y = CumVar, group = 1), color = "#ff7f0e", linewidth = 1.2) +
      geom_point(aes(y = CumVar), color = "#ff7f0e", size = 2) +
      
      # Elbow point (se visível no intervalo)
      {if(elbow_x <= max_pc) 
        list(
          annotate("segment",
                   x = 1, xend = elbow_x,
                   y = predict(fit1, data.frame(PC = 1)),
                   yend = elbow_y,
                   color = "red", linetype = "dashed"),
          annotate("point",
                   x = elbow_x,
                   y = elbow_y,
                   color = "#2ca02c", size = 4, shape = 18)
        )} +
      
      # Labels and theme
      labs(title = paste("Component Analysis -", label, label_suffix),
           x = "Principal Components",
           y = "Normalized R² Error / Cumulative Variance",
           caption = ifelse(max_pc == nrow(df),
                            paste("Elbow point at", round(elbow_x, 1), "PCs (", 
                                  round(elbow_y*100, 1), "% variance)"),
                            "")) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        panel.grid.minor = element_blank(),
        plot.caption = element_text(hjust = 0.5, face = "italic"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,1), "cm")
      ) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      scale_x_discrete(breaks = seq(0, max_pc, by = ifelse(max_pc <= 100, 10, 1000)))
  }
  
  # Criar todos os plots
  plots <- lapply(names(pc_breaks), function(name) {
    max_pc <- min(pc_breaks[[name]], nrow(df))
    create_single_plot(max_pc, paste0("(", name, ")"))
  })
  
  # Combinar os plots em uma única figura
  combined_plot <- ggpubr::ggarrange(plotlist = plots, ncol = 2, nrow = 2)
  
  # Salvar individualmente cada plot com alta resolução
  lapply(names(pc_breaks), function(name) {
    max_pc <- min(pc_breaks[[name]], nrow(df))
    p <- create_single_plot(max_pc, paste0("(", name, ")"))
    ggsave(paste0("images/elbow_plot_", label, "_", gsub(" ", "", name), ".png"), 
           plot = p, width = 8, height = 6, dpi = 300, bg = "white")
  })
  
  # Retornar o plot combinado e individual
  list(
    combined_plot = combined_plot,
    individual_plots = plots
  )
}


## ==========================================
## 7. MAIN EXECUTION
## ==========================================
load_data()
cat("\nPreparing matrices...\n")
original_matrix_R0 <- prepare_original_matrix(t_matrix_R0)
original_matrix_R30 <- prepare_original_matrix(t_matrix_R30)
rm(t_matrix_R0, t_matrix_R30)
gc()

## Defina o número de núcleos que deseja usar
available_cores <- parallel::detectCores()
n_cores_to_use <- min(50, available_cores)  # Exemplo: até 50 cores

## Análises paralelas
cat("\nRunning analyses...\n")
analysis_R0 <- analyze_pc_components_mclapply(pca_result_R0, original_matrix_R0, "R0", n_cores = n_cores_to_use)
save(analysis_R0, file = "analysis_R0.RData")

analysis_R30 <- analyze_pc_components_mclapply(pca_result_R30, original_matrix_R30, "R30", n_cores = n_cores_to_use)
save(analysis_R30, file = "analysis_R30.RData")

## Avaliação da reconstrução completa
full_recon_r2_R0 <- calculate_r_squared(original_matrix_R0, reconstructed_matrix_R0)
full_recon_r2_R30 <- calculate_r_squared(original_matrix_R30, reconstructed_matrix_R30)

## Salvamento dos resultados
save(analysis_R0, analysis_R30, full_recon_r2_R0, full_recon_r2_R30,
     file = "transcriptogram_validation_results_parallel.RData")


## ==========================================
## Atualização da execução principal - gráficos
## ==========================================

elbow_results_R0 <- generate_elbow_plots(pca_result_R0, analysis_R0$results, "R0")
elbow_results_R30 <- generate_elbow_plots(pca_result_R30, analysis_R30$results, "R30")

# Salvar o plot combinado
ggsave("images/elbow_plot_combined_R0.png", plot = elbow_results_R0$combined_plot, 
       width = 12, height = 10, dpi = 300, bg = "white")
ggsave("images/elbow_plot_combined_R30.png", plot = elbow_results_R30$combined_plot, 
       width = 12, height = 10, dpi = 300, bg = "white")


# Save each plot individually with explicit filenames
ggsave("images/reconstruction_error_R0.png", analysis_R0$plot, width = 8, height = 6)
ggsave("imagesreconstruction_error_R30.png", analysis_R30$plot, width = 8, height = 6)

# Alternative approach if you want to keep the list structure
plot_list <- list(
  reconstruction_error_R0 = analysis_R0$plot,
  reconstruction_error_R30 = analysis_R30$plot,
  elbow_plot_R0 = elbow_R0,
  elbow_plot_R30 = elbow_R30
)

# Save plots using a loop
for (plot_name in names(plot_list)) {
  ggsave(
    filename = paste0("images/",plot_name, ".png"),
    plot = plot_list[[plot_name]],
    width = 8,
    height = 6
  )
}
