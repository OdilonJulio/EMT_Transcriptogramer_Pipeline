# Carregar pacotes com suppressMessages para evitar poluição
suppressMessages({
  library(transcriptogramer)
  library(ggplot2)
  library(gridExtra)
  library(ggfortify)
  library(reshape2)
  library(dplyr)
  library(matrixStats)  # Para operações matriciais rápidas
  library(doParallel)   # Para paralelização
  library(progress)     # Para barra de progresso
})

# Configurar núcleos para paralelização (deixar 1 core livre)
registerDoParallel(cores = max(1, detectCores() - 100))

# Função otimizada para carregar dados
load_required_data <- function() {
  required_files <- c("t_matrix_R0.RData", "pca_result_R0.RData", 
                      "t_matrix_R30.RData", "pca_result_R30.RData")
  
  if (!all(file.exists(required_files))) {
    missing <- setdiff(required_files, list.files())
    stop("Arquivos faltando: ", paste(missing, collapse = ", "))
  }
  
  # Barra de progresso para carregamento
  pb <- progress_bar$new(
    format = "Carregando dados [:bar] :percent (:eta estimado)",
    total = length(required_files), clear = FALSE, width = 60)
  
  # Carregar com suppressMessages para limpeza
  suppressMessages({
    load("t_matrix_R0.RData", envir = .GlobalEnv)
    pb$tick()
    load("pca_result_R0.RData", envir = .GlobalEnv)
    pb$tick()
    load("t_matrix_R30.RData", envir = .GlobalEnv)
    pb$tick()
    load("pca_result_R30.RData", envir = .GlobalEnv)
    pb$tick()
  })
  
  message("\nDados carregados com sucesso")
}

# Função otimizada de análise com barra de progresso
analyze_transcriptogram_set <- function(t_matrix, pca_result, condition_label = "R30", 
                                        marker_genes = c("GENE1", "GENE2", "GENE3"),
                                        fast_mode = TRUE) {
  
  # Verificação rápida de classes
  if (!inherits(t_matrix, "Transcriptogram")) {
    stop("Objeto t_matrix deve ser da classe Transcriptogram")
  }
  
  # Configurar barra de progresso principal
  main_pb <- progress_bar$new(
    format = paste0(condition_label, " [:bar] :percent (:eta estimado)"),
    total = 13, clear = FALSE, width = 60)  # 13 passos totais
  
  # Extrair dados de forma mais eficiente
  main_pb$tick(0, tokens = list(step = "Extraindo dados"))
  expr_data_full <- as.matrix(t_matrix@transcriptogramS2[, !colnames(t_matrix@transcriptogramS2) %in% c("Protein", "Position")])
  position <- t_matrix@transcriptogramS2$Position
  main_pb$tick(1)
  
  # Filtrar células com variância zero - método rápido
  main_pb$tick(0, tokens = list(step = "Filtrando células"))
  col_vars <- matrixStats::colVars(expr_data_full, na.rm = TRUE)
  valid_cells <- col_vars > 0 & !is.na(col_vars)
  expr_data <- expr_data_full[, valid_cells, drop = FALSE]
  main_pb$tick(1)
  
  if (ncol(expr_data) < 3) {
    stop("Número insuficiente de células válidas (n < 3)")
  }
  
  # 1. Cálculo de variabilidade (mais rápido com matrixStats)
  main_pb$tick(0, tokens = list(step = "Calculando variabilidade"))
  variability <- setNames(matrixStats::colSds(expr_data, na.rm = TRUE), colnames(expr_data))
  main_pb$tick(1)
  
  # 2. Médias (otimizado)
  main_pb$tick(0, tokens = list(step = "Calculando médias"))
  means <- colMeans(expr_data, na.rm = TRUE)
  outliers_mean <- means[means > mean(means) + 2 * sd(means) | means < mean(means) - 2 * sd(means)]
  main_pb$tick(1)
  
  # 3. PCA (já pré-calculado)
  main_pb$tick(0, tokens = list(step = "Processando PCA"))
  scores <- pca_result$x + colMeans(pca_result$x) 
  pca_outliers <- rownames(scores)[abs(scores[, 1]) > 3 * sd(scores[, 1]) | abs(scores[, 2]) > 3 * sd(scores[, 2])]
  main_pb$tick(1)
  
  # ... (parte anterior da função)
  
  # 4. Clusterização (paralelizada)
  main_pb$tick(0, tokens = list(step = "Clusterizando"))
  d <- dist(t(expr_data[sample(nrow(expr_data), min(1000, nrow(expr_data))), ]))  # fechar parênteses do dist
  hc <- hclust(d)
  clust_groups <- cutree(hc, k = 5)
  small_clusters <- which(table(clust_groups) < 0.1 * length(clust_groups))
  clust_outliers <- names(clust_groups)[clust_groups %in% small_clusters]
  main_pb$tick(1)
  
  # 5. Marcadores genéticos
  main_pb$tick(0, tokens = list(step = "Analisando marcadores"))
  marker_expr <- colSums(expr_data[rownames(expr_data) %in% marker_genes, , drop = FALSE], na.rm = TRUE)
  top_marker_cells <- names(sort(marker_expr, decreasing = TRUE))[1:10]
  main_pb$tick(1)
  
  # 6. Score combinado (PARTE CRÍTICA - SELECIONA AS 6 CÉLULAS)
  main_pb$tick(0, tokens = list(step = "Calculando score combinado"))
  combined_score <- scale(variability) + scale(means) + 
    scale(abs(scores[,1])) + scale(abs(scores[,2]))
  combined_score <- setNames(as.vector(combined_score), colnames(expr_data))
  
  # Ordenar por score combinado e selecionar as 6 células
  sorted_cells <- names(sort(combined_score, decreasing = TRUE))
  selected_cells <- c(
    head(sorted_cells, 2),          # 2 com maior score
    tail(sorted_cells, 2),          # 2 com menor score
    sorted_cells[round(length(sorted_cells)/2):(round(length(sorted_cells)/2)+1)] # 2 do meio
  )
  main_pb$tick(1)
  
  # 7. Preparação de resultados finais
  main_pb$tick(0, tokens = list(step = "Preparando resultados"))
  final_results <- list(
    variability = variability,
    mean_expression = means,
    pca_outliers = pca_outliers,
    cluster_outliers = clust_outliers,
    top_marker_cells = top_marker_cells,
    top_combined_outliers = combined_score,
    selected_cells = selected_cells,  # Inclui as 6 células selecionadas
    top_variance_cell = names(which.max(matrixStats::colVars(expr_data))),
    position = position,
    expr_data = expr_data
  )
  main_pb$tick(1)
  
  # 8. Plotagens (opcional no modo rápido)
  if (!fast_mode) {
    main_pb$tick(0, tokens = list(step = "Gerando gráficos"))
    plots <- generate_plots(final_results, condition_label)
    final_results$plots <- plots
    main_pb$tick(1)
  } else {
    final_results$plots <- NULL
    main_pb$tick(1, tokens = list(step = "Pulando gráficos"))
  }
  
  return(final_results)
  
  # ... (continua)
  
}

# Função auxiliar para gerar gráficos
generate_plots <- function(results, condition_label) {
  # [Implementação dos gráficos como no código original]
  # Retorna lista de plots
  return(list(
    plot_variability = NULL,
    plot_means = NULL,
    # ... outros plots
  ))
}

# Função de exportação otimizada
export_selected_cells <- function(t_matrix, analysis_results, label = "R0") {
  pb <- progress_bar$new(
    format = paste0("Exportando ", label, " [:bar] :percent"),
    total = 4, clear = FALSE, width = 60)
  
  pb$tick(0, tokens = list(step = "Preparando dados"))
  expr_df <- t_matrix@transcriptogramS2
  
  # Usar as células já selecionadas na análise
  selected <- analysis_results$selected_cells
  pb$tick(1)
  
  # Verificar células disponíveis
  available_cells <- intersect(selected, colnames(expr_df))
  pb$tick(1)
  
  if (length(available_cells) == 0) {
    warning("Nenhuma célula selecionada encontrada nos dados")
    pb$tick(2)
    return(invisible(NULL))
  }
  
  # Preparar dados para exportação
  output <- cbind(Position = expr_df$Position, 
                  expr_df[, available_cells, drop = FALSE])
  pb$tick(1)
  
  # Exportar
  file_name <- paste0("celulas_selecionadas_", label, ".csv")
  tryCatch({
    write.csv(output, file = file_name, row.names = FALSE)
    message("\nArquivo salvo: ", file_name, 
            "\nCélulas exportadas: ", paste(available_cells, collapse = ", "))
    pb$tick(1)
  }, error = function(e) {
    warning("Falha ao exportar: ", e$message)
    pb$tick(1)
  })
  
  return(invisible(file_name))
}

# Função principal completa
main <- function(fast_mode = TRUE) {
  start_time <- Sys.time()
  
  # Barra de progresso geral
  overall_pb <- progress_bar$new(
    format = "Progresso total [:bar] :percent (:eta estimado)",
    total = 6, clear = FALSE, width = 60)
  
  tryCatch({
    # 1. Carregar dados
    overall_pb$tick(0, tokens = list(step = "Carregando dados"))
    load_required_data()
    overall_pb$tick(1)
    
    results <- list()
    
    # Processar R0
    overall_pb$tick(0, tokens = list(step = "Analisando R0"))
    results$R0 <- analyze_transcriptogram_set(t_matrix_R0, pca_result_R0$pca_result, "R0", fast_mode)
    overall_pb$tick(1)
    
    # Processar R30
    overall_pb$tick(0, tokens = list(step = "Analisando R30"))
    results$R30 <- analyze_transcriptogram_set(t_matrix_R30, pca_result_R30$pca_result, "R30", fast_mode)
    overall_pb$tick(1)
    
    # Exportar R0
    overall_pb$tick(0, tokens = list(step = "Exportando R0"))
    export_selected_cells(t_matrix_R0, results$R0, "R0")
    overall_pb$tick(1)
    
    # Exportar R30
    overall_pb$tick(0, tokens = list(step = "Exportando R30"))
    export_selected_cells(t_matrix_R30, results$R30, "R30")
    overall_pb$tick(1)
    
    # Finalização
    elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 1)
    message("\nAnálise concluída em ", elapsed, " minutos")
    
    return(results)
    
  }, error = function(e) {
    message("\nErro na execução: ", e$message)
    return(NULL)
  })
}

# Executar análise completa
analysis_results <- main(fast_mode = TRUE)

# Para R0
scores_r0 <- analysis_results[["R0"]][["top_combined_outliers"]]
ordered_r0 <- sort(scores_r0, decreasing = TRUE)

# Todas as células ordenadas por score
all_cells_r0 <- names(ordered_r0)

# Células selecionadas (já sabemos quais são)
selected_r0 <- analysis_results[["R0"]][["selected_cells"]]

# Identificar posições no ranking
positions_r0 <- sapply(selected_r0, function(x) which(all_cells_r0 == x))

# Classificar
classification_r0 <- ifelse(positions_r0 <= 2, "Top",
                            ifelse(positions_r0 >= (length(all_cells_r0)-1), "Bottom", "Middle"))

# Resultado para R0
result_r0 <- data.frame(
  Cell = selected_r0,
  Position_in_Ranking = positions_r0,
  Category = classification_r0,
  Score = scores_r0[selected_r0]
)


# Para R30
scores_r30 <- analysis_results[["R30"]][["top_combined_outliers"]]
ordered_r30 <- sort(scores_r30, decreasing = TRUE)

# Todas as células ordenadas por score
all_cells_r30 <- names(ordered_r30)

# Células selecionadas
selected_r30 <- analysis_results[["R30"]][["selected_cells"]]

# Identificar posições no ranking
positions_r30 <- sapply(selected_r30, function(x) which(all_cells_r30 == x))

# Classificar
classification_r30 <- ifelse(positions_r30 <= 2, "Top",
                             ifelse(positions_r30 >= (length(all_cells_r30)-1), "Bottom", "Middle"))

# Resultado para R30
result_r30 <- data.frame(
  Cell = selected_r30,
  Position_in_Ranking = positions_r30,
  Category = classification_r30,
  Score = scores_r30[selected_r30]
)
