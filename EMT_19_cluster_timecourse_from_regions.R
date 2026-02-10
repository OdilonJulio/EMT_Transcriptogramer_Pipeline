############################################
### IDENTIFICAÇÃO DE CLUSTERS E ANÁLISE TEMPORAL ###
############################################

# Configuração inicial
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(stringr)
library(RColorBrewer)
library(zoo)
library(scales)

# Diretórios
output_dir <- "cluster_analysis"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "images"), showWarnings = FALSE, recursive = TRUE)

############################################
### 1. CARREGAR E PREPARAR DADOS BASE ####
############################################

# Carregar dados necessários
load("~/mestrado/pca_result_R30.RData")
load("~/mestrado/t_matrix_R30.RData")

# Função robusta para carregar dados com fallbacks
carregar_dados_regioes <- function() {
  regions_dir <- "regions_data"
  
  # Verificar se o diretório existe
  if (!dir.exists(regions_dir)) {
    stop("Diretório 'regions_data' não encontrado. Execute primeiro o pipeline de regiões.")
  }
  
  # Tentar encontrar o arquivo mais recente combinado
  arquivos_combinados <- list.files(regions_dir, 
                                    pattern = "high_variation_autoRegions_combined.*\\.csv$",
                                    full.names = TRUE)
  
  if (length(arquivos_combinados) > 0) {
    # Usar o arquivo mais recente
    arquivo_mais_recente <- arquivos_combinados[which.max(file.info(arquivos_combinados)$mtime)]
    message("Carregando arquivo combinado mais recente: ", basename(arquivo_mais_recente))
    return(read.csv(arquivo_mais_recente, stringsAsFactors = FALSE))
  }
  
  # Fallback: tentar carregar arquivos individuais de regiões
  message("Arquivo combinado não encontrado. Buscando arquivos individuais de regiões...")
  arquivos_individuais <- list.files(regions_dir, 
                                     pattern = "high_variation_region_.*\\.csv$",
                                     full.names = TRUE)
  
  if (length(arquivos_individuais) == 0) {
    # Último fallback: gerar regiões automaticamente a partir dos dados PCA
    message("Nenhum arquivo de regiões encontrado. Gerando regiões automaticamente...")
    return(gerar_regioes_automaticas())
  }
  
  # Combinar arquivos individuais
  message("Combinando ", length(arquivos_individuais), " arquivos de regiões...")
  regions_combined <- combinar_arquivos_regioes(arquivos_individuais)
  return(regions_combined)
}

# Função para combinar arquivos individuais de regiões
combinar_arquivos_regioes <- function(arquivos) {
  dados_regioes <- map_dfr(arquivos, function(arquivo) {
    df <- read.csv(arquivo, stringsAsFactors = FALSE)
    # Extrair nome da região do nome do arquivo
    nome_regiao <- str_extract(basename(arquivo), "AutoRegion_[^\\\\.]+")
    if (is.na(nome_regiao)) {
      nome_regiao <- str_replace(basename(arquivo), "\\.csv$", "")
    }
    df$region_name <- nome_regiao
    df$region_source_file <- basename(arquivo)
    return(df)
  })
  return(dados_regioes)
}

# Função de último recurso: gerar regiões automaticamente
gerar_regioes_automaticas <- function() {
  message("Gerando regiões de alta variação automaticamente a partir dos dados PCA...")
  
  # Preparar dados de variação
  gene_names <- t_matrix_R30@transcriptogramS2[, 1]
  gene_positions <- t_matrix_R30@transcriptogramS2[, 2]
  pca_rotation <- pca_result_R30[["pca_result"]][["rotation"]]
  rownames(pca_rotation) <- gene_names
  
  n_pcs <- 6
  variation_signal <- rowSums(abs(pca_rotation[, 1:n_pcs, drop = FALSE]))
  
  var_df <- data.frame(
    Gene = gene_names,
    Position = gene_positions,
    Variation = variation_signal,
    stringsAsFactors = FALSE
  )
  
  # Suavizar o sinal
  window_size <- 15
  var_df$VariationSmooth <- zoo::rollmean(var_df$Variation, k = window_size, 
                                          fill = NA, align = "center")
  
  # Identificar regiões acima da média
  threshold <- mean(var_df$VariationSmooth, na.rm = TRUE)
  var_df$AboveMean <- var_df$VariationSmooth > threshold
  
  # Encontrar regiões contíguas
  rle_above <- rle(var_df$AboveMean)
  ends <- cumsum(rle_above$lengths)
  starts <- c(1, head(ends, -1) + 1)
  
  region_tbl <- data.frame(
    start_idx = starts[rle_above$values],
    end_idx = ends[rle_above$values]
  )
  
  # Filtrar regiões muito pequenas
  min_size <- 20
  region_tbl <- region_tbl[region_tbl$end_idx - region_tbl$start_idx + 1 >= min_size, ]
  
  if (nrow(region_tbl) == 0) {
    stop("Não foi possível identificar regiões de alta variação automaticamente.")
  }
  
  # Construir dados finais
  regions_auto <- map_dfr(1:nrow(region_tbl), function(i) {
    start_pos <- var_df$Position[region_tbl$start_idx[i]]
    end_pos <- var_df$Position[region_tbl$end_idx[i]]
    
    region_data <- var_df %>%
      filter(Position >= start_pos & Position <= end_pos) %>%
      mutate(
        region_name = paste0("AutoRegion_", i),
        region_source_file = "generated_automatically.csv"
      )
    
    return(region_data)
  })
  
  message("Geradas ", nrow(region_tbl), " regiões automaticamente.")
  return(regions_auto)
}

# Função robusta para carregar dados do transcriptograma
carregar_dados_transcriptograma <- function() {
  transcriptogram_file <- "transcriptograms_reconstructed_by_PC1_weighted.csv"
  
  if (!file.exists(transcriptogram_file)) {
    stop("Arquivo do transcriptograma reconstruído não encontrado: ", transcriptogram_file,
         "\nExecute primeiro o script da animação para gerar este arquivo.")
  }
  
  transcriptogram_data <- read.csv(transcriptogram_file, stringsAsFactors = FALSE)
  
  # Validar estrutura do arquivo
  required_cols <- c("Gene", "Expression", "PC1_value", "IntervalID")
  missing_cols <- setdiff(required_cols, colnames(transcriptogram_data))
  
  if (length(missing_cols) > 0) {
    stop("Arquivo do transcriptograma incompleto. Colunas faltantes: ", 
         paste(missing_cols, collapse = ", "))
  }
  
  message("Dados do transcriptograma carregados: ", 
          nrow(transcriptogram_data), " registros, ",
          length(unique(transcriptogram_data$PC1_value)), " valores de PC1")
  
  return(transcriptogram_data)
}

# Carregar dados com tratamento de erro robusto
message("=== CARREGAMENTO DE DADOS ===")
tryCatch({
  regions_combined <- carregar_dados_regioes()
  message("✅ Dados das regiões carregados: ", nrow(regions_combined), " registros")
}, error = function(e) {
  stop("Falha no carregamento dos dados das regiões: ", e$message)
})

tryCatch({
  transcriptogram_data <- carregar_dados_transcriptograma()
}, error = function(e) {
  stop("Falha no carregamento dos dados do transcriptograma: ", e$message)
})

############################################
### 2. IDENTIFICAÇÃO ROBUSTA DE CLUSTERS ###
############################################

identificar_clusters <- function(regions_df, min_genes = 5, max_gap = 2) {
  
  # Validar dados de entrada
  if (nrow(regions_df) == 0) {
    stop("Dataframe de regiões vazio.")
  }
  
  # Garantir que temos as colunas necessárias
  required_cols <- c("Gene", "Position")
  missing_cols <- setdiff(required_cols, colnames(regions_df))
  if (length(missing_cols) > 0) {
    stop("Colunas necessárias não encontradas: ", paste(missing_cols, collapse = ", "))
  }
  
  # Preparar dados
  regions_clean <- regions_df %>%
    filter(!is.na(Position)) %>%
    arrange(Position) %>%
    distinct(Gene, .keep_all = TRUE)
  
  if (nrow(regions_clean) == 0) {
    stop("Nenhum dado válido após limpeza.")
  }
  
  message("Identificando clusters em ", nrow(regions_clean), " genes...")
  
  # Se não tiver VariationSmooth, calcular a partir da posição
  if (!"VariationSmooth" %in% colnames(regions_clean)) {
    message("VariationSmooth não encontrada. Usando densidade de posições para clustering...")
    regions_clean$VariationSmooth <- 1  # Valor constante para clustering por posição
  }
  
  # Estratégia 1: Agrupamento por densidade de posições
  perform_position_clustering <- function(data, n_clusters = NULL) {
    if (is.null(n_clusters)) {
      # Estimar número ideal de clusters
      n_genes <- nrow(data)
      n_clusters <- max(2, min(8, floor(n_genes / 30)))
    }
    
    if (nrow(data) <= n_clusters) {
      # Caso especial: poucos genes
      clusters <- rep(1, nrow(data))
    } else {
      # K-means nas posições
      set.seed(123)  # Para reproducibilidade
      position_matrix <- matrix(data$Position, ncol = 1)
      kmeans_result <- kmeans(position_matrix, centers = n_clusters)
      clusters <- kmeans_result$cluster
    }
    
    return(clusters)
  }
  
  # Estratégia 2: Agrupamento hierárquico
  perform_hierarchical_clustering <- function(data, max_gap) {
    if (nrow(data) <= 1) return(rep(1, nrow(data)))
    
    distance_matrix <- dist(matrix(data$Position, ncol = 1))
    hc <- hclust(distance_matrix, method = "complete")
    
    # Cortar a árvore baseado no gap máximo
    clusters <- cutree(hc, h = max_gap)
    return(clusters)
  }
  
  # Tentar diferentes estratégias
  message("Executando agrupamento por densidade...")
  position_clusters <- perform_position_clustering(regions_clean)
  
  # Criar dataframe com clusters
  regions_clustered <- regions_clean %>%
    mutate(
      cluster_id = paste0("Cluster_", position_clusters),
      Cluster = factor(cluster_id, levels = paste0("Cluster_", sort(unique(position_clusters))))
    )
  
  # Filtrar clusters muito pequenos
  cluster_summary <- regions_clustered %>%
    group_by(Cluster) %>%
    summarise(
      n_genes = n(),
      mean_position = mean(Position),
      span = max(Position) - min(Position),
      .groups = "drop"
    ) %>%
    filter(n_genes >= min_genes)
  
  if (nrow(cluster_summary) == 0) {
    warning("Nenhum cluster atende ao critério mínimo de genes. Reduzindo min_genes para 3.")
    cluster_summary <- regions_clustered %>%
      group_by(Cluster) %>%
      summarise(
        n_genes = n(),
        mean_position = mean(Position),
        span = max(Position) - min(Position),
        .groups = "drop"
      ) %>%
      filter(n_genes >= 3)
  }
  
  regions_final <- regions_clustered %>%
    filter(Cluster %in% cluster_summary$Cluster) %>%
    select(Gene, Position, any_of(c("Variation", "VariationSmooth")), 
           Cluster, any_of(c("region_name", "region_source_file")))
  
  # Reordenar clusters por posição
  cluster_order <- cluster_summary %>%
    arrange(mean_position) %>%
    pull(Cluster) %>%
    as.character()
  
  regions_final$Cluster <- factor(regions_final$Cluster, levels = cluster_order)
  
  message("✅ Clusters identificados: ", nrow(cluster_summary))
  for (i in 1:nrow(cluster_summary)) {
    message("   ", cluster_summary$Cluster[i], ": ", cluster_summary$n_genes, 
            " genes (pos: ", round(cluster_summary$mean_position[i], 1), ")")
  }
  
  return(regions_final)
}

# Aplicar identificação de clusters
message("\n=== IDENTIFICAÇÃO DE CLUSTERS ===")
gene_clusters <- identificar_clusters(regions_combined, min_genes = 8, max_gap = 3)

# Salvar clusters identificados
write.csv(gene_clusters, 
          file.path(output_dir, "identified_gene_clusters.csv"), 
          row.names = FALSE)
message("✅ Clusters salvos em: ", file.path(output_dir, "identified_gene_clusters.csv"))

############################################
### 3. PREPARAR DADOS TEMPORAIS #########
############################################

preparar_dados_temporais <- function(transcriptogram_df, clusters_df) {
  
  # Validar dados de entrada
  if (nrow(transcriptogram_df) == 0 || nrow(clusters_df) == 0) {
    stop("Dados de entrada vazios para preparação temporal.")
  }
  
  # Extrair valores únicos de PC1 (ordenados)
  pc1_values <- unique(transcriptogram_df$PC1_value) 
  if (length(pc1_values) == 0) {
    stop("Nenhum valor de PC1 encontrado no transcriptograma.")
  }
  pc1_values <- sort(pc1_values)
  
  message("Preparando dados temporais para ", length(pc1_values), " valores de PC1...")
  
  # Mapear genes para clusters
  gene_to_cluster <- setNames(as.character(clusters_df$Cluster), clusters_df$Gene)
  
  # Calcular importância por cluster ao longo do tempo
  importance_data <- map_dfr(pc1_values, function(pc1_val) {
    
    # Dados do transcriptograma para este valor de PC1
    frame_data <- transcriptogram_df %>%
      filter(PC1_value == pc1_val) %>%
      select(Gene, Expression)
    
    if (nrow(frame_data) == 0) {
      return(NULL)
    }
    
    # Adicionar informação de cluster e calcular importância
    frame_with_clusters <- frame_data %>%
      mutate(Cluster = gene_to_cluster[Gene]) %>%
      filter(!is.na(Cluster)) %>%
      group_by(Cluster) %>%
      summarise(
        Importance = sum(abs(Expression), na.rm = TRUE),
        n_genes = n(),
        .groups = "drop"
      ) %>%
      mutate(PC1_value = pc1_val)
    
    return(frame_with_clusters)
  })
  
  if (nrow(importance_data) == 0) {
    stop("Nenhum dado de importância calculado. Verifique o mapeamento genes-clusters.")
  }
  
  # Preencher missing values com 0
  all_combinations <- expand.grid(
    Cluster = unique(importance_data$Cluster),
    PC1_value = pc1_values,
    stringsAsFactors = FALSE
  )
  
  importance_complete <- all_combinations %>%
    left_join(importance_data, by = c("Cluster", "PC1_value")) %>%
    mutate(
      Importance = ifelse(is.na(Importance), 0, Importance),
      n_genes = ifelse(is.na(n_genes), 0, n_genes)
    )
  
  message("✅ Dados temporais preparados: ", 
          nrow(importance_complete), " registros, ",
          length(unique(importance_complete$Cluster)), " clusters")
  
  return(importance_complete)
}

message("\n=== PREPARAÇÃO DE DADOS TEMPORAIS ===")
temporal_data <- preparar_dados_temporais(transcriptogram_data, gene_clusters)

############################################
### 4. SUAVIZAÇÃO TEMPORAL ###############
############################################

suavizar_dados_temporais <- function(temporal_df, window_size = 5) {
  
  if (nrow(temporal_df) == 0) {
    stop("Dados temporais vazios para suavização.")
  }
  
  message("Aplicando suavização temporal (janela = ", window_size, ")...")
  
  temporal_smoothed <- temporal_df %>%
    group_by(Cluster) %>%
    arrange(PC1_value) %>%
    mutate(
      Importance_smooth = zoo::rollmean(Importance, k = window_size, 
                                        fill = NA, align = "center"),
      Importance_normalized = Importance / max(Importance, na.rm = TRUE),
      Importance_smooth_norm = ifelse(
        all(is.na(Importance_normalized)), 
        NA,
        zoo::rollmean(Importance_normalized, k = window_size, fill = NA, align = "center")
      )
    ) %>%
    ungroup() %>%
    filter(!is.na(Importance_smooth))
  
  message("✅ Dados suavizados: ", nrow(temporal_smoothed), " registros")
  
  return(temporal_smoothed)
}

# Aplicar suavização
temporal_smoothed <- suavizar_dados_temporais(temporal_data, window_size = 7)

# [CONTINUA... O RESTO DO SCRIPT PERMANECE IGUAL]

############################################
### 5. GRÁFICO TEMPORAL PRINCIPAL #######
############################################

criar_grafico_temporal <- function(temporal_df, clusters_info) {
  
  # Preparar paleta de cores
  n_clusters <- length(unique(temporal_df$Cluster))
  cluster_colors <- brewer.pal(min(9, max(3, n_clusters)), "Set1")
  
  if (n_clusters > length(cluster_colors)) {
    cluster_colors <- colorRampPalette(cluster_colors)(n_clusters)
  }
  
  names(cluster_colors) <- sort(unique(temporal_df$Cluster))
  
  # Calcular estatísticas para anotações
  cluster_stats <- temporal_df %>%
    group_by(Cluster) %>%
    summarise(
      peak_importance = max(Importance_smooth, na.rm = TRUE),
      peak_time = PC1_value[which.max(Importance_smooth)],
      mean_importance = mean(Importance_smooth, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Gráfico principal
  p <- ggplot(temporal_df, aes(x = PC1_value, y = Importance_smooth, color = Cluster)) +
    geom_line(size = 1.2, alpha = 0.8) +
    geom_point(data = cluster_stats, 
               aes(x = peak_time, y = peak_importance, color = Cluster),
               size = 3, shape = 18) +
    scale_color_manual(values = cluster_colors) +
    labs(
      title = "Evolução Temporal da Importância dos Clusters de Genes",
      subtitle = paste("Baseado na reconstrução do transcriptograma ao longo do PC1",
                       "\nPontos destacados indicam picos de importância máxima"),
      x = "Valor do PC1 (Progressão Temporal)",
      y = "Importância do Cluster (Soma dos Módulos dos Loadings)",
      color = "Cluster"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.position = "right",
      panel.grid.major = element_line(color = "grey90", size = 0.2),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black"),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    scale_x_continuous(labels = scales::scientific) +
    scale_y_continuous(labels = scales::scientific)
  
  return(p)
}

message("\n=== CRIAÇÃO DE GRÁFICOS ===")
temporal_plot <- criar_grafico_temporal(temporal_smoothed, gene_clusters)

# Salvar gráfico principal
ggsave(file.path(output_dir, "images", "evolucao_temporal_clusters.png"),
       temporal_plot, width = 16, height = 9, dpi = 600)
message("✅ Gráfico temporal salvo: evolucao_temporal_clusters.png")
############################################
### 6. GRÁFICOS ADICIONAIS - CORRIGIDO ###
############################################

# Função para criar gráfico normalizado por frame - CORRIGIDA
criar_grafico_normalizado <- function(df, 
                                      value_col = "Importance_smooth", 
                                      cluster_col = "Cluster", 
                                      frame_col = "PC1_value") {
  
  # Verificar colunas disponíveis
  available_cols <- colnames(df)
  message("Colunas disponíveis no dataframe: ", paste(available_cols, collapse = ", "))
  
  # Verificar se as colunas necessárias existem
  required_cols <- c(value_col, cluster_col, frame_col)
  missing_cols <- setdiff(required_cols, available_cols)
  
  if (length(missing_cols) > 0) {
    stop("O data.frame não contém as colunas esperadas: ", 
         paste(missing_cols, collapse = ", "),
         "\nColunas disponíveis: ", paste(available_cols, collapse = ", "))
  }
  
  # Fazer uma cópia para não modificar o original
  df_plot <- df %>%
    mutate(
      Value = as.numeric(.data[[value_col]]),
      Frame = as.numeric(.data[[frame_col]]),
      Cluster = as.factor(.data[[cluster_col]])
    ) %>%
    filter(!is.na(Value), !is.na(Frame), !is.na(Cluster))
  
  if (nrow(df_plot) == 0) {
    stop("Nenhum dado válido após filtragem para o gráfico normalizado.")
  }
  
  message("Preparando dados para normalização: ", nrow(df_plot), " registros, ",
          length(unique(df_plot$Cluster)), " clusters, ",
          length(unique(df_plot$Frame)), " frames")
  
  # Normalizar valores por frame (percentual)
  df_normalized <- df_plot %>%
    group_by(Frame) %>%
    mutate(
      Total_Per_Frame = sum(Value, na.rm = TRUE),
      Normalized = ifelse(Total_Per_Frame == 0, 0, 100 * Value / Total_Per_Frame)
    ) %>%
    ungroup()
  
  # Verificar se a normalização funcionou
  normalization_summary <- df_normalized %>%
    group_by(Frame) %>%
    summarise(
      Total_Normalized = sum(Normalized, na.rm = TRUE),
      .groups = "drop"
    )
  
  message("Verificação de normalização - Soma por frame deve ser ~100: ",
          paste(round(unique(normalization_summary$Total_Normalized), 1), collapse = ", "))
  
  # Preparar paleta de cores consistente
  n_clusters <- length(unique(df_normalized$Cluster))
  cluster_colors <- brewer.pal(min(9, max(3, n_clusters)), "Set1")
  if (n_clusters > length(cluster_colors)) {
    cluster_colors <- colorRampPalette(cluster_colors)(n_clusters)
  }
  
  # Criar gráfico ggplot
  p <- ggplot(df_normalized, aes(x = Frame, y = Normalized, color = Cluster, group = Cluster)) +
    geom_line(size = 1.2, alpha = 0.8) +
    geom_point(size = 1.5, alpha = 0.7) +
    scale_color_manual(values = cluster_colors) +
    labs(
      title = "Evolução Temporal Normalizada da Importância dos Clusters",
      subtitle = "Valores normalizados por frame (soma = 100%) para comparação relativa",
      x = "Valor do PC1 (Progressão Temporal)",
      y = "Importância Relativa (%)",
      color = "Cluster"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.position = "right",
      panel.grid.major = element_line(color = "grey90", size = 0.2),
      panel.grid.minor = element_blank()
    ) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))
  
  return(p)
}

# Função alternativa para normalização por cluster (diferente abordagem)
criar_grafico_normalizado_por_cluster <- function(df) {
  
  # Verificar colunas disponíveis
  available_cols <- colnames(df)
  message("Criando gráfico normalizado por cluster - colunas: ", paste(available_cols, collapse = ", "))
  
  # Usar Importance_smooth_norm se disponível, caso contrário calcular
  if ("Importance_smooth_norm" %in% available_cols) {
    message("Usando Importance_smooth_norm pré-calculada")
    df_plot <- df %>%
      filter(!is.na(Importance_smooth_norm)) %>%
      mutate(
        Normalized = Importance_smooth_norm * 100,  # Converter para percentual
        Frame = PC1_value,
        Cluster = as.factor(Cluster)
      )
  } else {
    message("Calculando normalização por cluster...")
    df_plot <- df %>%
      mutate(
        Value = Importance_smooth,
        Frame = PC1_value,
        Cluster = as.factor(Cluster)
      ) %>%
      filter(!is.na(Value), !is.na(Frame)) %>%
      group_by(Cluster) %>%
      mutate(
        Max_Value = max(Value, na.rm = TRUE),
        Normalized = ifelse(Max_Value == 0, 0, 100 * Value / Max_Value)
      ) %>%
      ungroup()
  }
  
  if (nrow(df_plot) == 0) {
    stop("Nenhum dado válido para gráfico normalizado por cluster.")
  }
  
  # Preparar paleta de cores
  n_clusters <- length(unique(df_plot$Cluster))
  cluster_colors <- brewer.pal(min(9, max(3, n_clusters)), "Set1")
  if (n_clusters > length(cluster_colors)) {
    cluster_colors <- colorRampPalette(cluster_colors)(n_clusters)
  }
  
  p <- ggplot(df_plot, aes(x = Frame, y = Normalized, color = Cluster, group = Cluster)) +
    geom_line(size = 1.2, alpha = 0.8) +
    geom_point(size = 1.5, alpha = 0.6) +
    scale_color_manual(values = cluster_colors) +
    labs(
      title = "Importância Relativa dos Clusters (Normalizada por Cluster)",
      subtitle = "Cada cluster normalizado pelo seu próprio valor máximo (0-100%)",
      x = "Valor do PC1 (Progressão Temporal)",
      y = "Importância Relativa (%)",
      color = "Cluster"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right"
    ) +
    ylim(0, 100)
  
  return(p)
}

############################################
### 7. EXECUÇÃO DOS GRÁFICOS ADICIONAIS ###
############################################

message("\n=== CRIAÇÃO DE GRÁFICOS ADICIONAIS ===")

# Verificar estrutura do temporal_smoothed
if (!exists("temporal_smoothed")) {
  stop("O objeto 'temporal_smoothed' não foi encontrado.")
}

message("Estrutura do temporal_smoothed:")
print(str(temporal_smoothed))
message("\nPrimeiras linhas do temporal_smoothed:")
print(head(temporal_smoothed))

# Tentar diferentes abordagens de normalização
tryCatch({
  message("\nTentativa 1: Normalização por frame...")
  normalized_plot_1 <- criar_grafico_normalizado(temporal_smoothed)
  ggsave(file.path(output_dir, "images", "importancia_normalizada_por_frame.png"),
         normalized_plot_1, width = 14, height = 8, dpi = 600)
  message("✅ Gráfico normalizado por frame salvo")
}, error = function(e) {
  warning("Falha na normalização por frame: ", e$message)
})

tryCatch({
  message("\nTentativa 2: Normalização por cluster...")
  normalized_plot_2 <- criar_grafico_normalizado_por_cluster(temporal_smoothed)
  ggsave(file.path(output_dir, "images", "importancia_normalizada_por_cluster.png"),
         normalized_plot_2, width = 14, height = 8, dpi = 600)
  message("✅ Gráfico normalizado por cluster salvo")
}, error = function(e) {
  warning("Falha na normalização por cluster: ", e$message)
})

# Criar heatmap temporal (função existente)
tryCatch({
  message("\nCriando heatmap temporal...")
  heatmap_plot <- criar_heatmap_temporal(temporal_smoothed)
  ggsave(file.path(output_dir, "images", "heatmap_temporal_clusters.png"),
         heatmap_plot, width = 12, height = 8, dpi = 600)
  message("✅ Heatmap temporal salvo")
}, error = function(e) {
  warning("Falha ao criar heatmap: ", e$message)
})

############################################
### 8. GRÁFICO DE ÁREA STACKED ###
############################################

# Gráfico adicional: área stacked para ver contribuição percentual
criar_grafico_area_stacked <- function(df) {
  
  df_area <- df %>%
    mutate(
      Value = Importance_smooth,
      Frame = PC1_value,
      Cluster = as.factor(Cluster)
    ) %>%
    filter(!is.na(Value), !is.na(Frame)) %>%
    group_by(Frame) %>%
    mutate(
      Total_Per_Frame = sum(Value, na.rm = TRUE),
      Percentual = ifelse(Total_Per_Frame == 0, 0, 100 * Value / Total_Per_Frame)
    ) %>%
    ungroup()
  
  # Preparar paleta de cores
  n_clusters <- length(unique(df_area$Cluster))
  cluster_colors <- brewer.pal(min(9, max(3, n_clusters)), "Set1")
  if (n_clusters > length(cluster_colors)) {
    cluster_colors <- colorRampPalette(cluster_colors)(n_clusters)
  }
  
  p_area <- ggplot(df_area, aes(x = Frame, y = Percentual, fill = Cluster)) +
    geom_area(alpha = 0.8, color = "white", size = 0.2) +
    scale_fill_manual(values = cluster_colors) +
    labs(
      title = "Contribuição Percentual dos Clusters ao Longo do Tempo",
      subtitle = "Área stacked mostrando a contribuição relativa de cada cluster",
      x = "Valor do PC1 (Progressão Temporal)",
      y = "Contribuição Percentual (%)",
      fill = "Cluster"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right"
    ) +
    ylim(0, 100)
  
  return(p_area)
}

tryCatch({
  message("\nCriando gráfico de área stacked...")
  area_plot <- criar_grafico_area_stacked(temporal_smoothed)
  ggsave(file.path(output_dir, "images", "area_stacked_clusters.png"),
         area_plot, width = 14, height = 8, dpi = 600)
  message("✅ Gráfico de área stacked salvo")
}, error = function(e) {
  warning("Falha ao criar gráfico de área: ", e$message)
})

message("\n✅ Todos os gráficos adicionais processados")

############################################
### 9. VERIFICAÇÃO E SALVAMENTO FINAL ###
############################################

message("\n=== VERIFICAÇÃO FINAL DOS OBJETOS ===")

# Lista de objetos esperados
objetos_esperados <- c("gene_clusters", "temporal_smoothed", "statistical_analysis")
objetos_presentes <- ls()

# Verificar quais objetos estão presentes
objetos_faltantes <- setdiff(objetos_esperados, objetos_presentes)
objetos_disponiveis <- intersect(objetos_esperados, objetos_presentes)

message("Objetos disponíveis para salvamento: ", paste(objetos_disponiveis, collapse = ", "))

if (length(objetos_faltantes) > 0) {
  message("Objetos faltantes: ", paste(objetos_faltantes, collapse = ", "))
  
  # Recriar objetos faltantes se possível
  if ("statistical_analysis" %in% objetos_faltantes) {
    message("Recriando 'statistical_analysis'...")
    
    # Função robusta para análise estatística
    realizar_analise_estatistica <- function(temporal_df, clusters_df) {
      
      # Validar dados de entrada
      if (nrow(temporal_df) == 0 || nrow(clusters_df) == 0) {
        stop("Dados vazios para análise estatística.")
      }
      
      message("Realizando análise estatística...")
      
      # Correlação entre clusters
      correlation_data <- temporal_df %>%
        select(Cluster, PC1_value, Importance_smooth) %>%
        pivot_wider(names_from = Cluster, values_from = Importance_smooth) %>%
        select(-PC1_value)
      
      # Verificar se há dados suficientes para correlação
      if (ncol(correlation_data) >= 2 && nrow(correlation_data) >= 3) {
        correlation_matrix <- cor(correlation_data, use = "complete.obs")
      } else {
        warning("Dados insuficientes para cálculo de correlação.")
        correlation_matrix <- matrix(NA, 
                                     nrow = length(unique(temporal_df$Cluster)), 
                                     ncol = length(unique(temporal_df$Cluster)))
        colnames(correlation_matrix) <- rownames(correlation_matrix) <- unique(temporal_df$Cluster)
      }
      
      # Estatísticas descritivas
      descriptive_stats <- temporal_df %>%
        group_by(Cluster) %>%
        summarise(
          Mean_Importance = mean(Importance_smooth, na.rm = TRUE),
          SD_Importance = sd(Importance_smooth, na.rm = TRUE),
          Max_Importance = max(Importance_smooth, na.rm = TRUE),
          Time_of_Max = PC1_value[which.max(Importance_smooth)],
          CV = ifelse(Mean_Importance > 0, SD_Importance / Mean_Importance, NA),
          .groups = "drop"
        )
      
      # Informações dos clusters
      cluster_info <- clusters_df %>%
        group_by(Cluster) %>%
        summarise(
          n_genes = n(),
          mean_position = mean(Position, na.rm = TRUE),
          sd_position = sd(Position, na.rm = TRUE),
          .groups = "drop"
        )
      
      # Adicionar VariationSmooth se disponível
      if ("VariationSmooth" %in% colnames(clusters_df)) {
        variation_info <- clusters_df %>%
          group_by(Cluster) %>%
          summarise(
            mean_variation = mean(VariationSmooth, na.rm = TRUE),
            .groups = "drop"
          )
        cluster_info <- cluster_info %>% left_join(variation_info, by = "Cluster")
      }
      
      # Combinar estatísticas
      final_stats <- descriptive_stats %>%
        left_join(cluster_info, by = "Cluster") %>%
        arrange(Time_of_Max)
      
      return(list(
        correlation_matrix = correlation_matrix,
        cluster_statistics = final_stats
      ))
    }
    
    # Executar análise estatística
    tryCatch({
      statistical_analysis <- realizar_analise_estatistica(temporal_smoothed, gene_clusters)
      message("✅ 'statistical_analysis' criado com sucesso")
    }, error = function(e) {
      warning("Falha ao criar statistical_analysis: ", e$message)
      # Criar objeto vazio como fallback
      statistical_analysis <- list(
        correlation_matrix = matrix(NA, nrow = 1, ncol = 1),
        cluster_statistics = data.frame(Cluster = "N/A", Message = "Análise não disponível")
      )
    })
  }
  
  # Verificar se gene_clusters está faltante
  if ("gene_clusters" %in% objetos_faltantes) {
    stop("Objeto 'gene_clusters' é essencial e não foi encontrado. Execute a identificação de clusters primeiro.")
  }
  
  # Verificar se temporal_smoothed está faltante
  if ("temporal_smoothed" %in% objetos_faltantes) {
    stop("Objeto 'temporal_smoothed' é essencial e não foi encontrado. Execute a preparação dos dados temporais primeiro.")
  }
}

# Verificação final antes do salvamento
objetos_para_salvar <- intersect(objetos_esperados, ls())
message("Objetos que serão salvos: ", paste(objetos_para_salvar, collapse = ", "))

if (length(objetos_para_salvar) == 0) {
  stop("Nenhum objeto disponível para salvamento. Verifique a execução do script.")
}

############################################
### 10. SALVAMENTO DOS DADOS FINAIS #######
############################################

message("\n=== SALVAMENTO DOS DADOS FINAIS ===")

# 1. Salvar dados temporais processados
tryCatch({
  write.csv(temporal_smoothed,
            file.path(output_dir, "dados_temporais_suavizados.csv"),
            row.names = FALSE)
  message("✅ Dados temporais suavizados salvos: dados_temporais_suavizados.csv")
}, error = function(e) {
  warning("Falha ao salvar dados temporais: ", e$message)
})

# 2. Salvar clusters de genes
tryCatch({
  write.csv(gene_clusters,
            file.path(output_dir, "clusters_genes_identificados.csv"),
            row.names = FALSE)
  message("✅ Clusters de genes salvos: clusters_genes_identificados.csv")
}, error = function(e) {
  warning("Falha ao salvar clusters de genes: ", e$message)
})

# 3. Salvar análise estatística se disponível
if (exists("statistical_analysis")) {
  tryCatch({
    # Salvar matriz de correlação
    if (!is.null(statistical_analysis$correlation_matrix)) {
      write.csv(statistical_analysis$correlation_matrix,
                file.path(output_dir, "matriz_correlacao_clusters.csv"),
                row.names = TRUE)
      message("✅ Matriz de correlação salva: matriz_correlacao_clusters.csv")
    }
    
    # Salvar estatísticas dos clusters
    if (!is.null(statistical_analysis$cluster_statistics)) {
      write.csv(statistical_analysis$cluster_statistics,
                file.path(output_dir, "estatisticas_descritivas_clusters.csv"),
                row.names = FALSE)
      message("✅ Estatísticas descritivas salvas: estatisticas_descritivas_clusters.csv")
    }
  }, error = function(e) {
    warning("Falha ao salvar análise estatística: ", e$message)
  })
}

# 4. Salvar objeto R com todos os resultados
tryCatch({
  # Criar lista com todos os objetos disponíveis
  objetos_salvar <- list()
  
  if (exists("gene_clusters")) {
    objetos_salvar$gene_clusters <- gene_clusters
  }
  if (exists("temporal_smoothed")) {
    objetos_salvar$temporal_smoothed <- temporal_smoothed
  }
  if (exists("statistical_analysis")) {
    objetos_salvar$statistical_analysis <- statistical_analysis
  }
  
  saveRDS(objetos_salvar, file.path(output_dir, "analise_temporal_completa.rds"))
  message("✅ Objeto R salvo: analise_temporal_completa.rds")
  
  # Também salvar como .RData para compatibilidade
  save(list = objetos_para_salvar, 
       file = file.path(output_dir, "analise_temporal_completa.RData"))
  message("✅ Arquivo RData salvo: analise_temporal_completa.RData")
  
}, error = function(e) {
  warning("Falha ao salvar objetos R: ", e$message)
})

# 5. Salvar metadados e informações da sessão
tryCatch({
  # Informações da sessão
  sink(file.path(output_dir, "session_info.txt"))
  print(sessionInfo())
  sink()
  
  # Metadados da análise
  metadata <- data.frame(
    Data_Execucao = as.character(Sys.time()),
    Numero_Clusters = if(exists("gene_clusters")) length(unique(gene_clusters$Cluster)) else NA,
    Numero_Genes_Clusterizados = if(exists("gene_clusters")) nrow(gene_clusters) else NA,
    Numero_Frames_Temporais = if(exists("temporal_smoothed")) length(unique(temporal_smoothed$PC1_value)) else NA,
    Script_Utilizado = "analise_temporal_clusters_corrigido.R"
  )
  
  write.csv(metadata, file.path(output_dir, "metadata_analise.csv"), row.names = FALSE)
  message("✅ Metadados salvos: metadata_analise.csv")
  
}, error = function(e) {
  warning("Falha ao salvar metadados: ", e$message)
})

############################################
### 11. RELATÓRIO FINAL DE EXECUÇÃO #######
############################################

message(paste0("\n", strrep("=", 60)))
message("🎯 RELATÓRIO FINAL DA ANÁLISE TEMPORAL DE CLUSTERS")
message(strrep("=", 60))

# Resumo estatístico
if (exists("gene_clusters") && exists("temporal_smoothed")) {
  cat("\n📊 RESUMO ESTATÍSTICO:\n")
  cat("• Clusters identificados:", length(unique(gene_clusters$Cluster)), "\n")
  cat("• Genes clusterizados:", nrow(gene_clusters), "\n")
  cat("• Frames temporais analisados:", length(unique(temporal_smoothed$PC1_value)), "\n")
  cat("• Período do PC1: [", min(temporal_smoothed$PC1_value, na.rm = TRUE), 
      ", ", max(temporal_smoothed$PC1_value, na.rm = TRUE), "]\n", sep = "")
}

# Arquivos gerados
if (dir.exists(output_dir)) {
  cat("\n📁 ARQUIVOS GERADOS:\n")
  arquivos_gerados <- list.files(output_dir, recursive = TRUE, full.names = FALSE)
  categorias <- list(
    "Dados" = arquivos_gerados[grepl("\\.(csv|rds|RData)$", arquivos_gerados)],
    "Imagens" = arquivos_gerados[grepl("\\.(png|jpg|pdf)$", arquivos_gerados)],
    "Relatórios" = arquivos_gerados[grepl("\\.(txt|log)$", arquivos_gerados)]
  )
  
  for (categoria in names(categorias)) {
    if (length(categorias[[categoria]]) > 0) {
      cat("  ", categoria, " (", length(categorias[[categoria]]), "):\n", sep = "")
      for (arquivo in head(categorias[[categoria]], 10)) { # Mostrar apenas os 10 primeiros
        cat("    - ", arquivo, "\n", sep = "")
      }
      if (length(categorias[[categoria]]) > 10) {
        cat("    ... e mais", length(categorias[[categoria]]) - 10, "arquivos\n")
      }
    }
  }
  cat("• Total de arquivos gerados:", length(arquivos_gerados), "\n")
}

# Próximos passos
cat("\n🚀 PRÓXIMOS PASSOS SUGERIDOS:\n")
cat("1. Verificar gráficos em: ", file.path(output_dir, "images"), "\n")
cat("2. Analisar correlações entre clusters\n")
cat("3. Realizar enriquecimento funcional por cluster\n")
cat("4. Validar clusters com anotações biológicas conhecidas\n")
cat("5. Explorar padrões temporais em vias de sinalização específicas\n")

# Status final
cat("\n✅ ANÁLISE CONCLUÍDA COM SUCESSO!\n")
cat("📍 Diretório de resultados: ", normalizePath(output_dir), "\n")
cat("⏰ Timestamp: ", as.character(Sys.time()), "\n")

message(strrep("=", 60))
