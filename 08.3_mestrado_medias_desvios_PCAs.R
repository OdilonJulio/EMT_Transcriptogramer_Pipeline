# Carregar pacotes
library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)
library(magrittr)
library(transcriptogramer)
library(biomaRt)
library(dplyr)
library(vroom)
library(RcppCWB)
library(Ropj)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(RColorBrewer)

# Carregar os dataframes
load("t_matrix_batch_effect_correction_R0.RData")
load("t_matrix_batch_effect_correction_R30.RData")

##############
#### Funções de PCA
##############

# Função para preparar os dados para PCA
prepare_pca_data <- function(df) {
  if (any(is.na(df))) {
    stop("O dataset contém valores NA. Verifique e limpe antes de prosseguir.")
  }
  df <- df[, !is.na(colnames(df))]  # Remover colunas com rótulos NA
  rownames(df) <- df$Protein        # Usar a coluna "Protein" como nomes das linhas
  df <- df[, -c(1, 2)] # Remover as colunas "Protein" e "Position"
  return(t(df))                     # Transpor a matriz para PCA
}

# Função para rodar PCA
run_pca <- function(pca_input) {
  # Remover colunas com variância zero
  variances <- apply(pca_input, 2, var)
  pca_input_cleaned <- pca_input[, variances > 0]
  
  if (ncol(pca_input_cleaned) == 0) {
    stop("Todas as colunas têm variância zero.")
  }
  
  # Rodar PCA
  pca_result <- prcomp(pca_input_cleaned, center = TRUE, scale. = FALSE)
  
  # Calcular autovalores e variâncias explicadas
  eigenvalues <- pca_result$sdev^2
  explained_variance <- eigenvalues / sum(eigenvalues)
  cumulative_variance <- cumsum(explained_variance)
  
  cat("Eigenvalues (autovalores):\n", eigenvalues, "\n")
  cat("Explained Variance:\n", explained_variance, "\n")
  cat("Sum of Explained Variance:\n", sum(explained_variance), "\n")
  cat("Cumulative Variance:\n", cumulative_variance, "\n")
  
  return(list(
    pca_result = pca_result,
    eigenvalues = eigenvalues,
    explained_variance = explained_variance,
    cumulative_variance = cumulative_variance
  ))
}

# Função para extrair as condições dos nomes das colunas
extract_condition <- function(colnames, conditions) {
  sapply(colnames, function(col) {
    match <- grep(paste(conditions, collapse = "|"), col, value = TRUE)
    if (length(match) > 0) {
      return(sub("_.*", "", match[1]))  # Extrai apenas o prefixo da condição
    } else {
      return("Unknown")  # Caso nenhuma condição seja encontrada
    }
  })
}



##############
#### Processo de PCA
##############

# Lista de dias para filtrar
Day <- c("notreated-batch1", "notreated-batch2", "TGFbeta1-1day-batch2", 
                "TGFbeta1-2day-batch2", "TGFbeta1-3day-batch2", 
                "TGFbeta1-4day-batch1", "TGFbeta1-8day-batch1")

# Preparar os dados para PCA
df_t_matrix_batch_effect_correction_R0_for_PCA <- prepare_pca_data(t_matrix_batch_effect_correction_R0@transcriptogramS2)
df_t_matrix_batch_effect_correction_R30_for_PCA <- prepare_pca_data(t_matrix_batch_effect_correction_R30@transcriptogramS2)

# Gerar labels de condições
condition_labels_R0 <- extract_condition(rownames(df_t_matrix_batch_effect_correction_R0_for_PCA), Day)
condition_labels_R30 <- extract_condition(rownames(df_t_matrix_batch_effect_correction_R30_for_PCA), Day)

# Garantir que os labels sejam fatores para um mapeamento correto de cores
condition_labels_R0 <- factor(condition_labels_R0, levels = Day)
condition_labels_R30 <- factor(condition_labels_R30, levels = Day)

# Definir paleta de cores fixa para 7 condições
color_palette <- c(
  "notreated-batch1" = "#984EA3",  # Roxo
  "notreated-batch2" = "#377EB8",  # Azul
  "TGFbeta1-1day-batch2" = "#4DAF4A",  # Verde
  "TGFbeta1-2day-batch2" = "#3A5400", # Verde Escuro
  "TGFbeta1-3day-batch2" = "#A65628",   # Marrom 
  "TGFbeta1-4day-batch1" = "#FF7F00",  # Laranja
  "TGFbeta1-8day-batch1" = "#E41A1C"  # Vermelho
)

names(color_palette) <- Day  # Mapear cores para condições

# Rodar PCA
pca_result_R0 <- run_pca(df_t_matrix_batch_effect_correction_R0_for_PCA)
pca_result_R30 <- run_pca(df_t_matrix_batch_effect_correction_R30_for_PCA)

add_day_column <- function(pca_df, rownames_data) {
  # Criar a coluna Day a partir dos nomes das linhas (barcodes)
  pca_df$Day <- case_when(
    grepl("notreated-batch1", rownames_data) ~ "notreated-batch1",
    grepl("notreated-batch2", rownames_data) ~ "notreated-batch2",
    grepl("TGFbeta1-1day-batch2", rownames_data) ~ "TGFbeta1-1day-batch2",
    grepl("TGFbeta1-2day-batch2", rownames_data) ~ "TGFbeta1-2day-batch2",
    grepl("TGFbeta1-3day-batch2", rownames_data) ~ "TGFbeta1-3day-batch2",
    grepl("TGFbeta1-4day-batch1", rownames_data) ~ "TGFbeta1-4day-batch1",
    grepl("TGFbeta1-8day-batch1", rownames_data) ~ "TGFbeta1-8day-batch1",
    TRUE ~ "Unknown"
  )
  
  # Verificar quantas amostras não foram mapeadas
  unmatched <- sum(pca_df$Day == "Unknown")
  if (unmatched > 0) {
    warning(unmatched, " amostras não foram mapeadas para nenhuma condição conhecida.")
  }
  
  return(pca_df)
}

# Aplicar a função corretamente usando rownames(pca_r0)
pca_r0 <- as.data.frame(pca_result_R0$pca_result$x)
pca_r0 <- add_day_column(pca_r0, rownames(pca_r0))  # Usando rownames(pca_r0)

pca_r30 <- as.data.frame(pca_result_R30$pca_result$x)
pca_r30 <- add_day_column(pca_r30, rownames(pca_r30))

### CENTRO DE MASSA DAS PCAS

# Função para calcular o centro de massa
calculate_centroid <- function(pca_df) {
  # Calcular a média das coordenadas das PCs para cada condição
  centroid <- aggregate(. ~ Day, data = pca_df, FUN = mean)
  
  # Retornar o centro de massa calculado
  return(centroid)
}

# Calcular o centro de massa para pca_r0 e pca_r30
centroid_r0 <- calculate_centroid(pca_r0)
centroid_r30 <- calculate_centroid(pca_r30)


# Gráfico de PCA para pca_r0 com centros de massa
pca_plot_r0 <- ggplot(pca_r0, aes(x = PC1, y = PC2, color = Day)) +
  geom_point(size = 1) +  # Pontos das amostras
  geom_point(data = centroid_r0, aes(x = PC1, y = PC2), 
             color = "black", size = 4, shape = 3) +  # Centros de massa
  scale_color_manual(values = color_palette) +  # Cores personalizadas
  theme_minimal() +
  labs(title = "PCA com Centros de Massa (R30)", 
       x = "PC1", 
       y = "PC2")

# Gráfico de PCA para pca_r30 com centros de massa
pca_plot_r30 <- ggplot(pca_r30, aes(x = PC1, y = PC2, color = Day)) +
  geom_point(size = 1) +  # Pontos das amostras
  geom_point(data = centroid_r30, aes(x = PC1, y = PC2), 
             color = "black", size = 4, shape = 3) +  # Centros de massa
  scale_color_manual(values = color_palette) +  # Cores personalizadas
  theme_minimal() +
  labs(title = "PCA com Centros de Massa (R0)", 
       x = "PC1", 
       y = "PC2")

# Exibir os gráficos
print(pca_plot_r0)
print(pca_plot_r30)

# Salvar gráfico pca_plot_r0 em alta qualidade
ggsave("images/center_of_mass_pca_plot_r0.png", plot = pca_plot_r0, 
       dpi = 300, width = 8, height = 6, units = "in")

# Salvar gráfico pca_plot_r30 em alta qualidade
ggsave("images/center_of_mass_pca_plot_r30.png", plot = pca_plot_r30, 
       dpi = 300, width = 8, height = 6, units = "in")



########


calculate_mean_sd_by_day <- function(pca_df) {
  pca_df %>%
    select(Day, starts_with("PC")) %>%  # Seleciona apenas as colunas 'Day' e as PCs
    group_by(Day) %>%  # Agrupa por 'Day'
    summarise(
      across(starts_with("PC"), 
             list(mean = ~mean(., na.rm = TRUE), sd = ~sd(., na.rm = TRUE)),
             .names = "{.col}_{.fn}"),  # Renomeia as colunas como 'PC1_mean', 'PC1_sd', etc.
      .groups = 'drop'
    ) %>%
    mutate(across(contains("_mean") | contains("_sd"), ~replace_na(., 0))) %>%  # Substitui NA por 0
    arrange(factor(Day, levels = c("notreated-batch1",
                                   "notreated-batch2", 
                                   "TGFbeta1-1day-batch2",
                                   "TGFbeta1-2day-batch2", 
                                   "TGFbeta1-3day-batch2", 
                                   "TGFbeta1-4day-batch1",
                                   "TGFbeta1-8day-batch1")))  # Organiza os dias
}




mean_sd_by_day_r0 <- calculate_mean_sd_by_day(pca_r0)
mean_sd_by_day_r30 <- calculate_mean_sd_by_day(pca_r30)

# Exportar médias
write.table(mean_sd_by_day_r0, "mean_sd_by_day_r0.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(mean_sd_by_day_r30, "mean_sd_by_day_r30.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

##############
#### Geração de Gráficos
##############

library(ggplot2)
library(patchwork)


pca_r0$Day <- factor(pca_r0$Day, levels = names(color_palette))
pca_r30$Day <- factor(pca_r30$Day, levels = names(color_palette))




generate_pca_scatter_plots_with_lines <- function(pca_df, line_data, color_palette) {
  # Verificar correspondência de 'Day' com os nomes da paleta
  if (!all(unique(pca_df$Day) %in% names(color_palette))) {
    stop("Os valores de 'Day' no pca_df não correspondem aos nomes em color_palette.")
  }
  
  # Certificar-se de que 'Day' seja um fator ordenado
  pca_df$Day <- factor(pca_df$Day, levels = unique(line_data$Day))
  line_data$Day <- factor(line_data$Day, levels = unique(line_data$Day))
  
  # Calcular os limites globais dos eixos para PC1 e todos os outros PCs
  global_x_limits <- range(c(pca_df$PC1, line_data$PC1_mean), na.rm = TRUE)
  
  # Encontrar o maior limite de todos os PCs para ter eixos com o mesmo tamanho
  global_y_limits <- range(c(pca_df[, grep("^PC", colnames(pca_df))], 
                             line_data[, grep("^PC", colnames(line_data))]), na.rm = TRUE)
  
  # Gerar os gráficos para PC2 a PC10
  plots <- lapply(2:10, function(i) {
    # Preparar dados para linha
    line_data_current <- line_data %>%
      select(Day, PC1_mean, PC_mean = paste0("PC", i, "_mean"))
    
    # Criar o gráfico
    p <- ggplot(pca_df, aes(x = PC1, y = .data[[paste0("PC", i)]], color = Day)) +
      geom_point(alpha = 0.6, size = 0.4) + # Pontos de dispersão menores
      geom_path(data = line_data_current, aes(x = PC1_mean, y = PC_mean, group = 1), 
                color = "gray", linewidth = 0.5, arrow = arrow(length = unit(0.2, "cm"))) + # Linha direcional
      geom_point(data = line_data_current, aes(x = PC1_mean, y = PC_mean), 
                 color = "black", size = 2, shape = 16) + # Pontos maiores para as médias
      scale_color_manual(values = color_palette) + # Aplicar paleta de cores
      scale_x_continuous(limits = global_x_limits) +
      scale_y_continuous(limits = global_y_limits) + # Ajustar Y para ter os mesmos limites
      labs(title = paste("PC1 vs", paste0("PC", i)),
           x = "PC1", y = paste0("PC", i),
           color = "Condition") + # Título da legenda
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)
      )
    
    # Remover legenda dos gráficos subsequentes
    if (i != 2) {
      p <- p + theme(legend.position = "none", legend.title = element_blank())
    }
    
    return(p)
  })
  
  # Combinar os gráficos
  combined_plot <- wrap_plots(plots, ncol = 3) + 
    plot_annotation(title = "PCA Scatter Plots with Directed Lines") &
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    )
  
  return(combined_plot)
}



# Gerar os gráficos
pca_scatter_with_lines_r0 <- generate_pca_scatter_plots_with_lines(pca_r0, mean_sd_by_day_r0, color_palette)
pca_scatter_with_lines_r30 <- generate_pca_scatter_plots_with_lines(pca_r30, mean_sd_by_day_r30, color_palette)

# Salvar os gráficos
ggsave("images/PCA_Scatter_R0_with_Lines.png", pca_scatter_with_lines_r0, width = 15, height = 15, dpi = 300)
ggsave("images/PCA_Scatter_R30_with_Lines.png", pca_scatter_with_lines_r30, width = 15, height = 15, dpi = 300)
