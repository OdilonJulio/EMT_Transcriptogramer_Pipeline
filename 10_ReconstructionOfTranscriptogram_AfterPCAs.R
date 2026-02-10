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

# Carregar os dataframes DOS TRANSCRIPTOGRAMAS
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


################## RECONSTRUÇÃO DO TRANSCRIPTOGRAMA ##########################

# Função para expandir a matriz de centro de massa com base na contagem de ocorrências
expand_center_of_mass <- function(pca_data, centroid_data) {
  # Contar quantas vezes cada valor de centroid_data$Day aparece em pca_data$Day
  counts <- table(pca_data$Day)
  
  # Garantir que todos os valores de centroid_data$Day estão presentes em counts
  counts <- counts[names(counts) %in% centroid_data$Day]
  
  # Repetir as linhas de centroid_data de acordo com as contagens em counts
  expanded_centroid <- centroid_data[rep(seq_len(nrow(centroid_data)), counts), , drop = FALSE]
  
  # Transpor a matriz para adequação à análise
  full_centroid <- t(expanded_centroid)
  
  # Remover a linha correspondente ao identificador "Day" para deixar apenas os valores numéricos
  center_of_mass <- full_centroid[-1, , drop = FALSE]
  
  return(center_of_mass)
}

# Expandir centros de massa para R0 e R30
center_of_mass_r0 <- expand_center_of_mass(pca_r0, centroid_r0)
center_of_mass_r0 <- apply(center_of_mass_r0, 2, as.numeric)

center_of_mass_r30 <- expand_center_of_mass(pca_r30, centroid_r30)
center_of_mass_r30 <- apply(center_of_mass_r30, 2, as.numeric)


# Verificar dimensões das matrizes processadas
cat("Dimensões de center_of_mass_r0:", dim(center_of_mass_r0), "\n")
cat("Dimensões de center_of_mass_r30:", dim(center_of_mass_r30), "\n")

reconstruct_transcriptogram <- function(pca_result, center_of_mass) {
  # Matriz de rotação (loadings das PCAs)
  rotation <- pca_result$pca_result$rotation
  print(dim(rotation))  # Deve ser (n_genes, n_PC)
  
  # Componentes principais
  components <- as.matrix(pca_result$pca_result$x)
  print(dim(components))  # Deve ser (n_células, n_PC)
  
  # Multiplicação entre 'rotation' e a transposta de 'components'
  reconstructed_transcriptogram <- rotation %*% t(components)
  print(dim(reconstructed_transcriptogram))  # Deve ser (n_genes, n_células)
  
  # Garantir que center_of_mass tem a dimensão correta
  if (!all(dim(center_of_mass) == dim(reconstructed_transcriptogram))) {
    stop("Erro: center_of_mass tem dimensão incompatível com o transcriptograma reconstruído.")
  }
  
  # Adicionar o centro de massa
  reconstructed_transcriptogram <- reconstructed_transcriptogram + center_of_mass
  reconstructed_transcriptogram <- t(reconstructed_transcriptogram)
  print(dim(reconstructed_transcriptogram))
  # Ajustar colnames se necessário
  # if (!is.null(colnames(components))) {
  #   if (ncol(reconstructed_transcriptogram) == length(colnames(components))) {
  #     colnames(reconstructed_transcriptogram) <- colnames(components)
  #   } else {
  #     warning("Número de colunas do transcriptograma reconstruído não bate com os colnames de 'components'.")
  #   }
  # }
  
  print("Reconstrução finalizada com sucesso.")
  return(reconstructed_transcriptogram)
}

# Reconstruindo os transcriptogramas
reconstructed_r0 <- reconstruct_transcriptogram(pca_result_R0, center_of_mass_r0)
reconstructed_r30 <- reconstruct_transcriptogram(pca_result_R30, center_of_mass_r30)

# Verificar dimensões finais das matrizes reconstruídas
cat("Dimensões de reconstructed_r0:", dim(reconstructed_r0), "\n")
cat("Dimensões de reconstructed_r30:", dim(reconstructed_r30), "\n")

# Suponha que você tenha dois dataframes: df1 e df2

# Verificar se os nomes das linhas (row names) são idênticos
identical(rownames(df_t_matrix_batch_effect_correction_R30_for_PCA), rownames(reconstructed_r30))  # Retorna TRUE se forem iguais, caso contrário FALSE

# Verificar se os nomes das colunas (column names) são idênticos
identical(colnames(df_t_matrix_batch_effect_correction_R30_for_PCA), colnames(reconstructed_r30))  # Retorna TRUE se forem iguais, caso contrário FALSE


length(setdiff(colnames(df_t_matrix_batch_effect_correction_R30_for_PCA), colnames(reconstructed_r30)))


####################### IMPLEMENTAÇÃO DE R² ############################

calculate_r2_cumulative <- function(original, pca_result, plot = TRUE) {
  # Converter `original` para matriz numérica
  original <- as.matrix(original)
  
  # Verificar se há valores NA ou infinitos
  if (any(is.na(original))) {
    stop("Erro: O dataset original contém valores NA.")
  }
  if (any(is.infinite(original))) {
    stop("Erro: O dataset original contém valores infinitos.")
  }
  
  # Verificar se há nomes de linhas e colunas
  if (is.null(rownames(original)) || is.null(colnames(original))) {
    stop("Erro: O dataframe original precisa ter nomes de genes (rownames) e células (colnames).")
  }
  
  # Obter nomes de genes e células do original
  common_genes <- rownames(original)
  common_cells <- colnames(original)
  
  # Obter loadings (rotation) e componentes principais (PCs)
  rotation <- pca_result$pca_result$rotation
  components <- as.matrix(pca_result$pca_result$x)
  
  # Verificar se `rotation` e `components` têm nomes para fazer o alinhamento correto
  if (is.null(rownames(rotation))) {
    stop("Erro: O objeto PCA não contém nomes de genes em `rotation`.")
  }
  if (is.null(rownames(components))) {
    stop("Erro: O objeto PCA não contém nomes de células em `components`.")
  }
  
  # Ajustar dimensões: manter apenas genes e células em comum
  common_genes <- intersect(common_genes, rownames(rotation))
  common_cells <- intersect(common_cells, rownames(components))
  
  # Mensagem de depuração
  cat("Número de genes em comum:", length(common_genes), "\n")
  cat("Número de células em comum:", length(common_cells), "\n")
  
  # Se não houver genes ou células em comum, interrompe a execução
  if (length(common_genes) == 0 || length(common_cells) == 0) {
    stop("Erro: Nenhum gene ou célula em comum entre `original` e `pca_result`.")
  }
  
  # Filtrar apenas os genes e células em comum
  rotation <- rotation[common_genes, , drop = FALSE]
  components <- components[common_cells, , drop = FALSE]
  original <- original[common_genes, common_cells, drop = FALSE]
  
  # Calcular a média das expressões gênicas (Tmed)
  Tmed <- rowMeans(original)
  Tmed[Tmed == 0] <- 1e-10  # Evitar divisão por zero
  
  # Número total de PCs
  num_pcs <- ncol(components)
  r2_values <- numeric(num_pcs)  # Vetor para armazenar os valores de R²
  
  # Inicializar a reconstrução
  reconstructed <- matrix(0, nrow(original), ncol(original))
  
  # Loop para calcular R² removendo PCs progressivamente
  for (i in 1:num_pcs) {
    cat("Calculando R² para", i, "PCs...\n")
    
    # Atualizar a reconstrução com a i-ésima PC
    reconstructed <- reconstructed + rotation[, i, drop = FALSE] %*% t(components[, i, drop = FALSE])
    
    # Calcular erro normalizado
    error_term <- (original - reconstructed)^2 / (Tmed^2)
    
    # Calcular R² para esse número de PCs
    r2_values[i] <- sum(error_term) / (ncol(original) * nrow(original))
  }
  
  cat("Cálculo de R² concluído!\n")
  
  # Plotar R², se solicitado
  if (plot) {
    plot_r2(r2_values)
  }
  
  # Retornar os valores de R²
  return(r2_values)
}

# Função auxiliar para plotar R²
plot_r2 <- function(r2_values) {
  plot(r2_values, type = "b", pch = 19, col = "steelblue",
       xlab = "Número de PCs", ylab = "R² Acumulado",
       main = "Qualidade da Reconstrução em Função do Número de PCs")
  grid()
}

# Calcular R² acumulativo para R0 e R30
r2_r0_cumulative <- calculate_r2_cumulative(t(df_t_matrix_batch_effect_correction_R0_for_PCA), pca_result_R0, plot = TRUE)
save(r2_r0_cumulative, file = "r2_r0_cumulative.RData")

r2_r30_cumulative <- calculate_r2_cumulative(t(df_t_matrix_batch_effect_correction_R30_for_PCA), pca_result_R30, plot = TRUE)
save(r2_r30_cumulative, file = "r2_r30_cumulative.RData")

# Criar gráfico de barras
library(ggplot2)
df_r2 <- data.frame(
  PCs = 1:length(r2_r0_cumulative),
  R2_R0 = r2_r0_cumulative,
  R2_R30 = r2_r30_cumulative
)

ggplot(df_r2, aes(x = PCs)) +
  geom_bar(aes(y = R2_R0, fill = "R0"), stat = "identity", alpha = 0.6) +
  geom_bar(aes(y = R2_R30, fill = "R30"), stat = "identity", alpha = 0.6) +
  labs(title = "Cumulative Variance by Number of PCs",
       x = "Number of PCs",
       y = "R²") +
  scale_fill_manual(values = c("R0" = "blue", "R30" = "red")) +
  theme_minimal()


##### PLOT DO R² #####

library(ggplot2)

# Carregar os resultados salvos
load("r2_r0_cumulative.RData")
load("r2_r30_cumulative.RData")

# Criar dataframes individuais para cada condição
df_r2_r0 <- data.frame(PC = 1:length(r2_r0_cumulative), R2 = r2_r0_cumulative)
df_r2_r30 <- data.frame(PC = 1:length(r2_r30_cumulative), R2 = r2_r30_cumulative)

# Gráfico para R0
plot_r0 <- ggplot(df_r2_r0, aes(x = PC, y = R2)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.6) +
  geom_point(size = 3, color = "blue") +
  geom_line(group = 1, color = "blue", linetype = "dashed") +
  labs(title = "Elbow Plot - R² Acumulado (R0)",
       x = "Componentes Principais", y = "R² Acumulado") +
  scale_x_continuous(breaks = seq(1, max(df_r2_r0$PC), by = round(max(df_r2_r0$PC) / 10))) +
  theme_minimal()

# Gráfico para R30
plot_r30 <- ggplot(df_r2_r30, aes(x = PC, y = R2)) +
  geom_bar(stat = "identity", fill = "firebrick", alpha = 0.6) +
  geom_point(size = 3, color = "red") +
  geom_line(group = 1, color = "red", linetype = "dashed") +
  labs(title = "Elbow Plot - R² Acumulado (R30)",
       x = "Componentes Principais", y = "R² Acumulado") +
  scale_x_continuous(breaks = seq(1, max(df_r2_r30$PC), by = round(max(df_r2_r30$PC) / 10))) +
  theme_minimal()

# Mostrar os gráficos
print(plot_r0)
print(plot_r30)
  



library(ggplot2)
library(patchwork)  # Para visualizar múltiplos gráficos juntos

# Carregar os resultados salvos
load("r2_r0_cumulative.RData")
load("r2_r30_cumulative.RData")

# Definir limites para Top 100 e Top 10 PCs
num_pcs_100 <- min(100, length(r2_r0_cumulative), length(r2_r30_cumulative))
num_pcs_10 <- min(10, length(r2_r0_cumulative), length(r2_r30_cumulative))

# Criar dataframes para cada condição e número de PCs
df_r2_r0_100 <- data.frame(PC = 1:num_pcs_100, R2 = r2_r0_cumulative[1:num_pcs_100])
df_r2_r30_100 <- data.frame(PC = 1:num_pcs_100, R2 = r2_r30_cumulative[1:num_pcs_100])

df_r2_r0_10 <- data.frame(PC = 1:num_pcs_10, R2 = r2_r0_cumulative[1:num_pcs_10])
df_r2_r30_10 <- data.frame(PC = 1:num_pcs_10, R2 = r2_r30_cumulative[1:num_pcs_10])

# Função para gerar gráficos
plot_r2 <- function(df, title, color) {
  ggplot(df, aes(x = PC, y = R2)) +
    geom_bar(stat = "identity", fill = color, alpha = 0.6) +
    geom_point(size = 3, color = color) +
    geom_line(group = 1, color = color, linetype = "dashed") +
    labs(title = title, x = "Componentes Principais", y = "R² Acumulado") +
    scale_x_continuous(breaks = seq(1, max(df$PC), by = max(1, round(max(df$PC) / 10)))) +
    theme_minimal()
}

# Gerar gráficos para cada condição e número de PCs
plot_r0_100 <- plot_r2(df_r2_r0_100, "Elbow Plot - R² Acumulado (R0, Top 100 PCs)", "steelblue")
plot_r30_100 <- plot_r2(df_r2_r30_100, "Elbow Plot - R² Acumulado (R30, Top 100 PCs)", "firebrick")

plot_r0_10 <- plot_r2(df_r2_r0_10, "Elbow Plot - R² Acumulado (R0, Top 10 PCs)", "steelblue")
plot_r30_10 <- plot_r2(df_r2_r30_10, "Elbow Plot - R² Acumulado (R30, Top 10 PCs)", "firebrick")

# Exibir os gráficos
print(plot_r0_100)
print(plot_r30_100)
print(plot_r0_10)
print(plot_r30_10)


#########################################################################################

library(ggplot2)

plot_transcriptogram <- function(df_original, df_reconstruido, max_pcs = 10) {
  # Converte para data frame
  df_original <- as.data.frame(df_original)
  df_reconstruido <- as.data.frame(df_reconstruido)
  
  # Garante que o número máximo de PCs não ultrapasse o total de colunas disponíveis
  max_pcs <- min(max_pcs, ncol(df_original))
  
  # Lista para armazenar diferenças individuais
  diff_data <- data.frame(PCs = integer(), Difference = numeric())
  
  for (num_pcs in 1:max_pcs) {
    # Seleciona apenas a PC correspondente (não acumulada)
    diff <- abs(df_original[, num_pcs] - df_reconstruido[, num_pcs])
    
    # Calcula a soma das diferenças dessa PC
    total_diff <- sum(diff, na.rm = TRUE)
    
    # Adiciona ao data frame
    diff_data <- rbind(diff_data, data.frame(PCs = num_pcs, Difference = total_diff))
  }
  
  # Plota
  ggplot(diff_data, aes(x = factor(PCs), y = Difference)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_minimal() +
    labs(title = "Diferença por PC entre Original e Reconstruído",
         x = "Número da PC",
         y = "Soma das Diferenças Absolutas")
}


plot_transcriptogram(df_t_matrix_batch_effect_correction_R0_for_PCA, reconstructed_r30, max_pcs = 10000)

################ GRÁFICO DE LINHAS DOS TRANSCRIPTOGRAMAS ######################

# Carregar pacotes necessários
library(ggplot2)

# Supondo que o dataframe tenha as colunas "PC1", "PC2", ... e uma variável de tempo "Day"
ggplot(df_t_matrix_batch_effect_correction_R0_for_PCA, aes(x = Day, y = PC1, group = 1)) +
  geom_line(color = "blue") +  # Adicionar a linha
  geom_point(color = "red") +  # Adicionar pontos para cada observação
  theme_minimal() +  # Tema minimalista
  labs(x = "Dia", y = "Componente Principal 1", title = "Gráfico de Linha de PCA")  # Rótulos dos eixos


plot(t(df_t_matrix_batch_effect_correction_R30_for_PCA)[, 2], type = "p", pch = 16, 
     col = "blue",  cex = 0.1, 
     xlab = "Índice", ylab = "Valor", main = "Gráfico de Pontos com Linhas")

# Adicionando linhas conectando os pontos
lines(t(reconstructed_r30)[, 2], type = "o",  cex = 0.1, col = "red")




