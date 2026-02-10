

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

# Caminhos para as amostras
paths <- c(
  "mtx_conversions/notreated_batch1/outs/filtered_feature_bc_matrix",
  "mtx_conversions/notreated_batch2/outs/filtered_feature_bc_matrix/",
  "mtx_conversions/TGFbeta1_1day_batch2/outs/filtered_feature_bc_matrix/",
  "mtx_conversions/TGFbeta1_2day_batch2/outs/filtered_feature_bc_matrix/",
  "mtx_conversions/TGFbeta1_3day_batch2/outs/filtered_feature_bc_matrix/",
  "mtx_conversions/TGFbeta1_4day_batch1/outs/filtered_feature_bc_matrix/",
  "mtx_conversions/TGFbeta1_8day_batch1/outs/filtered_feature_bc_matrix/"
)

# Nomes para os objetos
sample_names <- c(
  "notreated-batch1", "notreated-batch2", "TGFbeta1-1day-batch2",
  "TGFbeta1-2day-batch2", "TGFbeta1-3day-batch2", "TGFbeta1-4day-batch1", "TGFbeta1-8day-batch1"
)

# Criar objetos Seurat e adicionar anotações
seurat_objects <- lapply(1:length(paths), function(i) {
  data <- Read10X(paths[i])
  seurat <- CreateSeuratObject(counts = data, project = sample_names[i])
  seurat$batch <- ifelse(i %in% c(1, 6, 7), "Batch1", "Batch2")
  seurat$day <- c("0", "0", "1", "2", "3", "4", "8")[i]
  return(seurat)
})

# Identificar doublets e adicionar porcentagem de genes mitocondriais
seurat_objects <- lapply(seurat_objects, function(seurat) {
  # Converter para SingleCellExperiment
  sce <- as.SingleCellExperiment(seurat)
  sce <- scDblFinder(sce)
  
  # Adicionar a classificação de doublet ao metadado do objeto Seurat
  seurat$doublet <- sce$scDblFinder.class
  
  # Filtrar células não classificadas como "doublet"
  seurat <- subset(seurat, subset = doublet == "singlet")
  
  # Adicionar a porcentagem de genes mitocondriais
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
  
  return(seurat)
})

# Combinar todos os objetos Seurat
combined_seurat <- merge(
  x = seurat_objects[[1]],
  y = seurat_objects[-1],
  add.cell.ids = sample_names, # Adiciona os nomes para diferenciar os barcodes
  project = "Combined"
)

# Atualizar os colnames para o formato sample_name_barcode
colnames(combined_seurat) <- make.unique(
  paste0(
    sub("_.*", "", colnames(combined_seurat)), # Extrai o nome do sample
    "_",
    sub(".*_", "", colnames(combined_seurat)) # Extrai o barcode original
  )
)

# Conferir os colnames atualizados
head(colnames(combined_seurat))



# Aplicar controle de qualidade para cada objeto Seurat individual
filtered_seurat_objects <- lapply(seurat_objects, function(seurat) {
  subset(seurat, subset = nFeature_RNA > 500  & percent.mt < 20)
})

# Combinar os objetos filtrados
filtered_combined_seurat <- merge(
  filtered_seurat_objects[[1]],
  y = filtered_seurat_objects[-1],
  add.cell.ids = sample_names,
  project = "FilteredCombined"
)
# Conferir os colnames atualizados
head(colnames(filtered_combined_seurat))

########################################

# Visualizações antes do controle de qualidade
VlnPlot(
  combined_seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0,
  group.by = "orig.ident"
)


# Visualizações após o controle de qualidade
VlnPlot(
  filtered_combined_seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0,
  group.by = "orig.ident"
)

# Função para criar histogramas
plot_histogram <- function(seurat, title_prefix) {
  qc_data <- FetchData(seurat, vars = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
  p1 <- ggplot(qc_data, aes(x = nFeature_RNA)) +
    geom_histogram(bins = 100, fill = "blue", alpha = 0.7) +
    labs(title = paste(title_prefix, "- Distribuição de nFeature_RNA"), x = "nFeature_RNA", y = "Frequência") +
    theme_minimal()
  
  p2 <- ggplot(qc_data, aes(x = percent.mt)) +
    geom_histogram(bins = 100, fill = "orange", alpha = 0.7) +
    geom_vline(xintercept = 20, color = "red") +
    labs(title = paste(title_prefix, "- Percentagem de genes mitocondriais"), x = "percent.mt", y = "Frequência") +
    theme_minimal()
  
  list(nFeature_RNA = p1, percent_mt = p2)
}

# Histogramas antes do filtro
plots_before <- plot_histogram(combined_seurat, "Antes do Filtro")
plots_after <- plot_histogram(filtered_combined_seurat, "Depois do Filtro")

# Visualizar gráficos antes e depois
library(patchwork)
(plots_before$nFeature_RNA + plots_before$percent_mt) /
  (plots_after$nFeature_RNA + plots_after$percent_mt)

# Dados de QC antes do filtro
qc_before <- FetchData(combined_seurat, vars = c("nFeature_RNA", "nCount_RNA", "percent.mt", "orig.ident"))

# 1. Distribuição de nFeature_RNA
qc_before %>%
  ggplot(aes(x = nFeature_RNA)) +
  geom_histogram(bins = 200, fill = "blue", alpha = 0.7) +
  labs(title = "Antes do Filtro: Distribuição de nFeature_RNA", x = "nFeature_RNA", y = "Frequência") +
  theme_minimal()

# 2. Distribuição de percent.mt
qc_before %>%
  ggplot(aes(x = percent.mt)) +
  geom_histogram(bins = 100, fill = "orange", alpha = 0.7) +
  geom_vline(xintercept = 20, color = "red") +
  labs(title = "Antes do Filtro: Percentual de Genes Mitocondriais", x = "percent.mt", y = "Frequência") +
  theme_minimal()

# 3. Correlação nCount_RNA vs nFeature_RNA
qc_before %>%
  ggplot(aes(x = nCount_RNA, y = nFeature_RNA)) +
  geom_point(alpha = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "Antes do Filtro: nCount_RNA vs nFeature_RNA", x = "nCount_RNA (log10)", y = "nFeature_RNA (log10)") +
  theme_minimal()

# 4. Correlação nCount_RNA vs percent.mt
qc_before %>%
  ggplot(aes(x = nCount_RNA, y = percent.mt)) +
  geom_point(alpha = 0.5, color = "purple") +
  labs(title = "Antes do Filtro: nCount_RNA vs percent.mt", x = "nCount_RNA", y = "percent.mt") +
  theme_minimal()


# Dados de QC após o filtro
qc_after <- FetchData(filtered_combined_seurat, vars = c("nFeature_RNA", "nCount_RNA", "percent.mt", "orig.ident"))

# 1. Distribuição de nFeature_RNA
qc_after %>%
  ggplot(aes(x = nFeature_RNA)) +
  geom_histogram(bins = 200, fill = "blue", alpha = 0.7) +
  labs(title = "Depois do Filtro: Distribuição de nFeature_RNA", x = "nFeature_RNA", y = "Frequência") +
  theme_minimal()

# 2. Distribuição de percent.mt
qc_after %>%
  ggplot(aes(x = percent.mt)) +
  geom_histogram(bins = 100, fill = "orange", alpha = 0.7) +
  geom_vline(xintercept = 20, color = "red") +
  labs(title = "Depois do Filtro: Percentual de Genes Mitocondriais", x = "percent.mt", y = "Frequência") +
  theme_minimal()

# 3. Correlação nCount_RNA vs nFeature_RNA
qc_after %>%
  ggplot(aes(x = nCount_RNA, y = nFeature_RNA)) +
  geom_point(alpha = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "Depois do Filtro: nCount_RNA vs nFeature_RNA", x = "nCount_RNA (log10)", y = "nFeature_RNA (log10)") +
  theme_minimal()

# 4. Correlação nCount_RNA vs percent.mt
qc_after %>%
  ggplot(aes(x = nCount_RNA, y = percent.mt)) +
  geom_point(alpha = 0.5, color = "purple") +
  labs(title = "Depois do Filtro: nCount_RNA vs percent.mt", x = "nCount_RNA", y = "percent.mt") +
  theme_minimal()

# FeatureScatter: nCount_RNA vs nFeature_RNA
FeatureScatter(combined_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  labs(title = "Antes do Filtro: nCount_RNA vs nFeature_RNA")

FeatureScatter(filtered_combined_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  labs(title = "Depois do Filtro: nCount_RNA vs nFeature_RNA")

# FeatureScatter: nCount_RNA vs percent.mt
FeatureScatter(combined_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  labs(title = "Antes do Filtro: nCount_RNA vs percent.mt")

FeatureScatter(filtered_combined_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  labs(title = "Depois do Filtro: nCount_RNA vs percent.mt")

library(patchwork)

# Comparar distribuições de nFeature_RNA
before_nFeature <- qc_before %>%
  ggplot(aes(x = nFeature_RNA)) +
  geom_histogram(bins = 200, fill = "blue", alpha = 0.7) +
  labs(title = "Antes do Filtro", x = "nFeature_RNA", y = "Frequência") +
  theme_minimal()

after_nFeature <- qc_after %>%
  ggplot(aes(x = nFeature_RNA)) +
  geom_histogram(bins = 200, fill = "green", alpha = 0.7) +
  labs(title = "Depois do Filtro", x = "nFeature_RNA", y = "Frequência") +
  theme_minimal()

# Comparação lado a lado
before_nFeature + after_nFeature




# Extrair a lista de camadas
layers <- Layers(filtered_combined_seurat)

# Criar uma lista para armazenar as matrizes de expressão de cada camada
count_matrices <- lapply(layers, function(layer_name) {
  # Acessar os dados de contagem para cada camada
  GetAssayData(filtered_combined_seurat, assay = "RNA", slot = layer_name)
})

# Combinar as matrizes de contagem
combined_expression_matrix <- do.call(cbind, count_matrices)

# Verificar as dimensões da matriz combinada
dim(combined_expression_matrix)


combined_expression_matrix_full <- as.matrix(combined_expression_matrix)

###### Normalização Total Count 

# 1. Calcular a soma de cada coluna (Total Count por célula)
column_sums <- colSums(combined_expression_matrix_full)

# 2. Dividir cada elemento pela soma da respectiva coluna
normalized_matrix <- sweep(combined_expression_matrix_full, 2, column_sums, FUN = "/")

# 3. Verificar se a soma de cada coluna é 1
column_sums_normalized <- colSums(normalized_matrix)

# Exibir as primeiras somas para confirmação
head(column_sums_normalized)

# Retorna TRUE se todos os elementos do vetor lógico forem TRUE, ou FALSE caso contrário.
all(0.99999 < column_sums_normalized & column_sums_normalized < 1.000001)

# Salvando a matriz normalizada como arquivo
write.csv(normalized_matrix, "normalized_expression_matrix.csv")


## Means - Notreated Batch 1 e 2

# Filtrar colunas para "notreated_batch1" e "notreated_batch2"
cols_batch1 <- grep("^notreated-batch1", colnames(normalized_matrix))
cols_batch2 <- grep("^notreated-batch2", colnames(normalized_matrix))

# Submatrizes para cada batch
sub_matrix_batch1 <- normalized_matrix[, cols_batch1]
sub_matrix_batch2 <- normalized_matrix[, cols_batch2]

# Calcular as médias por linha
mean_batch1 <- rowMeans(sub_matrix_batch1)
mean_batch2 <- rowMeans(sub_matrix_batch2)

# Criar um data frame com os resultados
result <- data.frame(
  Gene = rownames(normalized_matrix),
  Mean_Batch1 = mean_batch1,
  Mean_Batch2 = mean_batch2,
  Ratio_Batch2_Batch1 = ifelse(mean_batch1 == 0, 0, mean_batch2 / mean_batch1) # Divisão das médias
)

head(result$Ratio_Batch2_Batch1)

## Multiplicando todas as colunas de normalized_matrix que possuem "batch1" no colname.

colnames_batch1 <- grep("batch1", colnames(normalized_matrix))

# criando cópia de normalized_matrix

copy_matrix <- normalized_matrix

# Multiplicando a razão pelas colunas selecionadas.

copy_matrix[, colnames_batch1] <- sweep(normalized_matrix[, colnames_batch1], 1, result$Ratio_Batch2_Batch1, "*") 

# Conectar ao banco de dados ENSEMBL
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

dictionary <- getBM(
  attributes = c("ensembl_peptide_id", "ensembl_gene_id"),
  mart = ensembl
)

# Limpar o dicionário, removendo valores vazios e NAs
dictionary <- dictionary %>%
  mutate(ensembl_peptide_id = ifelse(ensembl_peptide_id == "", NA, ensembl_peptide_id)) %>%
  na.omit()  # Remover linhas com NA


# Obter os mapeamentos para 'external_gene_name' e 'ensembl_gene_id'
gene_mapping <- getBM(
  attributes = c("external_gene_name", "ensembl_gene_id"),
  filters = "external_gene_name",
  values = rownames(copy_matrix),  # Usar os nomes dos genes
  mart = ensembl
)

# Limpar o mapeamento removendo valores vazios e NAs
gene_mapping <- gene_mapping %>%
  mutate(ensembl_gene_id = ifelse(ensembl_gene_id == "", NA, ensembl_gene_id)) %>%
  na.omit()  # Remover linhas com NA em 'ensembl_gene_id'

# Alinhar os nomes das linhas com o mapeamento
mapped_genes <- gene_mapping$ensembl_gene_id[match(
  rownames(copy_matrix),
  gene_mapping$external_gene_name
)]

# Filtrar para manter apenas as linhas com mapeamentos válidos
valid_indices <- !is.na(mapped_genes)
copy_matrix <- copy_matrix[valid_indices, , drop = FALSE]  # Filtrar linhas válidas
rownames(copy_matrix) <- mapped_genes[valid_indices]  # Substituir os nomes das linhas

# Verificar os primeiros nomes trocados
head(rownames(copy_matrix))  # Verifique os novos nomes de genes (ENSEMBL Gene IDs - ENSG)

# Ordering
ord <- vroom("ordering_HomoSapiensScore800-2024-C.txt")

# Association matrix
assoc <- read.opj("Associationmatrix.opj")
assoc <- assoc$associationMa

inner_join(assoc, ord, by = c("A" = "dim1")) %>% 
  dplyr::rename(protein1 = Protein) %>% 
  dplyr::inner_join(ord, by = c("B" = "dim1")) %>% 
  dplyr::rename(protein2 = Protein) %>% 
  dplyr::select(protein1, protein2) -> assoc


# Pré-processamento com raio = 0
t_matrix_batch_effect_correction_R0 <- transcriptogramPreprocess(association = assoc, 
                                                                 ordering = ord$Protein, 
                                                                 radius = 0)

# Rodar o Transcriptogramer
t_matrix_batch_effect_correction_R0 <- transcriptogramStep1(object = t_matrix_batch_effect_correction_R0, 
                                                            expression = copy_matrix, 
                                                            dictionary = dictionary)
t_matrix_batch_effect_correction_R0 <- transcriptogramStep2(object = t_matrix_batch_effect_correction_R0)

# Salvar os objetos processados
save(t_matrix_batch_effect_correction_R0, file = "t_matrix_batch_effect_correction_R0.RData")

# Pré-processamento com raio = 30
t_matrix_batch_effect_correction_R30 <- transcriptogramPreprocess(association = assoc, 
                                                                 ordering = ord$Protein, 
                                                                 radius = 30)

# Rodar o Transcriptogramer
t_matrix_batch_effect_correction_R30 <- transcriptogramStep1(object = t_matrix_batch_effect_correction_R30, 
                                                            expression = copy_matrix, 
                                                            dictionary = dictionary)
t_matrix_batch_effect_correction_R30 <- transcriptogramStep2(object = t_matrix_batch_effect_correction_R30)

# Salvar os objetos processados
save(t_matrix_batch_effect_correction_R30, file = "t_matrix_batch_effect_correction_R30.RData")

##################################################################
##### SUBSETs com 20 células de cada dia.

df_t_matrix_batch_effect_correction_R0 <- t_matrix_batch_effect_correction_R0@transcriptogramS2
df_t_matrix_batch_effect_correction_R30 <- t_matrix_batch_effect_correction_R30@transcriptogramS2

# Verificar se todos os elementos numéricos são não negativos no primeiro dataframe
all_non_negative_R0 <- all(df_t_matrix_batch_effect_correction_R0 >= 0)

# Verificar se todos os elementos numéricos são não negativos no segundo dataframe
all_non_negative_R30 <- all(df_t_matrix_batch_effect_correction_R30 >= 0)


# Exibir os resultados
cat("Todos os valores em df_t_matrix_batch_effect_correction_R0 são não negativos?", all_non_negative_R0, "\n")
cat("Todos os valores em df_t_matrix_batch_effect_correction_R30 são não negativos?", all_non_negative_R30, "\n")


# Lista de condições para filtrar
conditions <- c("notreated-batch1", "notreated-batch2", "TGFbeta1-1day-batch2", 
                "TGFbeta1-2day-batch2", "TGFbeta1-3day-batch2", 
                "TGFbeta1-4day-batch1", "TGFbeta1-8day-batch1")

filter_barcodes <- function(data, conditions) {
  # Separar as duas primeiras colunas
  fixed_columns <- data[, 1:2, drop = FALSE]
  
  # Filtrar barcodes por condição, excluindo as duas primeiras colunas
  filtered_data <- lapply(conditions, function(cond) {
    subset <- data[, grepl(cond, colnames(data)), drop = FALSE]
    if (ncol(subset) > 20) {
      subset <- subset[, sample(ncol(subset), 20), drop = FALSE]
    }
    return(subset)
  })
  
  # Combine os subconjuntos filtrados
  combined_data <- do.call(cbind, filtered_data)
  
  # Juntar as duas primeiras colunas fixas com os barcodes filtrados
  result <- cbind(fixed_columns, combined_data)
  return(result)
}



# Aplicar a função para os dois datasets
df_t_matrix_batch_effect_correction_R0_filtered <- filter_barcodes(df_t_matrix_batch_effect_correction_R0, conditions)
df_t_matrix_batch_effect_correction_R30_filtered <- filter_barcodes(df_t_matrix_batch_effect_correction_R30, conditions)

# Exportar os resultados para TSV

write.table(df_t_matrix_batch_effect_correction_R0_filtered, "df_t_matrix_batch_effect_correction_R0_filtered.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(df_t_matrix_batch_effect_correction_R30_filtered, "df_t_matrix_batch_effect_correction_R30_filtered.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)




##############
#### PCA
##############


# Supondo que os dataframes estejam limpos e prontos para a análise (remover o ENSP ou outras colunas irrelevantes)
df_t_matrix_batch_effect_correction_R0_for_PCA <- t_matrix_batch_effect_correction_R0@transcriptogramS2[,-2]  # Remover coluna de identificação
df_t_matrix_batch_effect_correction_R30_for_PCA <- t_matrix_batch_effect_correction_R30@transcriptogramS2[,-2]

# Para o data frame R0
rownames(df_t_matrix_batch_effect_correction_R0_for_PCA) <- df_t_matrix_batch_effect_correction_R0_for_PCA$Protein
df_t_matrix_batch_effect_correction_R0_for_PCA <- df_t_matrix_batch_effect_correction_R0_for_PCA[, -1]

# Para o data frame R30
rownames(df_t_matrix_batch_effect_correction_R30_for_PCA) <- df_t_matrix_batch_effect_correction_R30_for_PCA$Protein
df_t_matrix_batch_effect_correction_R30_for_PCA <- df_t_matrix_batch_effect_correction_R30_for_PCA[, -1]


########## PCA e plots


# Remover colunas com rótulos NA (se existirem)
df_t_matrix_batch_effect_correction_R0_for_PCA <- df_t_matrix_batch_effect_correction_R0_for_PCA[, !is.na(colnames(df_t_matrix_batch_effect_correction_R0_for_PCA))]

# Transpor a matriz para realizar o PCA
pca_input_R0 <- t(df_t_matrix_batch_effect_correction_R0_for_PCA)



# Verificar a variância de cada coluna 
variances_R0 <- apply(pca_input_R0, 2, var)

# Remover colunas com variância zero
pca_input_R0_cleaned <- pca_input_R0[, variances_R0 > 0]

# Rodar PCA
pca_result_R0 <- prcomp(pca_input_R0_cleaned, center = TRUE, scale. = FALSE)

# Variância acumulada
variance_R0 <- summary(pca_result_R0)$importance[2, ]
cumulative_variance_R0 <- cumsum(variance_R0)

# Criar data.frame com variância
variance_df_R0 <- data.frame(
  PC = seq_along(cumulative_variance_R0),
  Variance = variance_R0,
  CumulativeVariance = cumulative_variance_R0
)

# Salvar a matriz PCA como .tsv
pca_matrix_R0 <- as.data.frame(pca_result_R0$x)
write.table(pca_matrix_R0, file = paste0("R0_pca_matrix.tsv"), 
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

# Gráfico de variância acumulada
variance_plot_R0 <- ggplot(variance_df_R0, aes(x = PC)) +
  geom_line(aes(y = CumulativeVariance, color = "Cumulative Variance")) +
  geom_point(aes(y = CumulativeVariance, color = "Cumulative Variance")) +
  scale_color_manual(values = c("Cumulative Variance" = "blue")) +
  labs(title = "Cumulative Variance by Principal Components (R0)", 
       x = "Principal Component", y = "Variance Explained (%)") +
  theme_minimal()

# Salvar o gráfico
ggsave(filename = paste0("images/R0_variance_plot.png"), plot = variance_plot_R0, dpi = 200)

# Gráfico para as 10 primeiras PCs
top10_plot_R0 <- ggplot(variance_df_R0[1:10, ], aes(x = PC)) +
  geom_line(aes(y = CumulativeVariance, color = "Cumulative Variance")) +
  geom_point(aes(y = CumulativeVariance, color = "Cumulative Variance")) +
  scale_color_manual(values = c("Cumulative Variance" = "green")) +
  labs(title = "Cumulative Variance - Top 10 PCs (R0)", 
       x = "Principal Component", y = "Variance Explained (%)") +
  theme_minimal()

# Salvar o gráfico
ggsave(filename = paste0("images/R0_top10_variance_plot.png"), plot = top10_plot_R0, dpi = 300)


# Criar a função para extrair a parte anterior ao primeiro '_'
condition_labels_R0 <- sapply(colnames(df_t_matrix_batch_effect_correction_R0_for_PCA), function(col) {
  match <- grep(paste(conditions, collapse = "|"), col, value = TRUE)
  if (length(match) > 0) {
    # Encontrar a parte anterior ao primeiro '_'
    prefix_term <- sub("_.*", "", match[1])
    return(prefix_term)
  } else {
    return(NA)
  }
})

# Verificar as condições extraídas
unique(condition_labels_R0)


# Criar gráficos pairwise para PCs 1 a 4
pairwise_plots_R0 <- list()
pcs_to_plot_R0 <- 1:4
color_palette_R0 <- RColorBrewer::brewer.pal(n = length(unique(condition_labels_R0)), "Set1")

# Verificar alinhamento
pca_matrix_R0$Condition <- factor(condition_labels_R0[rownames(pca_matrix_R0)])

counter_R0 <- 1
for (i in pcs_to_plot_R0) {
  for (j in pcs_to_plot_R0) {
    if (i < j) {
      pairwise_plots_R0[[counter_R0]] <- ggplot(pca_matrix_R0, aes(x = !!sym(paste0("PC", i)), 
                                                                   y = !!sym(paste0("PC", j)), 
                                                                   color = Condition)) +
        geom_point(size = 1.5) +
        labs(title = paste0("PC", i, " vs PC", j), x = paste0("PC", i), y = paste0("PC", j)) +
        theme_minimal() +
        theme(legend.position = "right") +
        scale_color_manual(values = color_palette_R0)
      counter_R0 <- counter_R0 + 1
    }
  }
}

# Salvar o grid com pairwise plots
pairwise_plots_R0 <- pairwise_plots_R0[!sapply(pairwise_plots_R0, is.null)]

# Combinar os gráficos
combined_plot_R0 <- wrap_plots(pairwise_plots_R0, ncol = 2)

# Salvar o gráfico combinado
ggsave(
  filename = "images/R0_pairwise_plots_patchwork.png", 
  plot = combined_plot_R0, 
  dpi = 300,
  width = 12, 
  height = 8
)

# Gráfico de PCA 2D (PC1 vs PC2)
pca_2d_plot_R0 <- ggplot(pca_matrix_R0, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 2, alpha = 0.7) +
  stat_ellipse(aes(group = Condition), type = "norm", level = 0.95, linetype = "dashed", linewidth = 0.5) +
  labs(title = "PCA (PC1 vs PC2) - R0", x = "PC1", y = "PC2") +
  theme_minimal() +
  scale_color_manual(values = color_palette_R0) +
  theme(legend.position = "right")

ggsave(
  filename = "images/R0_pca_2d_plot_improved.png", 
  plot = pca_2d_plot_R0, 
  dpi = 300, 
  width = 10, 
  height = 8
)

# Selecionar as 10 primeiras PCs
top10_variance_df_R0 <- variance_df_R0[1:10, ]

# Gráfico para as 10 primeiras PCs com variância acumulada
explained_variance_top10_plot_R0 <- ggplot(top10_variance_df_R0, aes(x = PC)) +
  geom_bar(aes(y = Variance * 100), stat = "identity", fill = "steelblue", alpha = 0.8) +  # Variância explicada
  geom_line(aes(y = CumulativeVariance * 100, group = 1, color = "Cumulative Variance"), size = 1) +  # Variância acumulada
  geom_point(aes(y = CumulativeVariance * 100, color = "Cumulative Variance"), size = 2) +  # Pontos da variância acumulada
  geom_text(aes(y = Variance * 100, label = round(Variance * 100, 2)), vjust = -0.5, size = 3) +  # Rótulos das barras
  labs(title = "Explained and Cumulative Variance by Top 10 PCs (R0)", 
       x = "Principal Component", 
       y = "Variance Explained (%)") +
  scale_x_continuous(breaks = 1:10) +  # Apenas as 10 primeiras PCs
  scale_y_continuous(limits = c(0, 100)) +  # Ajustar o eixo Y para ir de 0 a 100%
  scale_color_manual(values = c("Cumulative Variance" = "darkorange")) +  # Cor da linha acumulada
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotação dos rótulos
        legend.title = element_blank(), 
        legend.position = "right")

# Salvar o gráfico
ggsave(filename = "images/R0_explained_cumulative_variance_top10_plot_fixed.png", 
       plot = explained_variance_top10_plot_R0, 
       dpi = 300, 
       width = 10, 
       height = 6)

# Selecionar as 100 primeiras PCs
top100_variance_df_R0 <- variance_df_R0[1:100, ]

# Gráfico para as 100 primeiras PCs com variância acumulada
explained_variance_top100_plot_R0 <- ggplot(top100_variance_df_R0, aes(x = PC)) +
  geom_bar(aes(y = Variance * 100), stat = "identity", fill = "steelblue", alpha = 0.8) +  # Variância explicada
  geom_line(aes(y = CumulativeVariance * 100, group = 1, color = "Cumulative Variance"), size = 1) +  # Variância acumulada
  geom_point(aes(y = CumulativeVariance * 100, color = "Cumulative Variance"), size = 2) +  # Pontos da variância acumulada
  geom_text(aes(y = Variance * 100, label = ifelse(Variance * 100 > 1, round(Variance * 100, 2), "")), 
            vjust = -0.5, size = 2) +  # Rótulos das barras apenas para variâncias >1%
  labs(title = "Explained and Cumulative Variance by Top 100 PCs (R0)", 
       x = "Principal Component", 
       y = "Variance Explained (%)") +
  scale_x_continuous(breaks = seq(1, 100, by = 10)) +  # Marcar PCs a cada 10
  scale_y_continuous(limits = c(0, 100)) +  # Ajustar o eixo Y para 0-100%
  scale_color_manual(values = c("Cumulative Variance" = "darkorange")) +  # Cor da linha acumulada
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotação dos rótulos
        legend.title = element_blank(), 
        legend.position = "right")

# Salvar o gráfico
ggsave(filename = "images/R0_explained_cumulative_variance_top100_plot_fixed.png", 
       plot = explained_variance_top100_plot_R0, 
       dpi = 300, 
       width = 12, 
       height = 6)

# Gráfico para todas as PCs com variância acumulada
explained_variance_all_plot_R0 <- ggplot(variance_df_R0, aes(x = PC)) +
  geom_bar(aes(y = Variance * 100), stat = "identity", fill = "steelblue", alpha = 0.8) +  # Variância explicada
  geom_line(aes(y = CumulativeVariance * 100, group = 1, color = "Cumulative Variance"), size = 1) +  # Variância acumulada
  geom_point(aes(y = CumulativeVariance * 100, color = "Cumulative Variance"), size = 1) +  # Pontos da variância acumulada
  geom_text(aes(y = Variance * 100, label = ifelse(Variance * 100 > 2, round(Variance * 100, 2), "")), 
            vjust = -0.5, size = 2) +  # Rótulos das barras apenas para variâncias >2%
  labs(title = "Explained and Cumulative Variance by All PCs (R0)", 
       x = "Principal Component", 
       y = "Variance Explained (%)") +
  scale_x_continuous(breaks = seq(1, nrow(variance_df_R0), by = 10)) +  # Marcar PCs a cada 10
  scale_y_continuous(limits = c(0, 100)) +  # Ajustar o eixo Y para 0-100%
  scale_color_manual(values = c("Cumulative Variance" = "darkorange")) +  # Cor da linha acumulada
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotação dos rótulos
        legend.title = element_blank(), 
        legend.position = "right")

# Salvar o gráfico
ggsave(filename = "images/R0_explained_cumulative_variance_all_plot_fixed.png", 
       plot = explained_variance_all_plot_R0, 
       dpi = 300, 
       width = 15, 
       height = 6)



# Continuando o processo com os mesmos passos, mas para o dataset R30.

# Criar a função para extrair a parte anterior ao primeiro '_'
condition_labels_R30 <- sapply(colnames(df_t_matrix_batch_effect_correction_R30_for_PCA), function(col) {
  match <- grep(paste(conditions, collapse = "|"), col, value = TRUE)
  if (length(match) > 0) {
    # Encontrar a parte anterior ao primeiro '_'
    prefix_term <- sub("_.*", "", match[1])
    return(prefix_term)
  } else {
    return(NA)
  }
})

# Verificar as condições extraídas
unique(condition_labels_R30)

# # Substituir colnames por condições
# colnames(df_t_matrix_batch_effect_correction_R30_for_PCA) <- condition_labels_R30

# Remover colunas com rótulos NA (se existirem)
df_t_matrix_batch_effect_correction_R30_for_PCA <- df_t_matrix_batch_effect_correction_R30_for_PCA[, !is.na(colnames(df_t_matrix_batch_effect_correction_R30_for_PCA))]

# Transpor a matriz para realizar o PCA
pca_input_R30 <- t(df_t_matrix_batch_effect_correction_R30_for_PCA)

# Verificar a variância de cada coluna
variances_R30 <- apply(pca_input_R30, 2, var)

# Remover colunas com variância zero
pca_input_R30_cleaned <- pca_input_R30[, variances_R30 > 0]

# Rodar PCA
pca_result_R30 <- prcomp(pca_input_R30_cleaned, center = TRUE, scale. = FALSE)

# Variância acumulada
variance_R30 <- summary(pca_result_R30)$importance[2, ]
cumulative_variance_R30 <- cumsum(variance_R30)

# Criar data.frame com variância
variance_df_R30 <- data.frame(
  PC = seq_along(cumulative_variance_R30),
  Variance = variance_R30,
  CumulativeVariance = cumulative_variance_R30
)

# Salvar a matriz PCA como .tsv
pca_matrix_R30 <- as.data.frame(pca_result_R30$x)
write.table(pca_matrix_R30, file = paste0("R30_pca_matrix.tsv"), 
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

# Gráfico de variância acumulada
variance_plot_R30 <- ggplot(variance_df_R30, aes(x = PC)) +
  geom_line(aes(y = CumulativeVariance, color = "Cumulative Variance")) +
  geom_point(aes(y = CumulativeVariance, color = "Cumulative Variance")) +
  scale_color_manual(values = c("Cumulative Variance" = "blue")) +
  labs(title = "Cumulative Variance by Principal Components (R30)", 
       x = "Principal Component", y = "Variance Explained (%)") +
  theme_minimal()

# Salvar o gráfico
ggsave(filename = paste0("images/R30_variance_plot.png"), plot = variance_plot_R30, dpi = 200)

# Gráfico para as 10 primeiras PCs
top10_plot_R30 <- ggplot(variance_df_R30[1:10, ], aes(x = PC)) +
  geom_line(aes(y = CumulativeVariance, color = "Cumulative Variance")) +
  geom_point(aes(y = CumulativeVariance, color = "Cumulative Variance")) +
  scale_color_manual(values = c("Cumulative Variance" = "green")) +
  labs(title = "Cumulative Variance - Top 10 PCs (R30)", 
       x = "Principal Component", y = "Variance Explained (%)") +
  theme_minimal()

# Salvar o gráfico
ggsave(filename = paste0("images/R30_top10_variance_plot.png"), plot = top10_plot_R30, dpi = 300)

# Criar gráficos pairwise para PCs 1 a 4
pairwise_plots_R30 <- list()
pcs_to_plot_R30 <- 1:4
color_palette_R30 <- RColorBrewer::brewer.pal(n = length(unique(condition_labels_R30)), "Set1")

# Verificar alinhamento
pca_matrix_R30$Condition <- factor(condition_labels_R30[rownames(pca_matrix_R30)])

counter_R30 <- 1
for (i in pcs_to_plot_R30) {
  for (j in pcs_to_plot_R30) {
    if (i < j) {
      pairwise_plots_R30[[counter_R30]] <- ggplot(pca_matrix_R30, aes(x = !!sym(paste0("PC", i)), 
                                                                      y = !!sym(paste0("PC", j)), 
                                                                      color = Condition)) +
        geom_point(size = 1.5) +
        labs(title = paste0("PC", i, " vs PC", j), x = paste0("PC", i), y = paste0("PC", j)) +
        theme_minimal() +
        theme(legend.position = "right") +
        scale_color_manual(values = color_palette_R30)
      counter_R30 <- counter_R30 + 1
    }
  }
}

# Salvar o grid com pairwise plots
pairwise_plots_R30 <- pairwise_plots_R30[!sapply(pairwise_plots_R30, is.null)]

# Combinar os gráficos
combined_plot_R30 <- wrap_plots(pairwise_plots_R30, ncol = 2)

# Salvar o gráfico combinado
ggsave(
  filename = "images/R30_pairwise_plots_patchwork.png", 
  plot = combined_plot_R30, 
  dpi = 300,
  width = 12, 
  height = 8
)



pca_2d_plot_R30 <- ggplot(pca_matrix_R30, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 2, alpha = 0.7) +
  stat_ellipse(aes(group = Condition), type = "norm", level = 0.95, linetype = "dashed", linewidth = 0.5) +
  labs(title = "PCA (PC1 vs PC2) - R30", x = "PC1", y = "PC2") +
  theme_minimal() +
  scale_color_manual(values = color_palette_R30) +
  theme(legend.position = "right")

ggsave(
  filename = "images/R30_pca_2d_plot_improved.png", 
  plot = pca_2d_plot_R30, 
  dpi = 300, 
  width = 10, 
  height = 8
)
# Filtrar para as 10 primeiras PCs
top10_variance_df_R30 <- variance_df_R30[1:10, ]

# Gráfico para as 10 primeiras PCs com variância acumulada
explained_variance_top10_plot_R30 <- ggplot(top10_variance_df_R30, aes(x = PC)) +
  geom_bar(aes(y = Variance * 100), stat = "identity", fill = "steelblue", alpha = 0.8) +  # Variância explicada
  geom_line(aes(y = CumulativeVariance * 100, group = 1, color = "Cumulative Variance"), size = 1) +  # Variância acumulada
  geom_point(aes(y = CumulativeVariance * 100, color = "Cumulative Variance"), size = 2) +  # Pontos da variância acumulada
  geom_text(aes(y = Variance * 100, label = round(Variance * 100, 2)), vjust = -0.5, size = 3) +  # Rótulos das barras
  labs(title = "Explained and Cumulative Variance by Top 10 PCs (R30)", 
       x = "Principal Component", 
       y = "Variance Explained (%)") +
  scale_x_continuous(breaks = 1:10) +  # Apenas as 10 primeiras PCs
  scale_y_continuous(limits = c(0, 100)) +  # Ajustar o eixo Y para ir de 0 a 100%
  scale_color_manual(values = c("Cumulative Variance" = "darkorange")) +  # Cor da linha acumulada
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotação dos rótulos
        legend.title = element_blank(), 
        legend.position = "right")

# Salvar o gráfico
ggsave(filename = "images/R30_explained_cumulative_variance_top10_plot_fixed.png", 
       plot = explained_variance_top10_plot_R30, 
       dpi = 300, 
       width = 10, 
       height = 6)

# Selecionar as 100 primeiras PCs
top100_variance_df_R30 <- variance_df_R30[1:100, ]

# Gráfico para as 100 primeiras PCs com variância acumulada
explained_variance_top100_plot_R30 <- ggplot(top100_variance_df_R30, aes(x = PC)) +
  geom_bar(aes(y = Variance * 100), stat = "identity", fill = "steelblue", alpha = 0.8) +  # Variância explicada
  geom_line(aes(y = CumulativeVariance * 100, group = 1, color = "Cumulative Variance"), size = 1) +  # Variância acumulada
  geom_point(aes(y = CumulativeVariance * 100, color = "Cumulative Variance"), size = 2) +  # Pontos da variância acumulada
  geom_text(aes(y = Variance * 100, label = ifelse(Variance * 100 > 1, round(Variance * 100, 2), "")), 
            vjust = -0.5, size = 2) +  # Rótulos das barras apenas para variâncias >1%
  labs(title = "Explained and Cumulative Variance by Top 100 PCs (R30)", 
       x = "Principal Component", 
       y = "Variance Explained (%)") +
  scale_x_continuous(breaks = seq(1, 100, by = 10)) +  # Marcar PCs a cada 10
  scale_y_continuous(limits = c(0, 100)) +  # Ajustar o eixo Y para 0-100%
  scale_color_manual(values = c("Cumulative Variance" = "darkorange")) +  # Cor da linha acumulada
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotação dos rótulos
        legend.title = element_blank(), 
        legend.position = "right")

# Salvar o gráfico
ggsave(filename = "images/R30_explained_cumulative_variance_top100_plot_fixed.png", 
       plot = explained_variance_top100_plot_R30, 
       dpi = 300, 
       width = 12, 
       height = 6)



# Gráfico para todas as PCs com variância acumulada
explained_variance_all_plot_R30 <- ggplot(variance_df_R30, aes(x = PC)) +
  geom_bar(aes(y = Variance * 100), stat = "identity", fill = "steelblue", alpha = 0.8) +  # Variância explicada
  geom_line(aes(y = CumulativeVariance * 100, group = 1, color = "Cumulative Variance"), size = 1) +  # Variância acumulada
  geom_point(aes(y = CumulativeVariance * 100, color = "Cumulative Variance"), size = 2) +  # Pontos da variância acumulada
  geom_text(aes(y = Variance * 100, label = ifelse(Variance * 100 > 2, round(Variance * 100, 2), "")), 
            vjust = -0.5, size = 2) +  # Rótulos das barras apenas para variâncias >2%
  labs(title = "Explained and Cumulative Variance by All PCs (R30)", 
       x = "Principal Component", 
       y = "Variance Explained (%)") +
  scale_x_continuous(breaks = seq(1, nrow(variance_df_R30), by = 5)) +  # Marcação em intervalos de 5
  scale_y_continuous(limits = c(0, 100)) +  # Ajustar o eixo Y para 0-100%
  scale_color_manual(values = c("Cumulative Variance" = "darkorange")) +  # Cor da linha acumulada
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotação dos rótulos
        legend.title = element_blank(), 
        legend.position = "right")

# Salvar o gráfico
ggsave(filename = "images/R30_explained_cumulative_variance_all_plot_fixed.png", 
       plot = explained_variance_all_plot_R30, 
       dpi = 300, 
       width = 12, 
       height = 6)

#################
# PLOTS AGRUPADOS
#################

# Combinar os dois gráficos
combined_plot_top10 <- explained_variance_top10_plot_R30 / explained_variance_top10_plot_R0

# Salvar o gráfico combinado
ggsave(filename = "images/combined_explained_cumulative_variance_top10.png", 
       plot = combined_plot_top10, 
       dpi = 300, 
       width = 12, 
       height = 6)


#########################################################
########### ROTATION & COMPONENT
#########################################################

# Extrair as 10 primeiras componentes de rotação e componente principal de R0
rotation_R0_10PCs <- as.data.frame(pca_result_R0$rotation[, 1:10])
component_R0_10PCs <- as.data.frame(pca_result_R0$x[, 1:10])

# Salvar os resultados de R0
write.table(rotation_R0_10PCs, "rotation_R0_10PCs.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(component_R0_10PCs, "component_R0_10PCs.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

# Extrair as 10 primeiras componentes de rotação e componente principal de R30
rotation_R30_10PCs <- as.data.frame(pca_result_R30$rotation[, 1:10])
component_R30_10PCs <- as.data.frame(pca_result_R30$x[, 1:10])

# Salvar os resultados de R30
write.table(rotation_R30_10PCs, "rotation_R30_10PCs.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(component_R30_10PCs, "component_R30_10PCs.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

# Exibir uma amostra para verificação
head(rotation_R0_10PCs)
head(component_R0_10PCs)
head(rotation_R30_10PCs)
head(component_R30_10PCs)

## Colocando a coluna Position no dataframe rotation_R30_10PCs
# Verificando se as dimensões são compatíveis
if (nrow(rotation_R30_10PCs) == nrow(df_t_matrix_batch_effect_correction_R30_filtered)) {
  # Adicionando a coluna Position na primeira posição
  rotation_R30_10PCs <- cbind(Position = df_t_matrix_batch_effect_correction_R30_filtered$Position, 
                              rotation_R30_10PCs)
} else {
  stop("O número de linhas nos dataframes não corresponde!")
}

# Visualizando as primeiras linhas do novo dataframe
head(rotation_R30_10PCs)


library(ggplot2)
library(patchwork)

# Definir limites do eixo Y com base no intervalo geral
y_limits <- range(as.matrix(rotation_R30_10PCs[,-1]), na.rm = TRUE)

# Criar lista de gráficos
plots <- lapply(2:ncol(rotation_R30_10PCs), function(i) {
  pc_name <- colnames(rotation_R30_10PCs)[i]
  
  ggplot(rotation_R30_10PCs, aes(x = Position, y = .data[[pc_name]])) +
    geom_line(color = "blue", size = 0.3) + # Linha conectando os pontos
    scale_x_continuous(limits = c(0, 14480), breaks = seq(0, 14480, by = 2000)) + # Limite eixo X até 14480
    scale_y_continuous(limits = y_limits) + # Eixo Y consistente
    labs(title = paste("PCA:", pc_name), x = "Position", y = pc_name) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 10)
    )
})

# Combinar gráficos em layout de grade
output_plot <- wrap_plots(plots, ncol = 2) # Dois gráficos por linha

# Salvar em alta resolução
ggsave("images/PCA_Position_Line_HighRes.png", plot = output_plot, width = 20, height = 24, dpi = 600)


# Ordenando as proteínas pela contribuição na PC1
important_proteins_PC1 <- rotation_R30_10PCs[order(abs(rotation_R30_10PCs$PC1), decreasing = TRUE), ]
top_important_proteins <- head(important_proteins_PC1, 10)  # Top 10 proteínas
protein_list <- rownames(top_important_proteins)

library(biomaRt)

# Conectar ao banco de dados Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Obter Gene Symbols a partir dos ENSP IDs
converted_genes <- getBM(
  attributes = c("ensembl_peptide_id", "hgnc_symbol", "entrezgene_id"),
  filters = "ensembl_peptide_id",
  values = protein_list,  # Usar a lista de ENSP
  mart = ensembl
)

# Verifique os dados convertidos
head(converted_genes)



library(org.Hs.eg.db)
library(clusterProfiler)
# Para enriquecimento de GO
go_enrich <- enrichGO(
  gene = converted_genes$hgnc_symbol, 
  OrgDb = org.Hs.eg.db, 
  ont = "BP", 
  keyType = "SYMBOL"  # Especifica que os genes são símbolos HGNC
)


# Para enriquecimento de KEGG
# Remover valores NA em Entrez IDs
entrez_ids <- na.omit(converted_genes$entrezgene_id)

kegg_enrich <- enrichKEGG(
  gene = entrez_ids, 
  organism = "hsa"  # Humano
)

# Resultados de GO
if (!is.null(go_enrich)) {
  print(head(go_enrich))
} else {
  cat("Nenhum termo de GO enriquecido encontrado.\n")
}

# Resultados de KEGG
if (!is.null(kegg_enrich)) {
  print(head(kegg_enrich))
} else {
  cat("Nenhum termo KEGG enriquecido encontrado.\n")
}

library(enrichplot)

# Dot plot para GO
if (!is.null(go_enrich)) {
  dotplot(go_enrich) + ggtitle("Enriquecimento GO")
}

# Dot plot para KEGG
if (!is.null(kegg_enrich)) {
  dotplot(kegg_enrich) + ggtitle("Enriquecimento KEGG")
}


library(ReactomePA)
# Remover valores NA
entrez_ids <- na.omit(converted_genes$entrezgene_id)

# Verificar os primeiros valores
head(entrez_ids)

reactome_enrich <- enrichPathway(
  gene = entrez_ids, 
  organism = "human"
)

if (!is.null(reactome_enrich)) {
  print(head(reactome_enrich))
} else {
  cat("Nenhum termo Reactome enriquecido encontrado.\n")
}

# Dot plot para Reactome
if (!is.null(reactome_enrich)) {
  dotplot(reactome_enrich) + ggtitle("Enriquecimento Reactome")
}

# Bar plot para Reactome
if (!is.null(reactome_enrich)) {
  barplot(reactome_enrich, showCategory = 10, title = "Enriquecimento Reactome")
}

library(STRINGdb)
string_db <- STRINGdb$new(version = "11", species = 9606, score_threshold = 400)
string_network <- string_db$get_interactions(entrez_ids)

library(igraph)
gene_network <- graph_from_data_frame(string_network)
plot(gene_network, vertex.size = 5, vertex.label.cex = 0.7)


library(VennDiagram)
venn.plot <- venn.diagram(
  x = list(Reactome = reactome_enrich$ID, GO = go_enrich$ID), 
  filename = NULL, 
  fill = c("red", "blue")
)
grid.draw(venn.plot)

library(clusterProfiler)
install.packages("GOSemSim")

# Criar a matriz de similaridade semântica
sim_matrix <- pairwise_termsim(go_enrich, method = "Wang")

# Clusterizar termos redundantes
clusters <- hclust(as.dist(1 - sim_matrix))  # Clusterização hierárquica
cutoff <- 0.7  # Limite para corte
term_clusters <- cutree(clusters, h = cutoff)

# Selecionar termos representativos para cada cluster
selected_terms <- sapply(unique(term_clusters), function(cluster) {
  terms_in_cluster <- names(term_clusters[term_clusters == cluster])
  # Seleciona o termo com menor p.adjust
  terms_in_cluster[which.min(go_enrich@result$p.adjust[terms_in_cluster])]
})

# Filtrar resultados simplificados
simplified_results <- go_enrich@result[rownames(go_enrich@result) %in% selected_terms, ]

# Criar novo objeto de enriquecimento GO simplificado
go_enrich_simplified <- go_enrich
go_enrich_simplified@result <- simplified_results

# Visualizar os resultados simplificados
dotplot(go_enrich_simplified) + ggtitle("Enriquecimento GO Simplificado")




