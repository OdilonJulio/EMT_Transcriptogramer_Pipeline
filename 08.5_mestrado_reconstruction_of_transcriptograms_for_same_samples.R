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

##### PROCESS BATCH EFFECT CORRECT

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
# Ajuste na matriz copy_matrix
copy_matrix[, colnames_batch1] <- sweep(normalized_matrix[, colnames_batch1], 1, result$Ratio_Batch2_Batch1, "*") 

# Conectar ao banco de dados ENSEMBL
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Obter os mapeamentos de ENSP para ENSG
dictionary <- getBM(
  attributes = c("ensembl_peptide_id", "ensembl_gene_id"),
  mart = ensembl
)

# Limpar o dicionário, removendo valores vazios e NAs
dictionary <- dictionary %>%
  mutate(ensembl_peptide_id = ifelse(ensembl_peptide_id == "", NA, ensembl_peptide_id)) %>%
  na.omit()  # Remover linhas com NA

# Obter os mapeamentos para 'external_gene_name' e 'ensembl_peptide_id'
gene_mapping <- getBM(
  attributes = c("external_gene_name", "ensembl_peptide_id"),
  filters = "external_gene_name",
  values = rownames(copy_matrix),  # Usar os nomes dos genes
  mart = ensembl
)

# Limpar o mapeamento removendo valores vazios e NAs
gene_mapping <- gene_mapping %>%
  mutate(ensembl_peptide_id = ifelse(ensembl_peptide_id == "", NA, ensembl_peptide_id)) %>%
  na.omit()  # Remover linhas com NA em 'ensembl_peptide_id'

# Alinhar os nomes das linhas com o mapeamento ENSP
mapped_proteins <- gene_mapping$ensembl_peptide_id[match(
  rownames(copy_matrix),
  gene_mapping$external_gene_name
)]

# Filtrar para manter apenas as linhas com mapeamentos válidos
valid_indices <- !is.na(mapped_proteins)
copy_matrix <- copy_matrix[valid_indices, , drop = FALSE]  # Filtrar linhas válidas
rownames(copy_matrix) <- mapped_proteins[valid_indices]  # Substituir os nomes das linhas

# Verificar os primeiros nomes trocados
head(rownames(copy_matrix))  # Agora deve exibir ENSP IDs

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

### PCA ANTES DO TRANSCRIPTOGRAMER

# Transpor a matriz: agora células (colunas) ficam como linhas
t_copy_matrix <- t(copy_matrix)

# Realizar a PCA
pca_df_before_t <- prcomp(t_copy_matrix, center = TRUE, scale. = FALSE)

# Criar um dataframe com os scores das PCs
pca_df_before_t_scores <- as.data.frame(pca_df_before_t$x)
pca_df_before_t_scores$Condition <- sub("_.*", "", rownames(pca_df_before_t_scores))  # Extrai a condição do nome da célula

# Plot da PCA
ggplot(pca_df_before_t_scores, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "PCA das Células", x = "PC1", y = "PC2")


######
## TRANSCRIPTOGRAMAS 
######

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
  "notreated-batch1" = "#E41A1C",  # Vermelho
  "notreated-batch2" = "#377EB8",  # Azul
  "TGFbeta1-1day-batch2" = "#4DAF4A",  # Verde
  "TGFbeta1-2day-batch2" = "#984EA3",  # Roxo
  "TGFbeta1-3day-batch2" = "#FF7F00",  # Laranja
  "TGFbeta1-4day-batch1" = "#FFFF33",  # Amarelo
  "TGFbeta1-8day-batch1" = "#A65628"   # Marrom
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



