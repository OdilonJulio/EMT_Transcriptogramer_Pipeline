
# Instalar e carregar pacotes necessários

library(transcriptogramer)
library(rsparse)
library(Matrix)
library(biomaRt)
library(Ropj)
library(vroom)
library(tidyverse)
library(ggplot2)
library(parallel)
library(ComplexHeatmap)
library(magrittr)
library(patchwork)

# # Detectar o número total de núcleos disponíveis no sistema
detectCores()


# Para salvar vários objetos com nomes similares no ambiente do R, você pode usar uma função que busca e salva cada objeto em um arquivo separado com base em um padrão comum no nome. Por exemplo, se você tiver várias variáveis no ambiente com nomes começando com "dataset_" e quiser salvar cada uma delas em um arquivo .RData, aqui está como fazer:

salvar_objetos_rdata <- function(...) {
  # Captura os nomes dos objetos passados como argumentos
  objetos <- list(...)
  
  # Loop para salvar cada objeto no arquivo com seu próprio nome
  for (obj_nome in objetos) {
    save(list = obj_nome, file = paste0(obj_nome, ".RData"))
  }
  
  message("Objeto(s) ", objetos," salvo(s) com sucesso.")
}


############################## Após reunião de 27/11/2024 com a Prof Rita ######
# Não usar SEURAT.
# Ler arquivos manualmente.
# Controle de qualidade (MT e nFeatures)
# Normalização manual.
# Junção dos notreateds com o número de colunas pela de menor quantidade, excluindo o restante.
# Concatenação de todas os dias.
# Transcriptogramer de todos os dias.
# Quatro PCAs:
#   Scale T para R=0 e R=30
# Scale F para R=0 e R=30

matrix_notreated_batch1 <- readMM("mtx_conversions/notreated_batch1/outs/filtered_feature_bc_matrix/matrix.mtx")
matrix_notreated_batch2 <- readMM("mtx_conversions/notreated_batch2/outs/filtered_feature_bc_matrix/matrix.mtx")
matrix_TGFbeta1_1day_batch2 <- readMM("mtx_conversions/TGFbeta1_1day_batch2/outs/filtered_feature_bc_matrix/matrix.mtx")
matrix_TGFbeta1_2day_batch2 <- readMM("mtx_conversions/TGFbeta1_2day_batch2/outs/filtered_feature_bc_matrix/matrix.mtx")
matrix_TGFbeta1_3day_batch2 <- readMM("mtx_conversions/TGFbeta1_3day_batch2/outs/filtered_feature_bc_matrix/matrix.mtx")
matrix_TGFbeta1_4day_batch1 <- readMM("mtx_conversions/TGFbeta1_4day_batch1/outs/filtered_feature_bc_matrix/matrix.mtx")
matrix_TGFbeta1_8day_batch1 <- readMM("mtx_conversions/TGFbeta1_8day_batch1/outs/filtered_feature_bc_matrix/matrix.mtx")

genes_notreated_batch1 <- read.table("mtx_conversions/notreated_batch1/outs/filtered_feature_bc_matrix/genes.tsv", header = FALSE, sep = "\t")$V1
genes_notreated_batch2 <- read.table("mtx_conversions/notreated_batch2/outs/filtered_feature_bc_matrix/genes.tsv", header = FALSE, sep = "\t")$V1
genes_TGFbeta1_1day_batch2 <- read.table("mtx_conversions/TGFbeta1_1day_batch2/outs/filtered_feature_bc_matrix/genes.tsv", header = FALSE, sep = "\t")$V1
genes_TGFbeta1_2day_batch2 <- read.table("mtx_conversions/TGFbeta1_2day_batch2/outs/filtered_feature_bc_matrix/genes.tsv", header = FALSE, sep = "\t")$V1
genes_TGFbeta1_3day_batch2 <- read.table("mtx_conversions/TGFbeta1_3day_batch2/outs/filtered_feature_bc_matrix/genes.tsv", header = FALSE, sep = "\t")$V1
genes_TGFbeta1_4day_batch1 <- read.table("mtx_conversions/TGFbeta1_4day_batch1/outs/filtered_feature_bc_matrix/genes.tsv", header = FALSE, sep = "\t")$V1
genes_TGFbeta1_8day_batch1 <- read.table("mtx_conversions/TGFbeta1_8day_batch1/outs/filtered_feature_bc_matrix/genes.tsv", header = FALSE, sep = "\t")$V1
# 
# 
barcodes_notreated_batch1 <- read.table("mtx_conversions/notreated_batch1/outs/filtered_feature_bc_matrix/barcodes.tsv", header = FALSE, sep = "\t")$V1
barcodes_notreated_batch2 <- read.table("mtx_conversions/notreated_batch2/outs/filtered_feature_bc_matrix/barcodes.tsv", header = FALSE, sep = "\t")$V1
barcodes_TGFbeta1_1day_batch2 <- read.table("mtx_conversions/TGFbeta1_1day_batch2/outs/filtered_feature_bc_matrix/barcodes.tsv", header = FALSE, sep = "\t")$V1
barcodes_TGFbeta1_2day_batch2 <- read.table("mtx_conversions/TGFbeta1_2day_batch2/outs/filtered_feature_bc_matrix/barcodes.tsv", header = FALSE, sep = "\t")$V1
barcodes_TGFbeta1_3day_batch2 <- read.table("mtx_conversions/TGFbeta1_3day_batch2/outs/filtered_feature_bc_matrix/barcodes.tsv", header = FALSE, sep = "\t")$V1
barcodes_TGFbeta1_4day_batch1 <- read.table("mtx_conversions/TGFbeta1_4day_batch1/outs/filtered_feature_bc_matrix/barcodes.tsv", header = FALSE, sep = "\t")$V1
barcodes_TGFbeta1_8day_batch1 <- read.table("mtx_conversions/TGFbeta1_8day_batch1/outs/filtered_feature_bc_matrix/barcodes.tsv", header = FALSE, sep = "\t")$V1
#             
# 
rownames(matrix_notreated_batch1) <- genes_notreated_batch1
colnames(matrix_notreated_batch1) <- barcodes_notreated_batch1

rownames(matrix_notreated_batch2) <- genes_notreated_batch2
colnames(matrix_notreated_batch2) <- barcodes_notreated_batch2

rownames(matrix_TGFbeta1_1day_batch2) <- genes_TGFbeta1_1day_batch2
colnames(matrix_TGFbeta1_1day_batch2) <- barcodes_TGFbeta1_1day_batch2

rownames(matrix_TGFbeta1_2day_batch2) <- genes_TGFbeta1_2day_batch2
colnames(matrix_TGFbeta1_2day_batch2) <- barcodes_TGFbeta1_2day_batch2

rownames(matrix_TGFbeta1_3day_batch2) <- genes_TGFbeta1_3day_batch2
colnames(matrix_TGFbeta1_3day_batch2) <- barcodes_TGFbeta1_3day_batch2

rownames(matrix_TGFbeta1_4day_batch1) <- genes_TGFbeta1_4day_batch1
colnames(matrix_TGFbeta1_4day_batch1) <- barcodes_TGFbeta1_4day_batch1

rownames(matrix_TGFbeta1_8day_batch1) <- genes_TGFbeta1_8day_batch1
colnames(matrix_TGFbeta1_8day_batch1) <- barcodes_TGFbeta1_8day_batch1

matrix_notreated_batch1 <- as.matrix(matrix_notreated_batch1)
matrix_notreated_batch2 <- as.matrix(matrix_notreated_batch2)
matrix_TGFbeta1_1day_batch2 <- as.matrix(matrix_TGFbeta1_1day_batch2)
matrix_TGFbeta1_2day_batch2 <- as.matrix(matrix_TGFbeta1_2day_batch2)
matrix_TGFbeta1_3day_batch2 <- as.matrix(matrix_TGFbeta1_3day_batch2)
matrix_TGFbeta1_4day_batch1 <- as.matrix(matrix_TGFbeta1_4day_batch1)
matrix_TGFbeta1_8day_batch1 <- as.matrix(matrix_TGFbeta1_8day_batch1)

head(rownames(matrix_notreated_batch1))

library(biomaRt)

# Conectando ao Ensembl
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")

# Obter genes mitocondriais
mito_genes_info <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
  filters = "chromosome_name", ## DESCOBRIR EXATAMENTE
  values = "MT",  # Mitocôndria
  mart = mart
)

# Vetor de genes mitocondriais
mito_genes <- mito_genes_info$ensembl_gene_id

mito_genes_present_rownames_notreated_batch1 <- rownames(matrix_notreated_batch1) %in% mito_genes
mito_genes_present_rownames_notreated_batch2 <- rownames(matrix_notreated_batch2) %in% mito_genes
mito_genes_present_rownames_TGFbeta1_1day_batch2 <- rownames(matrix_TGFbeta1_1day_batch2) %in% mito_genes
mito_genes_present_rownames_TGFbeta1_2day_batch2 <- rownames(matrix_TGFbeta1_2day_batch2) %in% mito_genes
mito_genes_present_rownames_TGFbeta1_3day_batch2 <- rownames(matrix_TGFbeta1_3day_batch2) %in% mito_genes
mito_genes_present_rownames_TGFbeta1_4day_batch1 <- rownames(matrix_TGFbeta1_4day_batch1) %in% mito_genes
mito_genes_present_rownames_TGFbeta1_8day_batch1 <- rownames(matrix_TGFbeta1_8day_batch1) %in% mito_genes

# Subconjunto apenas com genes mitocondriais
mito_genes_present_matrix_notreated_batch1 <- matrix_notreated_batch1[mito_genes_present_rownames_notreated_batch1, ]
mito_genes_present_matrix_notreated_batch2 <- matrix_notreated_batch2[mito_genes_present_rownames_notreated_batch2, ]
mito_genes_present_matrix_TGFbeta1_1day_batch2 <- matrix_TGFbeta1_1day_batch2[mito_genes_present_rownames_TGFbeta1_1day_batch2, ]
mito_genes_present_matrix_TGFbeta1_2day_batch2 <- matrix_TGFbeta1_2day_batch2[mito_genes_present_rownames_TGFbeta1_2day_batch2, ]
mito_genes_present_matrix_TGFbeta1_3day_batch2 <- matrix_TGFbeta1_3day_batch2[mito_genes_present_rownames_TGFbeta1_3day_batch2, ]
mito_genes_present_matrix_TGFbeta1_4day_batch1 <- matrix_TGFbeta1_4day_batch1[mito_genes_present_rownames_TGFbeta1_4day_batch1, ]
mito_genes_present_matrix_TGFbeta1_8day_batch1 <- matrix_TGFbeta1_8day_batch1[mito_genes_present_rownames_TGFbeta1_8day_batch1, ]

# Total de expressão por célula
total_expression_matrix_notreated_batch1 <- colSums(matrix_notreated_batch1)
total_expression_matrix_notreated_batch2 <- colSums(matrix_notreated_batch2)
total_expression_matrix_TGFbeta1_1day_batch2 <- colSums(matrix_TGFbeta1_1day_batch2)
total_expression_matrix_TGFbeta1_2day_batch2 <- colSums(matrix_TGFbeta1_2day_batch2)
total_expression_matrix_TGFbeta1_3day_batch2 <- colSums(matrix_TGFbeta1_3day_batch2)
total_expression_matrix_TGFbeta1_4day_batch1 <- colSums(matrix_TGFbeta1_4day_batch1)
total_expression_matrix_TGFbeta1_8day_batch1 <- colSums(matrix_TGFbeta1_8day_batch1)

# Expressão mitocondrial por célula
mito_expression_per_cell_matrix_notreated_batch1 <- colSums(matrix_notreated_batch1[mito_genes_present_rownames_notreated_batch1, ])
mito_expression_per_cell_matrix_notreated_batch2 <- colSums(matrix_notreated_batch2[mito_genes_present_rownames_notreated_batch2, ])
mito_expression_per_cell_matrix_TGFbeta1_1day_batch2 <- colSums(matrix_TGFbeta1_1day_batch2[mito_genes_present_rownames_TGFbeta1_1day_batch2, ])
mito_expression_per_cell_matrix_TGFbeta1_2day_batch2 <- colSums(matrix_TGFbeta1_2day_batch2[mito_genes_present_rownames_TGFbeta1_2day_batch2, ])
mito_expression_per_cell_matrix_TGFbeta1_3day_batch2 <- colSums(matrix_TGFbeta1_3day_batch2[mito_genes_present_rownames_TGFbeta1_3day_batch2, ])
mito_expression_per_cell_matrix_TGFbeta1_4day_batch1 <- colSums(matrix_TGFbeta1_4day_batch1[mito_genes_present_rownames_TGFbeta1_4day_batch1, ])
mito_expression_per_cell_matrix_TGFbeta1_8day_batch1 <- colSums(matrix_TGFbeta1_8day_batch1[mito_genes_present_rownames_TGFbeta1_8day_batch1, ])

# Porcentagem de expressão mitocondrial
mito_percentage_matrix_notreated_batch1 <- (mito_expression_per_cell_matrix_notreated_batch1 / total_expression_matrix_notreated_batch1) * 100
mito_percentage_matrix_notreated_batch2 <- (mito_expression_per_cell_matrix_notreated_batch2 / total_expression_matrix_notreated_batch2) * 100
mito_percentage_matrix_TGFbeta1_1day_batch2 <- (mito_expression_per_cell_matrix_TGFbeta1_1day_batch2 / total_expression_matrix_TGFbeta1_1day_batch2) * 100
mito_percentage_matrix_TGFbeta1_2day_batch2 <- (mito_expression_per_cell_matrix_TGFbeta1_2day_batch2 / total_expression_matrix_TGFbeta1_2day_batch2) * 100
mito_percentage_matrix_TGFbeta1_3day_batch2 <- (mito_expression_per_cell_matrix_TGFbeta1_3day_batch2 / total_expression_matrix_TGFbeta1_3day_batch2) * 100
mito_percentage_matrix_TGFbeta1_4day_batch1 <- (mito_expression_per_cell_matrix_TGFbeta1_4day_batch1 / total_expression_matrix_TGFbeta1_4day_batch1) * 100
mito_percentage_matrix_TGFbeta1_8day_batch1 <- (mito_expression_per_cell_matrix_TGFbeta1_8day_batch1 / total_expression_matrix_TGFbeta1_8day_batch1) * 100

cels_high_mito_percentage_matrix_notreated_batch1 <- mito_percentage_matrix_notreated_batch1 > 20
cels_high_mito_percentage_matrix_notreated_batch2 <- mito_percentage_matrix_notreated_batch2 > 20
cels_high_mito_percentage_matrix_TGFbeta1_1day_batch2 <- mito_percentage_matrix_TGFbeta1_1day_batch2 > 20
cels_high_mito_percentage_matrix_TGFbeta1_2day_batch2 <- mito_percentage_matrix_TGFbeta1_2day_batch2 > 20
cels_high_mito_percentage_matrix_TGFbeta1_3day_batch2 <- mito_percentage_matrix_TGFbeta1_3day_batch2 > 20
cels_high_mito_percentage_matrix_TGFbeta1_4day_batch1 <- mito_percentage_matrix_TGFbeta1_4day_batch1 > 20
cels_high_mito_percentage_matrix_TGFbeta1_8day_batch1 <- mito_percentage_matrix_TGFbeta1_8day_batch1 > 20

# Subconjunto da matriz de expressão sem as células com alta expressão mitocondrial
filtered_matrix_notreated_batch1 <- matrix_notreated_batch1[, !cels_high_mito_percentage_matrix_notreated_batch1]
filtered_matrix_notreated_batch2 <- matrix_notreated_batch2[, !cels_high_mito_percentage_matrix_notreated_batch2]
filtered_matrix_TGFbeta1_1day_batch2 <- matrix_TGFbeta1_1day_batch2[, !cels_high_mito_percentage_matrix_TGFbeta1_1day_batch2]
filtered_matrix_TGFbeta1_2day_batch2 <- matrix_TGFbeta1_2day_batch2[, !cels_high_mito_percentage_matrix_TGFbeta1_2day_batch2]
filtered_matrix_TGFbeta1_3day_batch2 <- matrix_TGFbeta1_3day_batch2[, !cels_high_mito_percentage_matrix_TGFbeta1_3day_batch2]
filtered_matrix_TGFbeta1_4day_batch1 <- matrix_TGFbeta1_4day_batch1[, !cels_high_mito_percentage_matrix_TGFbeta1_4day_batch1]
filtered_matrix_TGFbeta1_8day_batch1 <- matrix_TGFbeta1_8day_batch1[, !cels_high_mito_percentage_matrix_TGFbeta1_8day_batch1]


## Agora, faremos o QC para nFeatures acima de 500.

# Calcular nFeatures (genes expressos por célula)
nFeatures_matrix_notreated_batch1 <- colSums(filtered_matrix_notreated_batch1 > 0)
nFeatures_matrix_notreated_batch2 <- colSums(filtered_matrix_notreated_batch2 > 0)
nFeatures_matrix_TGFbeta1_1day_batch2 <- colSums(filtered_matrix_TGFbeta1_1day_batch2 > 0)
nFeatures_matrix_TGFbeta1_2day_batch2 <- colSums(filtered_matrix_TGFbeta1_2day_batch2 > 0)
nFeatures_matrix_TGFbeta1_3day_batch2 <- colSums(filtered_matrix_TGFbeta1_3day_batch2 > 0)
nFeatures_matrix_TGFbeta1_4day_batch1 <- colSums(filtered_matrix_TGFbeta1_4day_batch1 > 0)
nFeatures_matrix_TGFbeta1_8day_batch1 <- colSums(filtered_matrix_TGFbeta1_8day_batch1 > 0)

# Identificar células com mais de 500 genes expressos
cells_with_high_features_matrix_notreated_batch1 <- nFeatures_matrix_notreated_batch1 > 500
cells_with_high_features_matrix_notreated_batch2 <- nFeatures_matrix_notreated_batch2 > 500
cells_with_high_features_matrix_TGFbeta1_1day_batch2 <- nFeatures_matrix_TGFbeta1_1day_batch2 > 500
cells_with_high_features_matrix_TGFbeta1_2day_batch2 <- nFeatures_matrix_TGFbeta1_2day_batch2 > 500
cells_with_high_features_matrix_TGFbeta1_3day_batch2 <- nFeatures_matrix_TGFbeta1_3day_batch2 > 500
cells_with_high_features_matrix_TGFbeta1_4day_batch1 <- nFeatures_matrix_TGFbeta1_4day_batch1 > 500
cells_with_high_features_matrix_TGFbeta1_8day_batch1 <- nFeatures_matrix_TGFbeta1_8day_batch1 > 500

# Filtrar a matriz para manter somente essas células
filtered_matrix_notreated_batch1 <- filtered_matrix_notreated_batch1[, cells_with_high_features_matrix_notreated_batch1]
filtered_matrix_notreated_batch2 <- filtered_matrix_notreated_batch2[, cells_with_high_features_matrix_notreated_batch2]
filtered_matrix_TGFbeta1_1day_batch2 <- filtered_matrix_TGFbeta1_1day_batch2[, cells_with_high_features_matrix_TGFbeta1_1day_batch2]
filtered_matrix_TGFbeta1_2day_batch2 <- filtered_matrix_TGFbeta1_2day_batch2[, cells_with_high_features_matrix_TGFbeta1_2day_batch2]
filtered_matrix_TGFbeta1_3day_batch2 <- filtered_matrix_TGFbeta1_3day_batch2[, cells_with_high_features_matrix_TGFbeta1_3day_batch2]
filtered_matrix_TGFbeta1_4day_batch1 <- filtered_matrix_TGFbeta1_4day_batch1[, cells_with_high_features_matrix_TGFbeta1_4day_batch1]
filtered_matrix_TGFbeta1_8day_batch1 <- filtered_matrix_TGFbeta1_8day_batch1[, cells_with_high_features_matrix_TGFbeta1_8day_batch1]

## Normalização

# Calcular a soma de cada coluna
column_sums_filtered_matrix_notreated_batch1 <- colSums(filtered_matrix_notreated_batch1)
column_sums_filtered_matrix_notreated_batch2 <- colSums(filtered_matrix_notreated_batch2)
column_sums_filtered_matrix_TGFbeta1_1day_batch2 <- colSums(filtered_matrix_TGFbeta1_1day_batch2)
column_sums_filtered_matrix_TGFbeta1_2day_batch2 <- colSums(filtered_matrix_TGFbeta1_2day_batch2)
column_sums_filtered_matrix_TGFbeta1_3day_batch2 <- colSums(filtered_matrix_TGFbeta1_3day_batch2)
column_sums_filtered_matrix_TGFbeta1_4day_batch1 <- colSums(filtered_matrix_TGFbeta1_4day_batch1)
column_sums_filtered_matrix_TGFbeta1_8day_batch1 <- colSums(filtered_matrix_TGFbeta1_8day_batch1)

# Dividir cada valor pela soma da coluna correspondente
normalized_matrix_notreated_batch1 <- sweep(filtered_matrix_notreated_batch1, 2, column_sums_filtered_matrix_notreated_batch1, "/")
normalized_matrix_notreated_batch2 <- sweep(filtered_matrix_notreated_batch2, 2, column_sums_filtered_matrix_notreated_batch2, "/")
normalized_matrix_TGFbeta1_1day_batch2 <- sweep(filtered_matrix_TGFbeta1_1day_batch2, 2, column_sums_filtered_matrix_TGFbeta1_1day_batch2, "/")
normalized_matrix_TGFbeta1_2day_batch2 <- sweep(filtered_matrix_TGFbeta1_2day_batch2, 2, column_sums_filtered_matrix_TGFbeta1_2day_batch2, "/")
normalized_matrix_TGFbeta1_3day_batch2 <- sweep(filtered_matrix_TGFbeta1_3day_batch2, 2, column_sums_filtered_matrix_TGFbeta1_3day_batch2, "/")
normalized_matrix_TGFbeta1_4day_batch1 <- sweep(filtered_matrix_TGFbeta1_4day_batch1, 2, column_sums_filtered_matrix_TGFbeta1_4day_batch1, "/")
normalized_matrix_TGFbeta1_8day_batch1 <- sweep(filtered_matrix_TGFbeta1_8day_batch1, 2, column_sums_filtered_matrix_TGFbeta1_8day_batch1, "/")

# Verificar a soma das colunas
col_sums_normalized_matrix_notreated_batch1 <- colSums(normalized_matrix_notreated_batch1)
col_sums_normalized_matrix_notreated_batch2 <- colSums(normalized_matrix_notreated_batch2)
col_sums_normalized_matrix_TGFbeta1_1day_batch2 <- colSums(normalized_matrix_TGFbeta1_1day_batch2)
col_sums_normalized_matrix_TGFbeta1_2day_batch2 <- colSums(normalized_matrix_TGFbeta1_2day_batch2)
col_sums_normalized_matrix_TGFbeta1_3day_batch2 <- colSums(normalized_matrix_TGFbeta1_3day_batch2)
col_sums_normalized_matrix_TGFbeta1_4day_batch1 <- colSums(normalized_matrix_TGFbeta1_4day_batch1)
col_sums_normalized_matrix_TGFbeta1_8day_batch1 <- colSums(normalized_matrix_TGFbeta1_8day_batch1)

print(col_sums_normalized_matrix_notreated_batch1)  # Deve ser aproximadamente igual a 1 para todas as colunas
print(col_sums_normalized_matrix_notreated_batch2)  # Deve ser aproximadamente igual a 1 para todas as colunas
print(col_sums_normalized_matrix_TGFbeta1_1day_batch2)  # Deve ser aproximadamente igual a 1 para todas as colunas
print(col_sums_normalized_matrix_TGFbeta1_2day_batch2)  # Deve ser aproximadamente igual a 1 para todas as colunas
print(col_sums_normalized_matrix_TGFbeta1_3day_batch2)  # Deve ser aproximadamente igual a 1 para todas as colunas
print(col_sums_normalized_matrix_TGFbeta1_4day_batch1)  # Deve ser aproximadamente igual a 1 para todas as colunas
print(col_sums_normalized_matrix_TGFbeta1_8day_batch1)  # Deve ser aproximadamente igual a 1 para todas as colunas

sum(0.999999 > col_sums_normalized_matrix_notreated_batch1 | 
      col_sums_normalized_matrix_notreated_batch1 > 1.00000001) # Deve ser 0.
sum(0.999999 > col_sums_normalized_matrix_notreated_batch2 | 
      col_sums_normalized_matrix_notreated_batch2 > 1.00000001) # Deve ser 0.
sum(0.999999 > col_sums_normalized_matrix_TGFbeta1_1day_batch2 | 
      col_sums_normalized_matrix_TGFbeta1_1day_batch2 > 1.00000001) # Deve ser 0.
sum(0.999999 > col_sums_normalized_matrix_TGFbeta1_2day_batch2 | 
      col_sums_normalized_matrix_TGFbeta1_2day_batch2 > 1.00000001) # Deve ser 0.
sum(0.999999 > col_sums_normalized_matrix_TGFbeta1_3day_batch2 | 
      col_sums_normalized_matrix_TGFbeta1_3day_batch2 > 1.00000001) # Deve ser 0.
sum(0.999999 > col_sums_normalized_matrix_TGFbeta1_4day_batch1 | 
      col_sums_normalized_matrix_TGFbeta1_4day_batch1 > 1.00000001) # Deve ser 0.
sum(0.999999 > col_sums_normalized_matrix_TGFbeta1_8day_batch1 | 
      col_sums_normalized_matrix_TGFbeta1_8day_batch1 > 1.00000001) # Deve ser 0.


View(normalized_matrix_notreated_batch1)
View(normalized_matrix_notreated_batch2)
View(normalized_matrix_TGFbeta1_1day_batch2)
View(normalized_matrix_TGFbeta1_2day_batch2)
View(normalized_matrix_TGFbeta1_3day_batch2)
View(normalized_matrix_TGFbeta1_4day_batch1)
View(normalized_matrix_TGFbeta1_8day_batch1)

## Vamos transformar os dados Notreateds em apenas uma matriz, realizando a média da expressão.

# Obtenha o número de colunas de cada matriz
n_col_batch1 <- ncol(normalized_matrix_notreated_batch1)
n_col_batch2 <- ncol(normalized_matrix_notreated_batch2)

# Identifique o menor número de colunas
min_cols <- min(n_col_batch1, n_col_batch2)

# Truncate as matrizes para terem o mesmo número de colunas
matrix_batch1_trimmed <- normalized_matrix_notreated_batch1[, 1:min_cols]
matrix_batch2_trimmed <- normalized_matrix_notreated_batch2[, 1:min_cols]

# Calcule a média elemento a elemento
normalized_matrix_notreated_mean <- (matrix_batch1_trimmed + matrix_batch2_trimmed) / 2

# Verifique a matriz resultante
View(normalized_matrix_notreated_mean)

# Crie uma lista com as matrizes e suas respectivas referências
matrices <- list(
  notreated_mean = normalized_matrix_notreated_mean,
  TGFbeta1_1day_batch2 = normalized_matrix_TGFbeta1_1day_batch2,
  TGFbeta1_2day_batch2 = normalized_matrix_TGFbeta1_2day_batch2,
  TGFbeta1_3day_batch2 = normalized_matrix_TGFbeta1_3day_batch2,
  TGFbeta1_4day_batch1 = normalized_matrix_TGFbeta1_4day_batch1,
  TGFbeta1_8day_batch1 = normalized_matrix_TGFbeta1_8day_batch1
)

# Adicione o sufixo aos nomes das colunas de cada matriz
renamed_matrices <- lapply(names(matrices), function(name) {
  mat <- matrices[[name]]
  colnames(mat) <- paste0(colnames(mat), "_", name) # Adiciona o sufixo com o nome da origem
  return(mat)
})

# Combine as matrizes lado a lado
combined_matrix <- do.call(cbind, renamed_matrices)

# Verifique a matriz resultante
dim(combined_matrix) # Mostra dimensões para conferir
head(colnames(combined_matrix)) # Mostra os novos nomes das colunas


# Transcriptogramas

library(biomaRt)

# Lista os biomarts disponíveis no espelho dos EUA
listMarts(host = "https://us.ensembl.org")

# Espelho principal (www)
listMarts(host = "https://www.ensembl.org")

# Espelho da Ásia
listMarts(host = "https://asia.ensembl.org")

# Conectar ao biomart de genes humanos 
# Mapear os nomes dos genes para ENSEMBL
gene_mapping <- getBM(
  attributes = c("external_gene_name", "ensembl_gene_id", "ensembl_peptide_id"),
  filters = "external_gene_name",
  values = rownames(combined_matrix),
  mart = mart
)

# # Mapear os nomes dos genes na matriz para `ensembl_gene_id`
mapped_genes <- gene_mapping$ensembl_gene_id[match(
  rownames(combined_matrix),
  gene_mapping$external_gene_name
)]

dictionary <- getBM(attributes = c("ensembl_peptide_id", "ensembl_gene_id"),
                    mart = mart)
dictionary %>%
  mutate(ensembl_peptide_id = ifelse(ensembl_peptide_id == "", NA, ensembl_peptide_id)) %>%
  na.omit() -> dictionary

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

##########


############ PCA ANTES DO TRANSCRIPTOGRAMER - A PETIDO DE RAFAELA 05.12.2024

# Obter os identificadores ENSP correspondentes aos ENSG
genes_mapping <- getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id"),
                       filters = "ensembl_gene_id",
                       values = rownames(combined_matrix),
                       mart = mart)

# Remover as linhas com ENSP vazio
genes_mapping <- genes_mapping[!is.na(genes_mapping$ensembl_peptide_id), ]

# Substituir ENSG por ENSP nos nomes das linhas da matriz
rownames(combined_matrix) <- genes_mapping$ensembl_peptide_id[match(rownames(combined_matrix), genes_mapping$ensembl_gene_id)]

# Excluir linhas sem ENSP
combined_matrix <- combined_matrix[!is.na(rownames(combined_matrix)), ]

# Verifique as primeiras linhas da matriz
head(combined_matrix)


# Extraindo o dataframe do slot transcriptogramS2
pca_Rafa_combined_matrix <- combined_matrix

# Mantendo apenas colunas numéricas
pca_Rafa_combined_matrix <- pca_Rafa_combined_matrix[, sapply(pca_Rafa_combined_matrix, is.numeric)]

# Verificando se sobrou algum dado
if (ncol(pca_Rafa_combined_matrix) == 0) {
  stop("Nenhuma coluna numérica encontrada no dataframe!")
}

# Realizando a padronização dos dados numéricos
pca_Rafa_combined_matrix <- scale(pca_Rafa_combined_matrix)

# Continuando com o PCA
pca_Rafa_combined_matrix <- prcomp(pca_Rafa_combined_matrix, center = TRUE, scale. = TRUE)


# Scree plot da variância explicada
plot(pca_Rafa_combined_matrix, type = "l", main = "Scree Plot - PCA - Matriz Combinada - ANTES DO TRANSCRIPTOGRAMER")

# Plotando os componentes principais (PC1 vs. PC2)
library(ggplot2)
pca_data <- data.frame(pca_Rafa_combined_matrix$x)
ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point() +
  labs(title = "PCA - PC1 vs PC2", x = "PC1", y = "PC2") +
  theme_minimal()

##############



# radius = 0

t_combined_matrix_R0 <- transcriptogramPreprocess(association = assoc, 
                                                   ordering = ord$Protein, 
                                                   radius = 0)

t_combined_matrix_R0 <- transcriptogramStep1(object = t_combined_matrix_R0, 
                                              expression = combined_matrix, 
                                              dictionary = dictionary)
t_combined_matrix_R0 <- transcriptogramStep2(object = t_combined_matrix_R0)

salvar_objetos_rdata("t_combined_matrix_R0")

# radius = 30

t_combined_matrix_R30 <- transcriptogramPreprocess(association = assoc, 
                                                   ordering = ord$Protein, 
                                                   radius = 30)

t_combined_matrix_R30 <- transcriptogramStep1(object = t_combined_matrix_R30, 
                                              expression = combined_matrix, 
                                              dictionary = dictionary)
t_combined_matrix_R30 <- transcriptogramStep2(object = t_combined_matrix_R30)


salvar_objetos_rdata("t_combined_matrix_R30")

#### Fazendo a PCA de t_combined_matrix_R0 - Scaled = TRUE ###

# Extraindo o dataframe do slot transcriptogramS2
df_t_combined_matrix_R0 <- t_combined_matrix_R0@transcriptogramS2[,-2]

# Verificando o tipo de cada coluna
sapply(df_t_combined_matrix_R0, class)  # Para identificar as colunas que não são numéricas

# Mantendo apenas colunas numéricas
df_t_combined_matrix_R0 <- df_t_combined_matrix_R0[, sapply(df_t_combined_matrix_R0, is.numeric)]

# Verificando se sobrou algum dado
if (ncol(df_t_combined_matrix_R0) == 0) {
  stop("Nenhuma coluna numérica encontrada no dataframe!")
}

# Realizando a padronização dos dados numéricos
df_scaled_t_combined_matrix_R0 <- scale(df_t_combined_matrix_R0)

# Continuando com o PCA
pca_t_combined_matrix_R0 <- prcomp(df_scaled_t_combined_matrix_R0, center = TRUE, scale. = TRUE)


summary(pca_t_combined_matrix_R0)

# Scree plot da variância explicada
plot(pca_t_combined_matrix_R0, type = "l", main = "Scree Plot - PCA - Matriz Combinada - R = 0")

## Fazendo a PCA de t_combined_matrix_R30

## Fazendo a PCA de t_combined_matrix_R30 - Scaled = TRUE

# Extraindo o dataframe do slot transcriptogramS2
df_t_combined_matrix_R30 <- t_combined_matrix_R30@transcriptogramS2[,-2]

# Verificando o tipo de cada coluna
sapply(df_t_combined_matrix_R30, class)  # Para identificar as colunas que não são numéricas

# Mantendo apenas colunas numéricas
df_t_combined_matrix_R30 <- df_t_combined_matrix_R30[, sapply(df_t_combined_matrix_R30, is.numeric)]

# Verificando se sobrou algum dado
if (ncol(df_t_combined_matrix_R30) == 0) {
  stop("Nenhuma coluna numérica encontrada no dataframe!")
}

# Realizando a padronização dos dados numéricos
df_scaled_t_combined_matrix_R30 <- scale(df_t_combined_matrix_R30)

# Continuando com o PCA
pca_t_combined_matrix_R30 <- prcomp(df_scaled_t_combined_matrix_R30, center = TRUE, scale. = TRUE)


summary(pca_t_combined_matrix_R30)

# Scree plot da variância explicada
plot(pca_t_combined_matrix_R30, type = "l", main = "Scree Plot - PCA - Matriz Combinada - R = 30")

library(ggplot2)

# Variância explicada para o PCA radius = 0
explained_variance_combined_matrix_R0 <- pca_t_combined_matrix_R0$sdev^2 / sum(pca_t_combined_matrix_R0$sdev^2)

# Variância explicada para o PCA radius = 30
explained_variance_combined_matrix_R30 <- pca_t_combined_matrix_R30$sdev^2 / sum(pca_t_combined_matrix_R30$sdev^2)

# Verificar o número de PCs em cada PCA
length(explained_variance_combined_matrix_R0)  # Quantidade de PCs no PCA radius = 0
length(explained_variance_combined_matrix_R30)  # Quantidade de PCs no PCA radius = 30

# Ajustar o comprimento de ambos os vetores de variância explicada (preenchendo com NAs ou truncando)
# Vamos supor que você queira alinhar os PCs até o número máximo de componentes entre os dois batches.

max_pcs <- max(length(explained_variance_combined_matrix_R0), length(explained_variance_combined_matrix_R30))

# Preencher com NAs para garantir que ambos tenham o mesmo número de componentes
explained_variance_combined_matrix_R0 <- c(explained_variance_combined_matrix_R0, rep(NA, max_pcs - length(explained_variance_combined_matrix_R0)))
explained_variance_combined_matrix_R30 <- c(explained_variance_combined_matrix_R30, rep(NA, max_pcs - length(explained_variance_combined_matrix_R30)))

# Criar o dataframe novamente com os dados ajustados
scree_data <- data.frame(
  PC = rep(1:max_pcs, 2),
  Variance = c(explained_variance_combined_matrix_R0, explained_variance_combined_matrix_R30),
  Radius = rep(c("Radius = 0", "Radius = 30"), each = max_pcs)
)

# Variância cumulativa
scree_data$cumulative_variance <- ave(scree_data$Variance, scree_data$Radius, FUN = cumsum)

# Plot do Scree Plot com variância cumulativa
png("images/PCAeVarianciaAcumuladaTranscriptograma_Combined_Matrix_R0_R30.png", width = 1920, height = 1080, res = 200)
ggplot(scree_data, aes(x = PC, y = Variance, color = Radius)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Scree Plot com Variância Explicada - PCA - Scaled = TRUE",
    x = "Principal Component",
    y = "Explained Variance"
  ) +
  theme_minimal() +
  theme(legend.position = "top") +
  geom_line(aes(x = PC, y = cumulative_variance, color = Radius), linetype = "dashed") +
  scale_color_manual(values = c("blue", "red")) +
  scale_x_continuous(breaks = 1:max(scree_data$PC))
dev.off()

# Filtrando para mostrar apenas as TOP 10 PCs
scree_data_top10 <- scree_data[scree_data$PC <= 10, ]

# Plot do Scree Plot com variância cumulativa para as TOP 10 PCs
# Plot do Scree Plot com variância cumulativa

png("images/PCAeVarianciaAcumuladaTop10Transcriptograma_Combined_Matrix_R0_R30.png", width = 1920, height = 1080, res = 200)
ggplot(scree_data_top10, aes(x = PC, y = Variance, color = Radius)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Scree Plot com Variância Explicada - TOP 10 PCs - Scaled = TRUE",
    x = "Principal Component",
    y = "Explained Variance"
  ) +
  theme_minimal() +
  theme(legend.position = "top") +
  geom_line(aes(x = PC, y = cumulative_variance, color = Radius), linetype = "dashed") +
  scale_color_manual(values = c("blue", "red")) +
  scale_x_continuous(breaks = 1:10)  # Exibe apenas os 10 primeiros PCs no eixo X
dev.off()


## Repetindo o processo de PCA com Scaled = FALSE


## Fazendo a PCA de t_combined_matrix_R0 - Scaled = FALSE

# Extraindo o dataframe do slot transcriptogramS2
df_t_combined_matrix_R0_ScaledFALSE <- t_combined_matrix_R0@transcriptogramS2[,-2]

# Verificando o tipo de cada coluna
sapply(df_t_combined_matrix_R0_ScaledFALSE, class)  # Para identificar as colunas que não são numéricas

# Mantendo apenas colunas numéricas
df_t_combined_matrix_R0_ScaledFALSE <- df_t_combined_matrix_R0_ScaledFALSE[, sapply(df_t_combined_matrix_R0_ScaledFALSE, is.numeric)]

# Verificando se sobrou algum dado
if (ncol(df_t_combined_matrix_R0_ScaledFALSE) == 0) {
  stop("Nenhuma coluna numérica encontrada no dataframe!")
}

# Realizando a padronização dos dados numéricos
# df_scaled_t_combined_matrix_R0_ScaledFALSE <- scale(df_t_combined_matrix_R0_ScaledFALSE)

# Continuando com o PCA
pca_t_combined_matrix_R0_ScaledFALSE <- prcomp(df_t_combined_matrix_R0_ScaledFALSE, center = TRUE, scale. = FALSE)


summary(pca_t_combined_matrix_R0_ScaledFALSE)

# Scree plot da variância explicada
plot(pca_t_combined_matrix_R0_ScaledFALSE, type = "l", main = "Scree Plot - PCA - Matriz Combinada - R = 0 - Scaled = FALSE")


## Fazendo a PCA de t_combined_matrix_R30 - Scaled = TRUE

# Extraindo o dataframe do slot transcriptogramS2
df_t_combined_matrix_R30_ScaledFALSE <- t_combined_matrix_R30@transcriptogramS2[,-2]

# Verificando o tipo de cada coluna
sapply(df_t_combined_matrix_R30_ScaledFALSE, class)  # Para identificar as colunas que não são numéricas

# Mantendo apenas colunas numéricas
df_t_combined_matrix_R30_ScaledFALSE <- df_t_combined_matrix_R30_ScaledFALSE[, sapply(df_t_combined_matrix_R30_ScaledFALSE, is.numeric)]

# Verificando se sobrou algum dado
if (ncol(df_t_combined_matrix_R30_ScaledFALSE) == 0) {
  stop("Nenhuma coluna numérica encontrada no dataframe!")
}

# Realizando a padronização dos dados numéricos
# df_scaled_t_combined_matrix_R30 <- scale(df_t_combined_matrix_R30)

# Continuando com o PCA
pca_t_combined_matrix_R30_ScaledFALSE <- prcomp(df_t_combined_matrix_R30_ScaledFALSE, center = TRUE, scale. = FALSE)


summary(pca_t_combined_matrix_R30_ScaledFALSE)

# Scree plot da variância explicada
plot(pca_t_combined_matrix_R30_ScaledFALSE, type = "l", main = "Scree Plot - PCA - Matriz Combinada - R = 30 - Scaled = FALSE")

library(ggplot2)

# Variância explicada para o PCA radius = 0
explained_variance_combined_matrix_R0_ScaledFALSE <- pca_t_combined_matrix_R0_ScaledFALSE$sdev^2 / sum(pca_t_combined_matrix_R0_ScaledFALSE$sdev^2)

# Variância explicada para o PCA radius = 30
explained_variance_combined_matrix_R30_ScaledFALSE <- pca_t_combined_matrix_R30_ScaledFALSE$sdev^2 / sum(pca_t_combined_matrix_R30_ScaledFALSE$sdev^2)

# Verificar o número de PCs em cada PCA
length(explained_variance_combined_matrix_R0_ScaledFALSE)  # Quantidade de PCs no PCA radius = 0
length(explained_variance_combined_matrix_R30_ScaledFALSE)  # Quantidade de PCs no PCA radius = 30

# Ajustar o comprimento de ambos os vetores de variância explicada (preenchendo com NAs ou truncando)
# Vamos supor que você queira alinhar os PCs até o número máximo de componentes entre os dois batches.

max_pcs_ScaledFALSE <- max(length(explained_variance_combined_matrix_R0_ScaledFALSE), length(explained_variance_combined_matrix_R30_ScaledFALSE))

# Preencher com NAs para garantir que ambos tenham o mesmo número de componentes
explained_variance_combined_matrix_R0_ScaledFALSE <- c(explained_variance_combined_matrix_R0_ScaledFALSE, rep(NA, max_pcs_ScaledFALSE - length(explained_variance_combined_matrix_R0_ScaledFALSE)))
explained_variance_combined_matrix_R30_ScaledFALSE <- c(explained_variance_combined_matrix_R30_ScaledFALSE, rep(NA, max_pcs_ScaledFALSE - length(explained_variance_combined_matrix_R30_ScaledFALSE)))

# Criar o dataframe novamente com os dados ajustados
scree_data_ScaledFALSE <- data.frame(
  PC = rep(1:max_pcs_ScaledFALSE, 2),
  Variance = c(explained_variance_combined_matrix_R0_ScaledFALSE, explained_variance_combined_matrix_R30_ScaledFALSE),
  Radius = rep(c("Radius = 0", "Radius = 30"), each = max_pcs_ScaledFALSE)
)

# Variância cumulativa  
scree_data_ScaledFALSE$cumulative_variance <- ave(scree_data_ScaledFALSE$Variance, scree_data_ScaledFALSE$Radius, FUN = cumsum)

# Plot do Scree Plot com variância cumulativa
png("images/PCAeVarianciaAcumuladaTranscriptograma_Combined_Matrix_R0_R30_ScaledFALSE.png", width = 1920, height = 1080, res = 200)
ggplot(scree_data_ScaledFALSE, aes(x = PC, y = Variance, color = Radius)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Scree Plot com Variância Explicada - PCA - Scaled = FALSE",
    x = "Principal Component",
    y = "Explained Variance"
  ) +
  theme_minimal() +
  theme(legend.position = "top") +
  geom_line(aes(x = PC, y = cumulative_variance, color = Radius), linetype = "dashed") +
  scale_color_manual(values = c("blue", "red")) +
  scale_x_continuous(breaks = 1:max(scree_data_ScaledFALSE$PC))
dev.off()

# Filtrando para mostrar apenas as TOP 10 PCs
scree_data_top10_ScaledFALSE <- scree_data_ScaledFALSE[scree_data_ScaledFALSE$PC <= 10, ]

# Plot do Scree Plot com variância cumulativa para as TOP 10 PCs
# Plot do Scree Plot com variância cumulativa

png("images/PCAeVarianciaAcumuladaTop10Transcriptograma_Combined_Matrix_R0_R30_ScaledFALSE.png", width = 1920, height = 1080, res = 200)
ggplot(scree_data_top10_ScaledFALSE, aes(x = PC, y = Variance, color = Radius)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Scree Plot com Variância Explicada - TOP 10 PCs - Scaled = FALSE",
    x = "Principal Component",
    y = "Explained Variance"
  ) +
  theme_minimal() +
  theme(legend.position = "top") +
  geom_line(aes(x = PC, y = cumulative_variance, color = Radius), linetype = "dashed") +
  scale_color_manual(values = c("blue", "red")) +
  scale_x_continuous(breaks = 1:10)  # Exibe apenas os 10 primeiros PCs no eixo X
dev.off()


### Plotando PCA em duas dimensões (PC1 e PC2).


# Coordenadas dos dois primeiros PCs para R0
pc1_pc2_R0 <- as.data.frame(pca_t_combined_matrix_R0$x[, 1:2])
colnames(pc1_pc2_R0) <- c("PC1", "PC2")

# Adicionar informação sobre o grupo (se houver), caso você tenha algum fator de agrupamento, por exemplo, batch ou tipo de amostra
# pc1_pc2_R0$group <- ...

# Plotando PC1 vs PC2
png("images/PCA_PC1_PC2_R0.png", width = 1920, height = 1080, res = 200)
ggplot(pc1_pc2_R0, aes(x = PC1, y = PC2)) +
  geom_point(color = "blue") +  # Você pode adicionar coloração conforme algum critério, como o grupo
  labs(
    title = "Distribuição das Amostras - PCA - R0",
    x = "PC1",
    y = "PC2"
  ) +
  theme_minimal()
dev.off()

# Coordenadas dos dois primeiros PCs para R30
pc1_pc2_R30 <- as.data.frame(pca_t_combined_matrix_R30$x[, 1:2])
colnames(pc1_pc2_R30) <- c("PC1", "PC2")

# Adicionar informação sobre o grupo (se houver), caso você tenha algum fator de agrupamento, por exemplo, batch ou tipo de amostra
# pc1_pc2_R30$group <- ...

# Plotando PC1 vs PC2
png("images/PCA_PC1_PC2_R30.png", width = 1920, height = 1080, res = 200)
ggplot(pc1_pc2_R30, aes(x = PC1, y = PC2)) +
  geom_point(color = "red") +  # Você pode adicionar coloração conforme algum critério, como o grupo
  labs(
    title = "Distribuição das Amostras - PCA - R30",
    x = "PC1",
    y = "PC2"
  ) +
  theme_minimal()
dev.off()

###############


# Coordenadas dos dois primeiros PCs para R30
pc1_pc2_R30 <- as.data.frame(pca_t_combined_matrix_R30$rotation[, 1:2])
colnames(pc1_pc2_R30) <- c("PC1", "PC2")

# Adicionar informação sobre o grupo (se houver), caso você tenha algum fator de agrupamento, por exemplo, batch ou tipo de amostra
# pc1_pc2_R30$group <- ...

# Plotando PC1 vs PC2
png("images/PCA_PC1_PC2_R30.png", width = 1920, height = 1080, res = 200)
ggplot(pc1_pc2_R30, aes(x = PC1, y = PC2)) +
  geom_point(color = "red") +  # Você pode adicionar coloração conforme algum critério, como o grupo
  labs(
    title = "Distribuição das Amostras - PCA - R30",
    x = "PC1",
    y = "PC2"
  ) +
  theme_minimal()
dev.off()


##############
# CORRIGINDO DE X PARA ROTATION
# Coordenadas de PC1 e PC2 para R0 e R30
pc1_pc2_combined <- rbind(
  cbind(pc1_pc2_R0, Radius = "R0"),
  cbind(pc1_pc2_R30, Radius = "R30")
)

# Plotando ambos os dados no mesmo gráfico
png("images/PCA_PC1_PC2_Combined_R0_R30.png", width = 1920, height = 1080, res = 200)
ggplot(pc1_pc2_combined, aes(x = PC1, y = PC2, color = Radius)) +
  geom_point() +
  labs(
    title = "Distribuição das Amostras - PCA - R0 e R30",
    x = "PC1",
    y = "PC2"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "red"))
dev.off()



# Coordenadas de dois PCs para R0 (PC2 e PC3)
pc2_pc3_R0 <- as.data.frame(pca_t_combined_matrix_R0$x[, c(2, 3)])  # PC2 e PC3
colnames(pc2_pc3_R0) <- c("PC2", "PC3")

# Plotando PC2 vs PC3 para R0
png("images/PCA_PC2_PC3_R0.png", width = 1920, height = 1080, res = 200)
ggplot(pc2_pc3_R0, aes(x = PC2, y = PC3)) +
  geom_point(color = "blue") +  # Cor azul para R0
  labs(
    title = "Distribuição das Amostras - PCA - PC2 vs PC3 - R0",
    x = "PC2",
    y = "PC3"
  ) +
  theme_minimal()
dev.off()

# Coordenadas de dois PCs para R30 (PC2 e PC3)
pc2_pc3_R30 <- as.data.frame(pca_t_combined_matrix_R30$x[, c(2, 3)])  # PC2 e PC3
colnames(pc2_pc3_R30) <- c("PC2", "PC3")

# Plotando PC2 vs PC3 para R30
png("images/PCA_PC2_PC3_R30.png", width = 1920, height = 1080, res = 200)
ggplot(pc2_pc3_R30, aes(x = PC2, y = PC3)) +
  geom_point(color = "red") +  # Cor vermelha para R30
  labs(
    title = "Distribuição das Amostras - PCA - PC2 vs PC3 - R30",
    x = "PC2",
    y = "PC3"
  ) +
  theme_minimal()
dev.off()

# Coordenadas de PC2 e PC3 para R0 e R30
pc2_pc3_combined <- rbind(
  cbind(pc2_pc3_R0, Radius = "R0"),
  cbind(pc2_pc3_R30, Radius = "R30")
)

# Plotando ambos os dados no mesmo gráfico (PC2 vs PC3)
png("images/PCA_PC2_PC3_Combined_R0_R30.png", width = 1920, height = 1080, res = 200)
ggplot(pc2_pc3_combined, aes(x = PC2, y = PC3, color = Radius)) +
  geom_point() +
  labs(
    title = "Distribuição das Amostras - PCA - PC2 vs PC3 - R0 e R30",
    x = "PC2",
    y = "PC3"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "red"))
dev.off()


# Coordenadas de dois PCs para R0 (PC1 e PC3)
pc1_pc3_R0 <- as.data.frame(pca_t_combined_matrix_R0$x[, c(1, 3)])  # PC1 e PC3
colnames(pc1_pc3_R0) <- c("PC1", "PC3")

# Plotando PC1 vs PC3 para R0
png("images/PCA_PC1_PC3_R0.png", width = 1920, height = 1080, res = 200)
ggplot(pc1_pc3_R0, aes(x = PC1, y = PC3)) +
  geom_point(color = "blue") +  # Cor azul para R0
  labs(
    title = "Distribuição das Amostras - PCA - PC1 vs PC3 - R0",
    x = "PC1",
    y = "PC3"
  ) +
  theme_minimal()
dev.off()

# Coordenadas de dois PCs para R30 (PC1 e PC3)
pc1_pc3_R30 <- as.data.frame(pca_t_combined_matrix_R30$x[, c(1, 3)])  # PC1 e PC3
colnames(pc1_pc3_R30) <- c("PC1", "PC3")

# Plotando PC1 vs PC3 para R30
png("images/PCA_PC1_PC3_R30.png", width = 1920, height = 1080, res = 200)
ggplot(pc1_pc3_R30, aes(x = PC1, y = PC3)) +
  geom_point(color = "red") +  # Cor vermelha para R30
  labs(
    title = "Distribuição das Amostras - PCA - PC1 vs PC3 - R30",
    x = "PC1",
    y = "PC3"
  ) +
  theme_minimal()
dev.off()

# Coordenadas de PC1 e PC3 para R0 e R30
pc1_pc3_combined <- rbind(
  cbind(pc1_pc3_R0, Radius = "R0"),
  cbind(pc1_pc3_R30, Radius = "R30")
)

# Plotando ambos os dados no mesmo gráfico (PC1 vs PC3)
png("images/PCA_PC1_PC3_Combined_R0_R30.png", width = 1920, height = 1080, res = 200)
ggplot(pc1_pc3_combined, aes(x = PC1, y = PC3, color = Radius)) +
  geom_point() +
  labs(
    title = "Distribuição das Amostras - PCA - PC1 vs PC3 - R0 e R30",
    x = "PC1",
    y = "PC3"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "red"))
dev.off()

#### Acessando o rotation

# Acessando a rotação para o objeto pca_t_combined_matrix_R0
rotation_R0 <- pca_t_combined_matrix_R0$rotation

# Acessando a rotação para o objeto pca_t_combined_matrix_R30
rotation_R30 <- pca_t_combined_matrix_R30$rotation


# Selecionar somente as colunas PC1, PC2 e PC3
rotation_R0_selected <- rotation_R0[, c("PC1", "PC2", "PC3")]
rotation_R30_selected <- rotation_R30[, c("PC1", "PC2", "PC3")]

# Exportar para arquivos CSV
write.csv(rotation_R0_selected, "rotation_R0_PC1_PC2_PC3.csv", row.names = TRUE)
write.csv(rotation_R30_selected, "rotation_R30_PC1_PC2_PC3.csv", row.names = TRUE)

#### Acessando o component

# Acessar os scores (componentes principais) para R0
components_R0 <- pca_t_combined_matrix_R0$x

# Visualizar as primeiras linhas e colunas
head(components_R0[, 1:3])  # Mostra PC1, PC2, e PC3

# Exportar os componentes principais para CSV
write.csv(components_R0[, 1:3], "components_R0_PC1_PC2_PC3.csv", row.names = TRUE)


# Acessar os scores (componentes principais) para R0
components_R30 <- pca_t_combined_matrix_R30$x

# Visualizar as primeiras linhas e colunas
head(components_R30[, 1:3])  # Mostra PC1, PC2, e PC3

# Exportar os componentes principais para CSV
write.csv(components_R30[, 1:3], "components_R30_PC1_PC2_PC3.csv", row.names = TRUE)

