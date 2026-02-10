# Instalar e carregar pacotes necessários

library(Seurat)
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

# Número total de núcleos disponíveis
num_cores <- detectCores()
cat("Núcleos disponíveis:", num_cores, "\n")

# Reservar 1 núcleo para o sistema (opcional)
num_cores_to_use <- max(1, num_cores - 1)
cat("Núcleos a serem usados:", num_cores_to_use, "\n")

# Criar cluster com núcleos disponíveis
cl <- makeCluster(num_cores_to_use)

# Exemplo de paralelização usando parLapply
results <- parLapply(cl, 1:10, function(x) x^2)
print(results)

# Encerrar o cluster
stopCluster(cl)



# Para salvar vários objetos com nomes similares no ambiente do R, você pode usar uma função que busca e salva cada objeto em um arquivo separado com base em um padrão comum no nome. Por exemplo, se você tiver várias variáveis no ambiente com nomes começando com "dataset_" e quiser salvar cada uma delas em um arquivo .RData, aqui está como fazer:
#   
# Exemplo: Salvando Objetos com Nomes Similares em Arquivos .RData
# r
# Copiar código
# Função para salvar objetos com um padrão no nome
salvar_objetos_rdata <- function(...) {
  # Captura os nomes dos objetos passados como argumentos
  objetos <- list(...)
  
  # Loop para salvar cada objeto no arquivo com seu próprio nome
  for (obj_nome in objetos) {
    save(list = obj_nome, file = paste0(obj_nome, ".RData"))
  }
  
  message("Objeto(s) ", objetos," salvo(s) com sucesso.")
}

# Exemplo de uso: salvar todos os objetos que começam com "dataset_"
# salvar_objetos("dataset_")


setwd("~/mestrado-single-cell/mestrado-single-cell")

# 1. Pré-processamento dos Dados de Single-Cell com Seurat

# Carregar dados de scRNA-seq (exemplo usando dados integrados no Seurat)

# data_teste <- Read10X(data.dir = "mtx_conversions/notreated_batch1/outs/filtered_feature_bc_matrix/")

##################### Pedido da Profa Rita em 18.11.2024 ########################

process_single_cell_data <- function(data_dir_name) {
  # Montar o caminho do diretório
  data_dir <- paste0("mtx_conversions/", data_dir_name, "/outs/filtered_feature_bc_matrix")
  
  # Ler os dados de expressão
  expression_data <- Read10X(data.dir = data_dir)
  
  # Criar o objeto Seurat
  seurat_object <- CreateSeuratObject(counts = expression_data)
  
  # Estatísticas básicas
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")  # Se houver genes mitocondriais
  
  # Normalização dos Dados
  seurat_object <- NormalizeData(seurat_object)
  
  # Identificação de Genes Variáveis
  seurat_object <- FindVariableFeatures(seurat_object)
  
  # Escalonar os Dados
  seurat_object <- ScaleData(seurat_object)
  
  # Redução de Dimensionalidade (PCA)
  seurat_object_with_PCA <- RunPCA(seurat_object)
  
  # Criar o Grafo de Vizinhança
  seurat_object_with_PCA <- FindNeighbors(seurat_object_with_PCA, dims = 1:10)  # Use as primeiras 10 PCs
  
  # Clustering de Células
  seurat_object_with_PCA <- FindClusters(seurat_object_with_PCA)
  
  # Rodar o UMAP
  seurat_object_with_PCA <- RunUMAP(seurat_object_with_PCA, dims = 1:10)
  
  # Visualizar os clusters com UMAP
  DimPlot(seurat_object_with_PCA, reduction = "umap", group.by = "seurat_clusters")
  
  # Identificar genes diferencialmente expressos entre todos os clusters
  markers <- FindAllMarkers(seurat_object_with_PCA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  # Exibir os principais genes marcadores
  print(head(markers, 30))
  
  # Exibir o número de células em cada cluster
  print(table(seurat_object_with_PCA$seurat_clusters))
  
  # Retornar o objeto Seurat final para possíveis usos futuros
  return(seurat_object_with_PCA)
}

# Chamar a função passando o nome do diretório (sem o caminho completo)
notreated_batch1_process_single_cell_data_with_PCA <- process_single_cell_data("notreated_batch1")

### Transcriptogramer

data_notreated_batch1 <- Read10X("mtx_conversions/notreated_batch1/outs/filtered_feature_bc_matrix")
MCF10A_notreated_batch1 <- CreateSeuratObject(counts = data_notreated_batch1, project = "MCF10A_notreated_batch1")
MCF10A_notreated_batch1[["percent.mt"]] <- PercentageFeatureSet(MCF10A_notreated_batch1, pattern = "^MT-")  # Para dados humanos (genes mitocondriais começam com "MT-")
MCF10A_notreated_batch1 <- NormalizeData(MCF10A_notreated_batch1)
FeatureScatter(MCF10A_notreated_batch1, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.5)
FeatureScatter(MCF10A_notreated_batch1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
MCF10A_notreated_batch1 <- subset(MCF10A_notreated_batch1, subset = nFeature_RNA > 500 & percent.mt < 20)

# Acessar a matriz de dados normalizados preservando os rownames e colnames
MCF10A_notreated_batch1_matrix_normalized <- as.matrix(GetAssayData(MCF10A_notreated_batch1, slot = "data"))

## Transcriptograma de ambos -  MCF10A_notreated_batch1 e MCF10A_notreated_batch1

# MCF10A_notreated_batch1 - R=0

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# Mapear os nomes dos genes para ENSEMBL
gene_mapping <- getBM(
  attributes = c("external_gene_name", "ensembl_gene_id", "ensembl_peptide_id"),
  filters = "external_gene_name",
  values = rownames(MCF10A_notreated_batch1_matrix_normalized),
  mart = ensembl
)
#
# # Mapear os nomes dos genes na matriz para `ensembl_gene_id`
mapped_genes <- gene_mapping$ensembl_gene_id[match(
  rownames(MCF10A_notreated_batch1_matrix_normalized),
  gene_mapping$external_gene_name
)]
#
# # Atualizar os rownames da matriz
rownames(MCF10A_notreated_batch1_matrix_normalized) <- mapped_genes
#
# Remover genes não mapeados (NA)
MCF10A_notreated_batch1_matrix_normalized <- MCF10A_notreated_batch1_matrix_normalized[!is.na(rownames(MCF10A_notreated_batch1_matrix_normalized)), ]

# Verificar os novos rownames
head(rownames(MCF10A_notreated_batch1_matrix_normalized))


# Visualizar o mapeamento
# head(gene_mapping)



dictionary <- getBM(attributes = c("ensembl_peptide_id", "ensembl_gene_id"),
                    mart = ensembl)
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

# 
# # Crie um dicionário único para cada ENSG -> ENSP usando as colunas corretas
# unique_dictionary <- dictionary[!duplicated(dictionary$ensembl_gene_id), c("ensembl_peptide_id", "ensembl_gene_id")]
# 
# # Verifique a interseção entre os `rownames` de `combined_matrix` e `dictionary`
# identifiers_expression <- rownames(MCF10A_notreated_batch1_matrix_normalized)
# identifiers_dictionary <- unique_dictionary$ensembl_gene_id
# 
# # Identificadores que estão em comum
# common_identifiers <- intersect(identifiers_expression, identifiers_dictionary)
# cat("Número de identificadores em comum:", length(common_identifiers), "\n")
# 
# # Filtre `combined_matrix` para manter apenas os identificadores em comum
# combined_matrix_gene_common <- combined_matrix[rownames(combined_matrix) %in% common_identifiers, ]
# 
# # Filtre `unique_dictionary` para manter apenas os identificadores em comum
# unique_dictionary <- unique_dictionary[unique_dictionary$ensembl_gene_id %in% common_identifiers, ]
# 
# # Crie um mapeamento ENSG -> ENSP e substitua os `rownames` da matriz
# rownames_mapping <- setNames(unique_dictionary$ensembl_peptide_id, unique_dictionary$ensembl_gene_id)
# 
# # Substituir os rownames da matriz `combined_matrix` para ENSP
# rownames(combined_matrix_gene_common) <- rownames_mapping[rownames(combined_matrix_gene_common)]
# 

t_MCF10A_notreated_batch1_matrix_normalized <- transcriptogramPreprocess(association = assoc, 
                                                                         ordering = ord$Protein, radius = 0)

t_MCF10A_notreated_batch1_matrix_normalized <- transcriptogramStep1(object = t_MCF10A_notreated_batch1_matrix_normalized, 
                                                                    expression = MCF10A_notreated_batch1_matrix_normalized, 
                                                                    dictionary = dictionary)
t_MCF10A_notreated_batch1_matrix_normalized <- transcriptogramStep2(object = t_MCF10A_notreated_batch1_matrix_normalized)

salvar_objetos_rdata("t_MCF10A_notreated_batch1_matrix_normalized")

# MCF10A_notreated_batch2 - R=0

data_notreated_batch2 <- Read10X("mtx_conversions/notreated_batch2/outs/filtered_feature_bc_matrix")
MCF10A_notreated_batch2 <- CreateSeuratObject(counts = data_notreated_batch2, project = "MCF10A_notreated_batch2")
MCF10A_notreated_batch2[["percent.mt"]] <- PercentageFeatureSet(MCF10A_notreated_batch2, pattern = "^MT-")  # Para dados humanos (genes mitocondriais começam com "MT-")
MCF10A_notreated_batch2 <- NormalizeData(MCF10A_notreated_batch2)
FeatureScatter(MCF10A_notreated_batch2, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.5)
FeatureScatter(MCF10A_notreated_batch2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
MCF10A_notreated_batch2 <- subset(MCF10A_notreated_batch2, subset = nFeature_RNA > 500 & percent.mt < 20)

# Acessar a matriz de dados normalizados preservando os rownames e colnames
MCF10A_notreated_batch2_matrix_normalized <- as.matrix(GetAssayData(MCF10A_notreated_batch2, slot = "data"))

# Transcriptograma de ambos
# 
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# Mapear os nomes dos genes para ENSEMBL
gene_mapping <- getBM(
  attributes = c("external_gene_name", "ensembl_gene_id", "ensembl_peptide_id"),
  filters = "external_gene_name",
  values = rownames(MCF10A_notreated_batch2_matrix_normalized),
  mart = ensembl
)
#
# # Mapear os nomes dos genes na matriz para `ensembl_gene_id`
mapped_genes <- gene_mapping$ensembl_gene_id[match(
  rownames(MCF10A_notreated_batch2_matrix_normalized),
  gene_mapping$external_gene_name
)]
#
# # Atualizar os rownames da matriz
rownames(MCF10A_notreated_batch2_matrix_normalized) <- mapped_genes

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

t_MCF10A_notreated_batch2_matrix_normalized <- transcriptogramPreprocess(association = assoc, 
                                                                         ordering = ord$Protein, radius = 0)

t_MCF10A_notreated_batch2_matrix_normalized <- transcriptogramStep1(object = t_MCF10A_notreated_batch2_matrix_normalized, 
                                                                    expression = MCF10A_notreated_batch2_matrix_normalized, 
                                                                    dictionary = dictionary)
t_MCF10A_notreated_batch2_matrix_normalized <- transcriptogramStep2(object = t_MCF10A_notreated_batch2_matrix_normalized)

salvar_objetos_rdata("t_MCF10A_notreated_batch2_matrix_normalized ")

###############################################################################

### PCAs da Dia 0 Batchs 1 e 2

# Extraindo o dataframe do slot transcriptogramS2
df_notreated_batch1 <- t_MCF10A_notreated_batch1_matrix_normalized@transcriptogramS2

# Verificando o tipo de cada coluna
sapply(df_notreated_batch1, class)  # Para identificar as colunas que não são numéricas

# Mantendo apenas colunas numéricas
df_notreated_batch1 <- df_notreated_batch1[, sapply(df_notreated_batch1, is.numeric)]

# Verificando se sobrou algum dado
if (ncol(df_notreated_batch1) == 0) {
  stop("Nenhuma coluna numérica encontrada no dataframe!")
}

# Realizando a padronização dos dados numéricos
df_scaled_notreated_batch1 <- scale(df_notreated_batch1)

# Continuando com o PCA
pca_notreated_batch1 <- prcomp(df_scaled, center = TRUE, scale. = TRUE)

# Resumo dos resultados
summary(pca_notreated_batch1)

# Scree plot da variância explicada
plot(pca_notreated_batch1, type = "l", main = "Scree Plot - PCA - Notreated Batch 1")

# Extraindo o dataframe do slot transcriptogramS2
df_notreated_batch2 <- t_MCF10A_notreated_batch2_matrix_normalized@transcriptogramS2

# Verificando o tipo de cada coluna
sapply(df_notreated_batch2, class)  # Para identificar as colunas que não são numéricas

# Mantendo apenas colunas numéricas
df_notreated_batch2 <- df_notreated_batch2[, sapply(df_notreated_batch2, is.numeric)]

# Verificando se sobrou algum dado
if (ncol(df_notreated_batch2) == 0) {
  stop("Nenhuma coluna numérica encontrada no dataframe!")
}

# Realizando a padronização dos dados numéricos
df_scaled_notreated_batch2 <- scale(df_notreated_batch2)

# Continuando com o PCA
pca_notreated_batch2 <- prcomp(df_scaled_notreated_batch2, center = TRUE, scale. = TRUE)

# Resumo dos resultados
summary(pca_notreated_batch2)

# Scree plot da variância explicada
plot(pca_notreated_batch2, type = "l", main = "Scree Plot - PCA - Notreated Batch 2")


library(ggplot2)

# Variância explicada para o PCA Batch 1
explained_variance_batch1 <- pca_notreated_batch1$sdev^2 / sum(pca_notreated_batch1$sdev^2)

# Variância explicada para o PCA Batch 2
explained_variance_batch2 <- pca_notreated_batch2$sdev^2 / sum(pca_notreated_batch2$sdev^2)

# Verificar o número de PCs em cada PCA
length(explained_variance_batch1)  # Quantidade de PCs no PCA do Batch 1
length(explained_variance_batch2)  # Quantidade de PCs no PCA do Batch 2

# Ajustar o comprimento de ambos os vetores de variância explicada (preenchendo com NAs ou truncando)
# Vamos supor que você queira alinhar os PCs até o número máximo de componentes entre os dois batches.

max_pcs <- max(length(explained_variance_batch1), length(explained_variance_batch2))

# Preencher com NAs para garantir que ambos tenham o mesmo número de componentes
explained_variance_batch1 <- c(explained_variance_batch1, rep(NA, max_pcs - length(explained_variance_batch1)))
explained_variance_batch2 <- c(explained_variance_batch2, rep(NA, max_pcs - length(explained_variance_batch2)))

# Criar o dataframe novamente com os dados ajustados
scree_data <- data.frame(
  PC = rep(1:max_pcs, 2),
  Variance = c(explained_variance_batch1, explained_variance_batch2),
  Batch = rep(c("Batch 1", "Batch 2"), each = max_pcs)
)

# Variância cumulativa
scree_data$cumulative_variance <- ave(scree_data$Variance, scree_data$Batch, FUN = cumsum)

# Plot do Scree Plot com variância cumulativa
png("images/PCAeVarianciaAcumuladaTranscriptograma_NotreatedBatchs1e2.png", width = 1920, height = 1080, res = 200)
ggplot(scree_data, aes(x = PC, y = Variance, color = Batch)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Scree Plot com Variância Explicada - PCA",
    x = "Principal Component",
    y = "Explained Variance"
  ) +
  theme_minimal() +
  theme(legend.position = "top") +
  geom_line(aes(x = PC, y = cumulative_variance, color = Batch), linetype = "dashed") +
  scale_color_manual(values = c("blue", "red")) +
  scale_x_continuous(breaks = 1:max(scree_data$PC))
dev.off()

# Filtrando para mostrar apenas as TOP 10 PCs
scree_data_top10 <- scree_data[scree_data$PC <= 10, ]

# Plot do Scree Plot com variância cumulativa para as TOP 10 PCs
# Plot do Scree Plot com variância cumulativa

png("images/PCAeVarianciaAcumuladaTop10Transcriptograma_NotreatedBatchs1e2.png", width = 1920, height = 1080, res = 200)
ggplot(scree_data_top10, aes(x = PC, y = Variance, color = Batch)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Scree Plot com Variância Explicada - TOP 10 PCs",
    x = "Principal Component",
    y = "Explained Variance"
  ) +
  theme_minimal() +
  theme(legend.position = "top") +
  geom_line(aes(x = PC, y = cumulative_variance, color = Batch), linetype = "dashed") +
  scale_color_manual(values = c("blue", "red")) +
  scale_x_continuous(breaks = 1:10)  # Exibe apenas os 10 primeiros PCs no eixo X
dev.off()


#######################

# MCF10A_notreated_batch1 - R=30

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# Mapear os nomes dos genes para ENSEMBL
gene_mapping <- getBM(
  attributes = c("external_gene_name", "ensembl_gene_id", "ensembl_peptide_id"),
  filters = "external_gene_name",
  values = rownames(MCF10A_notreated_batch1_matrix_normalized),
  mart = ensembl
)
#
# # Mapear os nomes dos genes na matriz para `ensembl_gene_id`
mapped_genes <- gene_mapping$ensembl_gene_id[match(
  rownames(MCF10A_notreated_batch1_matrix_normalized),
  gene_mapping$external_gene_name
)]
#
# # Atualizar os rownames da matriz
rownames(MCF10A_notreated_batch1_matrix_normalized) <- mapped_genes

# Atualizar os rownames da matriz
rownames(MCF10A_notreated_batch1_matrix_normalized) <- mapped_genes
#
# Remover genes não mapeados (NA)
MCF10A_notreated_batch1_matrix_normalized <- MCF10A_notreated_batch1_matrix_normalized[!is.na(rownames(MCF10A_notreated_batch1_matrix_normalized)), ]

# Verificar os novos rownames
head(rownames(MCF10A_notreated_batch1_matrix_normalized))


# Visualizar o mapeamento
# head(gene_mapping)


dictionary <- getBM(attributes = c("ensembl_peptide_id", "ensembl_gene_id"),
                    mart = ensembl)
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


t_MCF10A_notreated_batch1_matrix_normalized_R30 <- transcriptogramPreprocess(association = assoc, 
                                                                         ordering = ord$Protein, radius = 30)

t_MCF10A_notreated_batch1_matrix_normalized_R30 <- transcriptogramStep1(object = t_MCF10A_notreated_batch1_matrix_normalized_R30, 
                                                                    expression = MCF10A_notreated_batch1_matrix_normalized, 
                                                                    dictionary = dictionary)
t_MCF10A_notreated_batch1_matrix_normalized_R30 <- transcriptogramStep2(object = t_MCF10A_notreated_batch1_matrix_normalized_R30)


salvar_objetos_rdata("t_MCF10A_notreated_batch1_matrix_normalized_R30")

# MCF10A_notreated_batch2 - R=30

data_notreated_batch2 <- Read10X("mtx_conversions/notreated_batch2/outs/filtered_feature_bc_matrix")
MCF10A_notreated_batch2 <- CreateSeuratObject(counts = data_notreated_batch2, project = "MCF10A_notreated_batch2")
MCF10A_notreated_batch2[["percent.mt"]] <- PercentageFeatureSet(MCF10A_notreated_batch2, pattern = "^MT-")  # Para dados humanos (genes mitocondriais começam com "MT-")
MCF10A_notreated_batch2 <- NormalizeData(MCF10A_notreated_batch2)
FeatureScatter(MCF10A_notreated_batch2, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.5)
FeatureScatter(MCF10A_notreated_batch2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
MCF10A_notreated_batch2 <- subset(MCF10A_notreated_batch2, subset = nFeature_RNA > 500 & percent.mt < 20)

# Acessar a matriz de dados normalizados preservando os rownames e colnames
MCF10A_notreated_batch2_matrix_normalized <- as.matrix(GetAssayData(MCF10A_notreated_batch2, slot = "data"))

# Transcriptograma de ambos
# 
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# Mapear os nomes dos genes para ENSEMBL
gene_mapping <- getBM(
  attributes = c("external_gene_name", "ensembl_gene_id", "ensembl_peptide_id"),
  filters = "external_gene_name",
  values = rownames(MCF10A_notreated_batch2_matrix_normalized),
  mart = ensembl
)
#
# # Mapear os nomes dos genes na matriz para `ensembl_gene_id`
mapped_genes <- gene_mapping$ensembl_gene_id[match(
  rownames(MCF10A_notreated_batch2_matrix_normalized),
  gene_mapping$external_gene_name
)]
#
# # Atualizar os rownames da matriz
rownames(MCF10A_notreated_batch2_matrix_normalized) <- mapped_genes

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

t_MCF10A_notreated_batch2_matrix_normalized_R30 <- transcriptogramPreprocess(association = assoc, 
                                                                         ordering = ord$Protein, radius = 30)

t_MCF10A_notreated_batch2_matrix_normalized_R30 <- transcriptogramStep1(object = t_MCF10A_notreated_batch2_matrix_normalized_R30, 
                                                                    expression = MCF10A_notreated_batch2_matrix_normalized, 
                                                                    dictionary = dictionary)
t_MCF10A_notreated_batch2_matrix_normalized_R30 <- transcriptogramStep2(object = t_MCF10A_notreated_batch2_matrix_normalized_R30)

salvar_objetos_rdata("t_MCF10A_notreated_batch2_matrix_normalized_R30")



###############################################################################

### PCAs da Dia 0 Batchs 1 e 2 --- R=30

# Extraindo o dataframe do slot transcriptogramS2
df_notreated_batch1 <- t_MCF10A_notreated_batch1_matrix_normalized_R30@transcriptogramS2

# Verificando o tipo de cada coluna
sapply(df_notreated_batch1, class)  # Para identificar as colunas que não são numéricas

# Mantendo apenas colunas numéricas
df_notreated_batch1 <- df_notreated_batch1[, sapply(df_notreated_batch1, is.numeric)]

# Verificando se sobrou algum dado
if (ncol(df_notreated_batch1) == 0) {
  stop("Nenhuma coluna numérica encontrada no dataframe!")
}

# Realizando a padronização dos dados numéricos
df_scaled_notreated_batch1 <- scale(df_notreated_batch1)

# Continuando com o PCA
pca_notreated_batch1 <- prcomp(df_scaled_notreated_batch1, center = TRUE, scale. = TRUE)

# Resumo dos resultados
summary(pca_notreated_batch1)

# Scree plot da variância explicada
plot(pca_notreated_batch1, type = "l", main = "Scree Plot - PCA - Notreated Batch 1 - R = 30")

# Extraindo o dataframe do slot transcriptogramS2
df_notreated_batch2 <- t_MCF10A_notreated_batch2_matrix_normalized_R30@transcriptogramS2

# Verificando o tipo de cada coluna
sapply(df_notreated_batch2, class)  # Para identificar as colunas que não são numéricas

# Mantendo apenas colunas numéricas
df_notreated_batch2 <- df_notreated_batch2[, sapply(df_notreated_batch2, is.numeric)]

# Verificando se sobrou algum dado
if (ncol(df_notreated_batch2) == 0) {
  stop("Nenhuma coluna numérica encontrada no dataframe!")
}

# Realizando a padronização dos dados numéricos
df_scaled_notreated_batch2 <- scale(df_notreated_batch2)

# Continuando com o PCA
pca_notreated_batch2 <- prcomp(df_scaled_notreated_batch2, center = TRUE, scale. = TRUE)

# Resumo dos resultados
summary(pca_notreated_batch2)

# Scree plot da variância explicada
plot(pca_notreated_batch2, type = "l", main = "Scree Plot - PCA - Notreated Batch 2 - R = 30")


library(ggplot2)

# Variância explicada para o PCA Batch 1
explained_variance_batch1 <- pca_notreated_batch1$sdev^2 / sum(pca_notreated_batch1$sdev^2)

# Variância explicada para o PCA Batch 2
explained_variance_batch2 <- pca_notreated_batch2$sdev^2 / sum(pca_notreated_batch2$sdev^2)

# Verificar o número de PCs em cada PCA
length(explained_variance_batch1)  # Quantidade de PCs no PCA do Batch 1
length(explained_variance_batch2)  # Quantidade de PCs no PCA do Batch 2

# Ajustar o comprimento de ambos os vetores de variância explicada (preenchendo com NAs ou truncando)
# Vamos supor que você queira alinhar os PCs até o número máximo de componentes entre os dois batches.

max_pcs <- max(length(explained_variance_batch1), length(explained_variance_batch2))

# Preencher com NAs para garantir que ambos tenham o mesmo número de componentes
explained_variance_batch1 <- c(explained_variance_batch1, rep(NA, max_pcs - length(explained_variance_batch1)))
explained_variance_batch2 <- c(explained_variance_batch2, rep(NA, max_pcs - length(explained_variance_batch2)))

# Criar o dataframe novamente com os dados ajustados
scree_data <- data.frame(
  PC = rep(1:max_pcs, 2),
  Variance = c(explained_variance_batch1, explained_variance_batch2),
  Batch = rep(c("Batch 1", "Batch 2"), each = max_pcs)
)

# Variância cumulativa
scree_data$cumulative_variance <- ave(scree_data$Variance, scree_data$Batch, FUN = cumsum)

# Plot do Scree Plot com variância cumulativa
png("images/PCAeVarianciaAcumuladaTranscriptograma_NotreatedBatchs1e2_R30.png", width = 1920, height = 1080, res = 200)
ggplot(scree_data, aes(x = PC, y = Variance, color = Batch)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Scree Plot com Variância Explicada - PCA - R = 30",
    x = "Principal Component",
    y = "Explained Variance"
  ) +
  theme_minimal() +
  theme(legend.position = "top") +
  geom_line(aes(x = PC, y = cumulative_variance, color = Batch), linetype = "dashed") +
  scale_color_manual(values = c("blue", "red")) +
  scale_x_continuous(breaks = 1:max(scree_data$PC))
dev.off()

# Filtrando para mostrar apenas as TOP 10 PCs
scree_data_top10 <- scree_data[scree_data$PC <= 10, ]

# Plot do Scree Plot com variância cumulativa para as TOP 10 PCs
# Plot do Scree Plot com variância cumulativa

png("images/PCAeVarianciaAcumuladaTop10Transcriptograma_NotreatedBatchs1e2_R30.png", width = 1920, height = 1080, res = 200)
ggplot(scree_data_top10, aes(x = PC, y = Variance, color = Batch)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Scree Plot com Variância Explicada - TOP 10 PCs - R = 30",
    x = "Principal Component",
    y = "Explained Variance"
  ) +
  theme_minimal() +
  theme(legend.position = "top") +
  geom_line(aes(x = PC, y = cumulative_variance, color = Batch), linetype = "dashed") +
  scale_color_manual(values = c("blue", "red")) +
  scale_x_continuous(breaks = 1:10)  # Exibe apenas os 10 primeiros PCs no eixo X
dev.off()

############## REUNIÃO 27/11/2024 #### PCAS sem a coluna 2


# Extraindo o dataframe do slot transcriptogramS2
df_notreated_batch1 <- t_MCF10A_notreated_batch1_matrix_normalized_R30@transcriptogramS2[,-2]

# Verificando o tipo de cada coluna
sapply(df_notreated_batch1, class)  # Para identificar as colunas que não são numéricas

# Mantendo apenas colunas numéricas
df_notreated_batch1 <- df_notreated_batch1[, sapply(df_notreated_batch1, is.numeric)]

# Verificando se sobrou algum dado
if (ncol(df_notreated_batch1) == 0) {
  stop("Nenhuma coluna numérica encontrada no dataframe!")
}

# Realizando a padronização dos dados numéricos
df_scaled_notreated_batch1 <- scale(df_notreated_batch1)

# Continuando com o PCA
pca_notreated_batch1 <- prcomp(df_scaled_notreated_batch1, center = TRUE, scale. = TRUE)


summary(pca_notreated_batch1)

# Scree plot da variância explicada
plot(pca_notreated_batch1, type = "l", main = "Scree Plot - PCA - Notreated Batch 1 - R = 30")

t_MCF10A_notreated_batch1_matrix_normalized_R30@transcriptogramS2# Extraindo o dataframe do slot transcriptogramS2
df_notreated_batch2 <- t_MCF10A_notreated_batch2_matrix_normalized_R30@transcriptogramS2[,-2]

# Verificando o tipo de cada coluna
sapply(df_notreated_batch2, class)  # Para identificar as colunas que não são numéricas

# Mantendo apenas colunas numéricas
df_notreated_batch2 <- df_notreated_batch2[, sapply(df_notreated_batch2, is.numeric)]

# Verificando se sobrou algum dado
if (ncol(df_notreated_batch2) == 0) {
  stop("Nenhuma coluna numérica encontrada no dataframe!")
}

# Realizando a padronização dos dados numéricos
df_scaled_notreated_batch2 <- scale(df_notreated_batch2)

# Continuando com o PCA
pca_notreated_batch2 <- prcomp(df_scaled_notreated_batch2, center = TRUE, scale. = TRUE)

# Resumo dos resultados
summary(pca_notreated_batch2)

# Scree plot da variância explicada
plot(pca_notreated_batch2, type = "l", main = "Scree Plot - PCA - Notreated Batch 2 - R = 30")


library(ggplot2)

# Variância explicada para o PCA Batch 1
explained_variance_batch1 <- pca_notreated_batch1$sdev^2 / sum(pca_notreated_batch1$sdev^2)

# Variância explicada para o PCA Batch 2
explained_variance_batch2 <- pca_notreated_batch2$sdev^2 / sum(pca_notreated_batch2$sdev^2)

# Verificar o número de PCs em cada PCA
length(explained_variance_batch1)  # Quantidade de PCs no PCA do Batch 1
length(explained_variance_batch2)  # Quantidade de PCs no PCA do Batch 2

# Ajustar o comprimento de ambos os vetores de variância explicada (preenchendo com NAs ou truncando)
# Vamos supor que você queira alinhar os PCs até o número máximo de componentes entre os dois batches.

max_pcs <- max(length(explained_variance_batch1), length(explained_variance_batch2))

# Preencher com NAs para garantir que ambos tenham o mesmo número de componentes
explained_variance_batch1 <- c(explained_variance_batch1, rep(NA, max_pcs - length(explained_variance_batch1)))
explained_variance_batch2 <- c(explained_variance_batch2, rep(NA, max_pcs - length(explained_variance_batch2)))

# Criar o dataframe novamente com os dados ajustados
scree_data <- data.frame(
  PC = rep(1:max_pcs, 2),
  Variance = c(explained_variance_batch1, explained_variance_batch2),
  Batch = rep(c("Batch 1", "Batch 2"), each = max_pcs)
)

# Variância cumulativa
scree_data$cumulative_variance <- ave(scree_data$Variance, scree_data$Batch, FUN = cumsum)

# Plot do Scree Plot com variância cumulativa
png("images/PCAeVarianciaAcumuladaTranscriptograma_NotreatedBatchs1e2_R30_WithoutCol2.png", width = 1920, height = 1080, res = 200)
ggplot(scree_data, aes(x = PC, y = Variance, color = Batch)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Scree Plot com Variância Explicada - PCA - R = 30",
    x = "Principal Component",
    y = "Explained Variance"
  ) +
  theme_minimal() +
  theme(legend.position = "top") +
  geom_line(aes(x = PC, y = cumulative_variance, color = Batch), linetype = "dashed") +
  scale_color_manual(values = c("blue", "red")) +
  scale_x_continuous(breaks = 1:max(scree_data$PC))
dev.off()

# Filtrando para mostrar apenas as TOP 10 PCs
scree_data_top10 <- scree_data[scree_data$PC <= 10, ]

# Plot do Scree Plot com variância cumulativa para as TOP 10 PCs
# Plot do Scree Plot com variância cumulativa

png("images/PCAeVarianciaAcumuladaTop10Transcriptograma_NotreatedBatchs1e2_R30_WithoutCol2.png", width = 1920, height = 1080, res = 200)
ggplot(scree_data_top10, aes(x = PC, y = Variance, color = Batch)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Scree Plot com Variância Explicada - TOP 10 PCs - R = 30",
    x = "Principal Component",
    y = "Explained Variance"
  ) +
  theme_minimal() +
  theme(legend.position = "top") +
  geom_line(aes(x = PC, y = cumulative_variance, color = Batch), linetype = "dashed") +
  scale_color_manual(values = c("blue", "red")) +
  scale_x_continuous(breaks = 1:10)  # Exibe apenas os 10 primeiros PCs no eixo X
dev.off()

################################################################################3

### Aqui, usando Seurat e Transcriptogramer. Depois, faremos a PCA e o plot.


# MCF10A_TGFbeta1_1day_batch2

data_TGFbeta1_1day_batch2 <- Read10X("mtx_conversions/TGFbeta1_1day_batch2/outs/filtered_feature_bc_matrix")
MCF10A_TGFbeta1_1day_batch2 <- CreateSeuratObject(counts = data_TGFbeta1_1day_batch2, project = "MCF10A_TGFbeta1_1day_batch2")
MCF10A_TGFbeta1_1day_batch2[["percent.mt"]] <- PercentageFeatureSet(MCF10A_TGFbeta1_1day_batch2, pattern = "^MT-")  # Para dados humanos (genes mitocondriais começam com "MT-")
MCF10A_TGFbeta1_1day_batch2 <- NormalizeData(MCF10A_TGFbeta1_1day_batch2)
FeatureScatter(MCF10A_TGFbeta1_1day_batch2, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.5)
FeatureScatter(MCF10A_TGFbeta1_1day_batch2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
MCF10A_TGFbeta1_1day_batch2 <- subset(MCF10A_TGFbeta1_1day_batch2, subset = nFeature_RNA > 500 & percent.mt < 20)

# Acessar a matriz de dados normalizados preservando os rownames e colnames
MCF10A_TGFbeta1_1day_batch2_matrix_normalized <- as.matrix(GetAssayData(MCF10A_TGFbeta1_1day_batch2, slot = "data"))

# Transcriptograma de ambos
# 
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# Mapear os nomes dos genes para ENSEMBL
gene_mapping <- getBM(
  attributes = c("external_gene_name", "ensembl_gene_id", "ensembl_peptide_id"),
  filters = "external_gene_name",
  values = rownames(MCF10A_TGFbeta1_1day_batch2_matrix_normalized),
  mart = ensembl
)
#
# # Mapear os nomes dos genes na matriz para `ensembl_gene_id`
mapped_genes <- gene_mapping$ensembl_gene_id[match(
  rownames(MCF10A_TGFbeta1_1day_batch2_matrix_normalized),
  gene_mapping$external_gene_name
)]
#
# # Atualizar os rownames da matriz
rownames(MCF10A_TGFbeta1_1day_batch2_matrix_normalized) <- mapped_genes

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

t_MCF10A_TGFbeta1_1day_batch2_matrix_normalized <- transcriptogramPreprocess(association = assoc, 
                                                                         ordering = ord$Protein, radius = 30)

t_MCF10A_TGFbeta1_1day_batch2_matrix_normalized <- transcriptogramStep1(object = t_MCF10A_TGFbeta1_1day_batch2_matrix_normalized, 
                                                                    expression = MCF10A_TGFbeta1_1day_batch2_matrix_normalized, 
                                                                    dictionary = dictionary)
t_MCF10A_TGFbeta1_1day_batch2_matrix_normalized <- transcriptogramStep2(object = t_MCF10A_TGFbeta1_1day_batch2_matrix_normalized)

salvar_objetos_rdata("t_MCF10A_TGFbeta1_1day_batch2_matrix_normalized")





#########################################
#########################################
expression_data_notreated_batch2 <- Read10X("mtx_conversions/notreated_batch2/outs/filtered_feature_bc_matrix")
MCF10A_notreated_batch2 <- CreateSeuratObject(counts = expression_data_notreated_batch2, project = "MCF10A_notreated_batch2")
MCF10A_notreated_batch2[["percent.mt"]] <- PercentageFeatureSet(MCF10A_notreated_batch2, pattern = "^MT-")
FeatureScatter(MCF10A_notreated_batch2, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.5)
FeatureScatter(MCF10A_notreated_batch2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
MCF10A_notreated_batch2 <- subset(MCF10A_notreated_batch2, subset = nFeature_RNA > 500 & percent.mt < 20)



# Converter para uma matriz densa (regular)
expression_matrix_dense_notreated_batch1 <- as.matrix(expression_data_notreated_batch1)
expression_matrix_dense_notreated_batch2 <- as.matrix(expression_data_notreated_batch2)



expression_data <- expression_matrix_dense_notreated_batch1 

seurat_object <- CreateSeuratObject(counts = expression_data)

# Rodar Transcriptogramer com o dicionário e a associação
seurat_object_with_transcriptogram <- transcriptogramStep1(seurat_object, 
                                                           dictionary = dictionary,
                                                           association = assoc, 
                                                           ncores = 255)

# 7. Verificar os resultados
# Verifique os resultados dos cálculos do Transcriptogramer
seurat_object_with_transcriptogram

##################################################################################


expression_data_notreated_batch1 <- Read10X("mtx_conversions/notreated_batch1/outs/filtered_feature_bc_matrix")
expression_data_notreated_batch2 <- Read10X("mtx_conversions/notreated_batch2/outs/filtered_feature_bc_matrix")
expression_data_TGFbeta1_1day_batch2 <- Read10X("mtx_conversions/TGFbeta1_1day_batch2/outs/filtered_feature_bc_matrix")
expression_data_TGFbeta1_2day_batch2 <- Read10X("mtx_conversions/TGFbeta1_2day_batch2/outs/filtered_feature_bc_matrix")
expression_data_TGFbeta1_3day_batch2 <- Read10X("mtx_conversions/TGFbeta1_3day_batch2/outs/filtered_feature_bc_matrix")
expression_data_TGFbeta1_4day_batch1 <- Read10X("mtx_conversions/TGFbeta1_4day_batch1/outs/filtered_feature_bc_matrix")
expression_data_TGFbeta1_8day_batch1 <- Read10X("mtx_conversions/TGFbeta1_8day_batch1/outs/filtered_feature_bc_matrix")


dense_matrix_notreated_batch1 <- as.matrix(expression_data_notreated_batch1)


# salvar_objetos_rdata("expression_data_notreated_batch1",
#                      "expression_data_notreated_batch2",
#                      "expression_data_TGFbeta1_1day_batch2",
#                      "expression_data_TGFbeta1_2day_batch2",
#                      "expression_data_TGFbeta1_3day_batch2",
#                      "expression_data_TGFbeta1_4day_batch1",
#                      "expression_data_TGFbeta1_8day_batch1")


# Para normalizar manualmente uma matriz esparsa de classe dgCMatrix no R (como 
# as que são frequentemente usadas em análises de scRNA-seq) com base no número 
# de counts de cada célula, você pode seguir os seguintes passos:
#   
# Obter a soma total de counts por célula: A soma de todas as expressões gênicas
# de cada célula será utilizada como fator de normalização.
# Dividir cada gene pelo total de counts da célula: Cada valor de expressão será
# dividido pela soma dos counts daquela célula.
# 
# # A soma dos counts por célula (colunas) é obtida usando a função colSums().
# cell_sums_notreated_batch1 <- Matrix::colSums(matrix_notreated_batch1)
# 
# 
# # Normalizar os dados dividindo os valores de cada célula pela soma dos counts:
# # Evita divisão por zero em células que podem ter 0 counts
# cell_sums_notreated_batch1[cell_sums_notreated_batch1 == 0] <- 1
# 
# # Normalizar dividindo cada coluna pela soma correspondente
# mat_normalized_notreated_batch1 <- t(t(matrix_notreated_batch1) / cell_sums_notreated_batch1)
# 
# # Verificando a soma das colunas após normalização (deverá ser próximo de 1 para todas as células)
# colSums(mat_normalized_notreated_batch1)


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

library(biomaRt)

# Conectando ao Ensembl
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Obter genes mitocondriais
mito_genes_info <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
  filters = "chromosome_name",
  values = "MT",  # Mitocôndria
  mart = mart
)

# Vetor de genes mitocondriais
mito_genes <- mito_genes_info$ensembl_gene_id

mito_genes_present_matrix_notreated_batch1 <- rownames(matrix_notreated_batch1) %in% mito_genes
mito_genes_present_matrix_notreated_batch2 <- rownames(matrix_notreated_batch2) %in% mito_genes
mito_genes_present_matrix_TGFbeta1_1day_batch2 <- rownames(matrix_TGFbeta1_1day_batch2) %in% mito_genes
mito_genes_present_matrix_TGFbeta1_2day_batch2 <- rownames(matrix_TGFbeta1_2day_batch2) %in% mito_genes
mito_genes_present_matrix_TGFbeta1_3day_batch2 <- rownames(matrix_TGFbeta1_3day_batch2) %in% mito_genes
mito_genes_present_matrix_TGFbeta1_4day_batch1 <- rownames(matrix_TGFbeta1_4day_batch1) %in% mito_genes
mito_genes_present_matrix_TGFbeta1_8day_batch1 <- rownames(matrix_TGFbeta1_8day_batch1) %in% mito_genes

# Subconjunto apenas com genes mitocondriais
mito_expression_matrix_notreated_batch1 <- matrix_notreated_batch1[mito_genes_present_matrix_notreated_batch1, ]
mito_expression_matrix_notreated_batch2 <- matrix_notreated_batch2[mito_genes_present_matrix_notreated_batch2, ]
mito_genes_present_matrix_TGFbeta1_1day_batch2 <- matrix_TGFbeta1_1day_batch2[mito_genes_present_matrix_TGFbeta1_1day_batch2, ]
mito_genes_present_matrix_TGFbeta1_2day_batch2 <- matrix_TGFbeta1_2day_batch2[mito_genes_present_matrix_TGFbeta1_2day_batch2, ]
mito_genes_present_matrix_TGFbeta1_3day_batch2 <- matrix_TGFbeta1_3day_batch2[mito_genes_present_matrix_TGFbeta1_3day_batch2, ]
mito_genes_present_matrix_TGFbeta1_4day_batch1 <- matrix_TGFbeta1_4day_batch1[mito_genes_present_matrix_TGFbeta1_4day_batch1, ]
mito_genes_present_matrix_TGFbeta1_8day_batch1 <- matrix_TGFbeta1_8day_batch1[mito_genes_present_matrix_TGFbeta1_8day_batch1, ]

# Total de expressão por célula
total_expression_matrix_notreated_batch1 <- colSums(matrix_notreated_batch1)
total_expression_matrix_notreated_batch2 <- colSums(matrix_notreated_batch2)
total_expression_matrix_TGFbeta1_1day_batch2 <- colSums(matrix_TGFbeta1_1day_batch2)
total_expression_matrix_TGFbeta1_2day_batch2 <- colSums(matrix_TGFbeta1_2day_batch2)
total_expression_matrix_TGFbeta1_3day_batch2 <- colSums(matrix_TGFbeta1_3day_batch2)
total_expression_matrix_TGFbeta1_4day_batch1 <- colSums(matrix_TGFbeta1_4day_batch1)
total_expression_matrix_TGFbeta1_8day_batch1 <- colSums(matrix_TGFbeta1_8day_batch1)

# Expressão mitocondrial por célula
mito_expression_per_cell_matrix_notreated_batch1 <- colSums(matrix_notreated_batch1[mito_genes_present_matrix_notreated_batch1, ])
mito_expression_per_cell_matrix_notreated_batch2 <- colSums(matrix_notreated_batch2[mito_genes_present_matrix_notreated_batch2, ])
mito_expression_per_cell_matrix_TGFbeta1_1day_batch2 <- colSums(matrix_TGFbeta1_1day_batch2[mito_genes_present_matrix_TGFbeta1_1day_batch2, ])
mito_expression_per_cell_matrix_TGFbeta1_2day_batch2 <- colSums(matrix_TGFbeta1_2day_batch2[mito_genes_present_matrix_TGFbeta1_2day_batch2, ])
mito_expression_per_cell_matrix_TGFbeta1_3day_batch2 <- colSums(matrix_TGFbeta1_3day_batch2[mito_genes_present_matrix_TGFbeta1_3day_batch2, ])
mito_expression_per_cell_matrix_TGFbeta1_4day_batch1 <- colSums(matrix_TGFbeta1_4day_batch1[mito_genes_present_matrix_TGFbeta1_4day_batch1, ])
mito_expression_per_cell_matrix_TGFbeta1_8day_batch1 <- colSums(matrix_TGFbeta1_8day_batch1[mito_genes_present_matrix_TGFbeta1_8day_batch1, ])

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


# Criar objeto Seurat
# objetoSeurat_teste <- CreateSeuratObject(counts = data_teste, project = "MCF10A", min.cells = 3, min.features = 200)

# MCF10A_notreated_batch1 <- CreateSeuratObject(counts = matrix_notreated_batch1, project = "MCF10A_notreated_batch1", min.cells = 3, min.features = 200)
# MCF10A_notreated_batch2 <- CreateSeuratObject(counts = matrix_notreated_batch2, project = "MCF10A_notreated_batch2", min.cells = 3, min.features = 200)
# MCF10A_TGFbeta1_1day_batch2 <- CreateSeuratObject(counts = matrix_TGFbeta1_1day_batch2, project = "MCF10A_TGFbeta1_1day_batch2", min.cells = 3, min.features = 200)
# MCF10A_TGFbeta1_2day_batch2 <- CreateSeuratObject(counts = matrix_TGFbeta1_2day_batch2, project = "MCF10A_TGFbeta1_2day_batch2", min.cells = 3, min.features = 200)
# MCF10A_TGFbeta1_3day_batch2 <- CreateSeuratObject(counts = matrix_TGFbeta1_3day_batch2, project = "MCF10A_TGFbeta1_3day_batch2", min.cells = 3, min.features = 200)
# MCF10A_TGFbeta1_4day_batch1 <- CreateSeuratObject(counts = matrix_TGFbeta1_4day_batch1, project = "MCF10A_TGFbeta1_4day_batch1", min.cells = 3, min.features = 200)
# MCF10A_TGFbeta1_8day_batch1 <- CreateSeuratObject(counts = matrix_TGFbeta1_8day_batch1, project = "MCF10A_TGFbeta1_8day_batch1", min.cells = 3, min.features = 200)

MCF10A_notreated_batch1 <- CreateSeuratObject(counts = expression_data_notreated_batch1, project = "MCF10A_notreated_batch1")
MCF10A_notreated_batch2 <- CreateSeuratObject(counts = expression_data_notreated_batch2, project = "MCF10A_notreated_batch2")
MCF10A_TGFbeta1_1day_batch2 <- CreateSeuratObject(counts = expression_data_TGFbeta1_1day_batch2, project = "MCF10A_TGFbeta1_1day_batch2")
MCF10A_TGFbeta1_2day_batch2 <- CreateSeuratObject(counts = expression_data_TGFbeta1_2day_batch2, project = "MCF10A_TGFbeta1_2day_batch2")
MCF10A_TGFbeta1_3day_batch2 <- CreateSeuratObject(counts = expression_data_TGFbeta1_3day_batch2, project = "MCF10A_TGFbeta1_3day_batch2")
MCF10A_TGFbeta1_4day_batch1 <- CreateSeuratObject(counts = expression_data_TGFbeta1_4day_batch1, project = "MCF10A_TGFbeta1_4day_batch1")
MCF10A_TGFbeta1_8day_batch1 <- CreateSeuratObject(counts = expression_data_TGFbeta1_8day_batch1, project = "MCF10A_TGFbeta1_8day_batch1")


# Calcular a porcentagem de genes mitocondriais

MCF10A_notreated_batch1[["percent.mt"]] <- PercentageFeatureSet(MCF10A_notreated_batch1, pattern = "^MT-")  # Para dados humanos (genes mitocondriais começam com "MT-")
MCF10A_notreated_batch2[["percent.mt"]] <- PercentageFeatureSet(MCF10A_notreated_batch2, pattern = "^MT-")
MCF10A_TGFbeta1_1day_batch2[["percent.mt"]] <- PercentageFeatureSet(MCF10A_TGFbeta1_1day_batch2, pattern = "^MT-")
MCF10A_TGFbeta1_2day_batch2[["percent.mt"]] <- PercentageFeatureSet(MCF10A_TGFbeta1_2day_batch2, pattern = "^MT-")
MCF10A_TGFbeta1_3day_batch2[["percent.mt"]] <- PercentageFeatureSet(MCF10A_TGFbeta1_3day_batch2, pattern = "^MT-")
MCF10A_TGFbeta1_4day_batch1[["percent.mt"]] <- PercentageFeatureSet(MCF10A_TGFbeta1_4day_batch1, pattern = "^MT-")
MCF10A_TGFbeta1_8day_batch1[["percent.mt"]] <- PercentageFeatureSet(MCF10A_TGFbeta1_8day_batch1, pattern = "^MT-")

# salvar_objetos_rdata("MCF10A_TGFbeta1_8day_batch1")

# Violin Plot das métricas de QC

library(Seurat)
library(patchwork)
library(ggplot2)
library(SeuratObject)

# save.image("AllObjectsSeuratQC.RData")

load("AllObjectsSeuratQC.RData")

# Combinar os objetos Seurat
combined_data <- merge(x = MCF10A_notreated_batch1, 
                       y = c(MCF10A_notreated_batch2, 
                             MCF10A_TGFbeta1_1day_batch2, 
                             MCF10A_TGFbeta1_2day_batch2, 
                             MCF10A_TGFbeta1_3day_batch2, 
                             MCF10A_TGFbeta1_4day_batch1, 
                             MCF10A_TGFbeta1_8day_batch1),
                       add.cell.ids = c("MCF10A_notreated_batch1",
                                        "MCF10A_notreated_batch2",
                                        "MCF10A_TGFbeta1_1day_batch2", 
                                        "MCF10A_TGFbeta1_2day_batch2",
                                        "MCF10A_TGFbeta1_3day_batch2",
                                        "MCF10A_TGFbeta1_4day_batch1",
                                        "MCF10A_TGFbeta1_8day_batch1"), 
                       project = "combined_data")

# Gerar o gráfico de violino para as métricas de QC
vln_plot <- VlnPlot(combined_data, 
                    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                    group.by = "orig.ident", 
                    ncol = 3) & 
                    theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

# Adicionar título e ajustar o tamanho da fonte
vln_plot <- vln_plot + 
  plot_annotation(title = "Gráficos de Violino para Métricas de QC") & 
  theme(plot.title = element_text(size = 14, hjust = 0.5))

# Exibir o gráfico
png(file = "images/GraficosdeViolinoparaMétricasDeQC.png", width = 1920, height = 1080, res = 100)
print(vln_plot)
dev.off()

# Salvar objeto, se necessário
salvar_objetos_rdata("vln_plot")


# Scatter Plots para identificar células problemáticas

FeatureScatter(MCF10A_notreated_batch1, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.5)
FeatureScatter(MCF10A_notreated_batch2, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.5)
FeatureScatter(MCF10A_TGFbeta1_1day_batch2, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.5)
FeatureScatter(MCF10A_TGFbeta1_2day_batch2, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.5)
FeatureScatter(MCF10A_TGFbeta1_3day_batch2, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.5)
FeatureScatter(MCF10A_TGFbeta1_4day_batch1, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.5)
FeatureScatter(MCF10A_TGFbeta1_8day_batch1, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.5)

png(file = "images/GraficosdeDispersãoCountMTFullDay.png", width = 1920, height = 1080, res = 100)
FeatureScatter(combined_data, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.5) + 
  plot_annotation(title = "Gráficos de Dispersão Count x Percent MT") & 
  theme(plot.title = element_text(size = 14, hjust = 0.5))
dev.off()

FeatureScatter(MCF10A_notreated_batch1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(MCF10A_notreated_batch2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(MCF10A_TGFbeta1_1day_batch2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(MCF10A_TGFbeta1_2day_batch2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(MCF10A_TGFbeta1_3day_batch2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(MCF10A_TGFbeta1_4day_batch1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(MCF10A_TGFbeta1_8day_batch1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

png(file = "images/GraficosdeDispersaoCountFeatureFullDay.png", width = 1920, height = 1080, res = 100)
FeatureScatter(combined_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.5) + 
  plot_annotation(title = "Gráficos de Dispersão Count x Feature RNA") & 
  theme(plot.title = element_text(size = 14, hjust = 0.5))
dev.off()


# Filtrar células com base nas métricas de QC
MCF10A_notreated_batch1 <- subset(MCF10A_notreated_batch1, subset = nFeature_RNA > 500 & percent.mt < 20)
MCF10A_notreated_batch2 <- subset(MCF10A_notreated_batch2, subset = nFeature_RNA > 500 & percent.mt < 20)
MCF10A_TGFbeta1_1day_batch2 <- subset(MCF10A_TGFbeta1_1day_batch2, subset = nFeature_RNA > 500 & percent.mt < 20)
MCF10A_TGFbeta1_2day_batch2 <- subset(MCF10A_TGFbeta1_2day_batch2, subset = nFeature_RNA > 500 & percent.mt < 20)
MCF10A_TGFbeta1_3day_batch2 <- subset(MCF10A_TGFbeta1_3day_batch2, subset = nFeature_RNA > 500 & percent.mt < 20)
MCF10A_TGFbeta1_4day_batch1 <- subset(MCF10A_TGFbeta1_4day_batch1, subset = nFeature_RNA > 500 & percent.mt < 20)
MCF10A_TGFbeta1_8day_batch1 <- subset(MCF10A_TGFbeta1_8day_batch1, subset = nFeature_RNA > 500 & percent.mt < 20)

# save.image("AllObjectsSeuratAfterQC.RData")
load("AllObjectsSeuratAfterQC.RData")

# Antes de normalizar com a função NormalizeData, a Prof Rita pediu que, após o QC,
# fosse realizada a soma de cada coluna para achar o total de count de cada célula.

## MCF10A_notreated_batch1
matrix_counts_notreated_batch1 <- MCF10A_notreated_batch1[["RNA"]]$counts
matrix_counts_notreated_batch1 <- as.data.frame(as.matrix(matrix_counts_notreated_batch1))
sum_notreated_batch1 <- colSums(matrix_counts_notreated_batch1)
matrix_counts_notreated_batch1_normalized <- sweep(matrix_counts_notreated_batch1, 2, sum_notreated_batch1, "/")
sum(colSums(matrix_counts_notreated_batch1_normalized) != 1)
# [1] 0

## MCF10A_notreated_batch2
matrix_counts_notreated_batch2 <- MCF10A_notreated_batch2[["RNA"]]$counts
matrix_counts_notreated_batch2 <- as.data.frame(as.matrix(matrix_counts_notreated_batch2))
sum_notreated_batch2 <- colSums(matrix_counts_notreated_batch2)
matrix_counts_notreated_batch2_normalized <- sweep(matrix_counts_notreated_batch2, 2, sum_notreated_batch2, "/")
sum(colSums(matrix_counts_notreated_batch2_normalized) != 1)
# [1] 10 // Descobrir o motivo de esar aparecendo 10. Deveria ser 0.

## MCF10A_TGFbeta1_1day_batch2
matrix_counts_TGFbeta1_1day_batch2 <- MCF10A_TGFbeta1_1day_batch2[["RNA"]]$counts
matrix_counts_TGFbeta1_1day_batch2 <- as.data.frame(as.matrix(matrix_counts_TGFbeta1_1day_batch2))
sum_TGFbeta1_1day_batch2 <- colSums(matrix_counts_TGFbeta1_1day_batch2)
matrix_counts_TGFbeta1_1day_batch2_normalized <- sweep(matrix_counts_TGFbeta1_1day_batch2, 2, sum_TGFbeta1_1day_batch2, "/")
sum(colSums(matrix_counts_TGFbeta1_1day_batch2_normalized) != 1)
# [1] 21 // Descobrir o motivo de esar aparecendo 10. Deveria ser 0.

## MCF10A_TGFbeta1_2day_batch2
matrix_counts_TGFbeta1_2day_batch2 <- MCF10A_TGFbeta1_2day_batch2[["RNA"]]$counts
matrix_counts_TGFbeta1_2day_batch2 <- as.data.frame(as.matrix(matrix_counts_TGFbeta1_2day_batch2))
sum_TGFbeta1_2day_batch2 <- colSums(matrix_counts_TGFbeta1_2day_batch2)
matrix_counts_TGFbeta1_2day_batch2_normalized <- sweep(matrix_counts_TGFbeta1_2day_batch2, 2, sum_TGFbeta1_2day_batch2, "/")
sum(colSums(matrix_counts_TGFbeta1_2day_batch2_normalized) != 1)
# [1] 2 // Descobrir o motivo de esar aparecendo 10. Deveria ser 0.

## MCF10A_TGFbeta1_3day_batch2
matrix_counts_TGFbeta1_3day_batch2 <- MCF10A_TGFbeta1_3day_batch2[["RNA"]]$counts
matrix_counts_TGFbeta1_3day_batch2 <- as.data.frame(as.matrix(matrix_counts_TGFbeta1_3day_batch2))
sum_TGFbeta1_3day_batch2 <- colSums(matrix_counts_TGFbeta1_3day_batch2)
matrix_counts_TGFbeta1_3day_batch2_normalized <- sweep(matrix_counts_TGFbeta1_3day_batch2, 2, sum_TGFbeta1_3day_batch2, "/")
sum(colSums(matrix_counts_TGFbeta1_3day_batch2_normalized) != 1)
# [1] 4 // Descobrir o motivo de esar aparecendo 10. Deveria ser 0.

## MCF10A_TGFbeta1_4day_batch1
matrix_counts_TGFbeta1_4day_batch1 <- MCF10A_TGFbeta1_4day_batch1[["RNA"]]$counts
matrix_counts_TGFbeta1_4day_batch1 <- as.data.frame(as.matrix(matrix_counts_TGFbeta1_4day_batch1))
sum_TGFbeta1_4day_batch1 <- colSums(matrix_counts_TGFbeta1_4day_batch1)
matrix_counts_TGFbeta1_4day_batch1_normalized <- sweep(matrix_counts_TGFbeta1_4day_batch1, 2, sum_TGFbeta1_4day_batch1, "/")
sum(colSums(matrix_counts_TGFbeta1_4day_batch1_normalized) != 1)
# [1] 0

## MCF10A_TGFbeta1_8day_batch1
matrix_counts_TGFbeta1_8day_batch1 <- MCF10A_TGFbeta1_8day_batch1[["RNA"]]$counts
matrix_counts_TGFbeta1_8day_batch1 <- as.data.frame(as.matrix(matrix_counts_TGFbeta1_8day_batch1))
sum_TGFbeta1_8day_batch1 <- colSums(matrix_counts_TGFbeta1_8day_batch1)
matrix_counts_TGFbeta1_8day_batch1_normalized <- sweep(matrix_counts_TGFbeta1_8day_batch1, 2, sum_TGFbeta1_8day_batch1, "/")
sum(0.9999 < colSums(matrix_counts_TGFbeta1_8day_batch1_normalized) < 1.0001)
 # [1] 0

### Criando as matrizes de contagem por meio do Pacote Matrix

## MCF10A_notreated_batch1

matrix_notreated_batch1 <- as.data.frame(as.matrix(matrix_notreated_batch1))
sum_notreated_batch1 <- colSums(matrix_notreated_batch1)
matrix_counts_notreated_batch1_normalized <- sweep(matrix_notreated_batch1, 2, sum_notreated_batch1, "/")
sum(colSums(matrix_counts_notreated_batch1_normalized) != 1)
# [1] 0

## MCF10A_notreated_batch2

matrix_notreated_batch2 <- as.data.frame(as.matrix(matrix_notreated_batch2))
sum_notreated_batch2 <- colSums(matrix_notreated_batch2)
matrix_counts_notreated_batch2_normalized <- sweep(matrix_notreated_batch2, 2, sum_notreated_batch2, "/")
# Verificar se a soma das colunas está fora do intervalo [0.999999, 1.000001]
sum(colSums(matrix_counts_notreated_batch2_normalized) < 0.999999 | colSums(matrix_counts_notreated_batch2_normalized) > 1.000001)


## MCF10A_TGFbeta1_1day_batch2

matrix_TGFbeta1_1day_batch2 <- as.data.frame(as.matrix(matrix_TGFbeta1_1day_batch2))
sum_TGFbeta1_1day_batch2 <- colSums(matrix_TGFbeta1_1day_batch2)
matrix_counts_TGFbeta1_1day_batch2_normalized <- sweep(matrix_TGFbeta1_1day_batch2, 2, sum_TGFbeta1_1day_batch2, "/")
# Verificar se a soma das colunas está fora do intervalo [0.999999, 1.000001]
sum(colSums(matrix_counts_TGFbeta1_1day_batch2_normalized) < 0.999999 | colSums(matrix_counts_TGFbeta1_1day_batch2_normalized) > 1.000001)

## MCF10A_TGFbeta1_2day_batch2

matrix_TGFbeta1_2day_batch2 <- as.data.frame(as.matrix(matrix_TGFbeta1_2day_batch2))
sum_TGFbeta1_2day_batch2 <- colSums(matrix_TGFbeta1_2day_batch2)
matrix_counts_TGFbeta1_2day_batch2_normalized <- sweep(matrix_TGFbeta1_2day_batch2, 2, sum_TGFbeta1_2day_batch2, "/")
# Verificar se a soma das colunas está fora do intervalo [0.999999, 1.000001]
sum(colSums(matrix_counts_TGFbeta1_2day_batch2_normalized) < 0.999999 | colSums(matrix_counts_TGFbeta1_2day_batch2_normalized) > 1.000001)

## MCF10A_TGFbeta1_3day_batch2

matrix_TGFbeta1_3day_batch2 <- as.data.frame(as.matrix(matrix_TGFbeta1_3day_batch2))
sum_TGFbeta1_3day_batch2 <- colSums(matrix_TGFbeta1_3day_batch2)
matrix_counts_TGFbeta1_3day_batch2_normalized <- sweep(matrix_TGFbeta1_3day_batch2, 2, sum_TGFbeta1_3day_batch2, "/")
# Verificar se a soma das colunas está fora do intervalo [0.999999, 1.000001]
sum(colSums(matrix_counts_TGFbeta1_3day_batch2_normalized) < 0.999999 | colSums(matrix_counts_TGFbeta1_3day_batch2_normalized) > 1.000001)

## MCF10A_TGFbeta1_4day_batch1

matrix_TGFbeta1_4day_batch1 <- as.data.frame(as.matrix(matrix_TGFbeta1_4day_batch1))
sum_TGFbeta1_4day_batch1 <- colSums(matrix_TGFbeta1_4day_batch1)
matrix_counts_TGFbeta1_4day_batch1_normalized <- sweep(matrix_TGFbeta1_4day_batch1, 2, sum_TGFbeta1_4day_batch1, "/")
# Verificar se a soma das colunas está fora do intervalo [0.999999, 1.000001]
sum(colSums(matrix_counts_TGFbeta1_4day_batch1_normalized) < 0.999999 | colSums(matrix_counts_TGFbeta1_4day_batch1_normalized) > 1.000001)

## MCF10A_TGFbeta1_8day_batch1

matrix_TGFbeta1_8day_batch1 <- as.data.frame(as.matrix(matrix_TGFbeta1_8day_batch1))
sum_TGFbeta1_8day_batch1 <- colSums(matrix_TGFbeta1_8day_batch1)
matrix_counts_TGFbeta1_8day_batch1_normalized <- sweep(matrix_TGFbeta1_8day_batch1, 2, sum_TGFbeta1_8day_batch1, "/")
# Verificar se a soma das colunas está fora do intervalo [0.999999, 1.000001]
sum(colSums(matrix_counts_TGFbeta1_8day_batch1_normalized) < 0.999999 | colSums(matrix_counts_TGFbeta1_8day_batch1_normalized) > 1.000001)

save.image("AllWithoutSeuratAfterNormalized.RData")

#Agora, com o NormalizeData

MCF10A_notreated_batch1 <- NormalizeData(MCF10A_notreated_batch1)
MCF10A_notreated_batch2 <- NormalizeData(MCF10A_notreated_batch2)
MCF10A_TGFbeta1_1day_batch2 <- NormalizeData(MCF10A_TGFbeta1_1day_batch2)
MCF10A_TGFbeta1_2day_batch2 <- NormalizeData(MCF10A_TGFbeta1_2day_batch2)
MCF10A_TGFbeta1_3day_batch2 <- NormalizeData(MCF10A_TGFbeta1_3day_batch2)
MCF10A_TGFbeta1_4day_batch1 <- NormalizeData(MCF10A_TGFbeta1_4day_batch1)
MCF10A_TGFbeta1_8day_batch1 <- NormalizeData(MCF10A_TGFbeta1_8day_batch1)

write.csv(matrix_counts_notreated_batch1_normalized, file = "matrix_counts_notreated_batch1_normalized.csv", row.names = FALSE)
write.csv(matrix_counts_notreated_batch2_normalized, file = "matrix_counts_notreated_batch2_normalized.csv", row.names = FALSE)
write.csv(matrix_counts_TGFbeta1_1day_batch2_normalized, file = "matrix_counts_TGFbeta1_1day_batch2_normalized.csv", row.names = FALSE)
write.csv(matrix_counts_TGFbeta1_2day_batch2_normalized, file = "matrix_counts_TGFbeta1_2day_batch2_normalized.csv", row.names = FALSE)
write.csv(matrix_counts_TGFbeta1_3day_batch2_normalized, file = "matrix_counts_TGFbeta1_3day_batch2_normalized.csv", row.names = FALSE)
write.csv(matrix_counts_TGFbeta1_4day_batch1_normalized, file = "matrix_counts_TGFbeta1_4day_batch1_normalized.csv", row.names = FALSE)
write.csv(matrix_counts_TGFbeta1_8day_batch1_normalized, file = "matrix_counts_TGFbeta1_8day_batch1_normalized.csv", row.names = FALSE)

# save.image("AllObjectsSeuratAfterNormalized.RData")
# load("AllObjectsSeuratAfterNormalized.RData")

# Remover nomes das colunas para desconsiderar a identificação
colnames(matrix_counts_notreated_batch1_normalized) <- NULL
colnames(matrix_counts_notreated_batch2_normalized) <- NULL

matrix_counts_notreated_batch1_normalized <- as.matrix(matrix_counts_notreated_batch1_normalized)
matrix_counts_notreated_batch2_normalized <- as.matrix(matrix_counts_notreated_batch2_normalized)

# Número de colunas das duas matrizes
num_cols_batch1 <- ncol(matrix_counts_notreated_batch1_normalized)
num_cols_batch2 <- ncol(matrix_counts_notreated_batch2_normalized)

# Média das células, coluna a coluna, até o limite de colunas da primeira matriz
matrix_counts_notreated_batch2_normalized_1.1342 <- matrix_counts_notreated_batch2_normalized[, 1:num_cols_batch1]

matrix_mean_notreated <- (matrix_counts_notreated_batch1_normalized + matrix_counts_notreated_batch2_normalized_1.1342) / 2

# Adicionar colunas adicionais da segunda matriz (sem média)
matrix_mean_notreated <- cbind(matrix_mean_notreated, matrix_counts_notreated_batch2_normalized[, (num_cols_batch1 + 1):num_cols_batch2])

# Exibir a matriz resultante
View(matrix_mean_notreated)


###### Juntando todas as 6 matrizes normalizadas em apenas uma.
# Primeiro, assegure-se de que todos os dataframes tenham o mesmo número de linhas
if (all(sapply(list(matrix_mean_notreated, 
                    matrix_counts_TGFbeta1_1day_batch2_normalized,
                    matrix_counts_TGFbeta1_2day_batch2_normalized,
                    matrix_counts_TGFbeta1_3day_batch2_normalized,
                    matrix_counts_TGFbeta1_4day_batch1_normalized,
                    matrix_counts_TGFbeta1_8day_batch1_normalized), 
               nrow) == nrow(matrix_mean_notreated))) {
  
  # Combine todos os dataframes por coluna
  combined_matrix <- cbind(matrix_mean_notreated,
                           matrix_counts_TGFbeta1_1day_batch2_normalized,
                           matrix_counts_TGFbeta1_2day_batch2_normalized,
                           matrix_counts_TGFbeta1_3day_batch2_normalized,
                           matrix_counts_TGFbeta1_4day_batch1_normalized,
                           matrix_counts_TGFbeta1_8day_batch1_normalized)
  
  # Renomeie as colunas
  colnames(combined_matrix) <- c(colnames(matrix_mean_notreated), "1day_batch2", "2day_batch2", "3day_batch2", "4day_batch1", "8day_batch1")
  
  # Salve como um arquivo CSV
  write.csv(combined_matrix, file = "combined_matrix.csv", row.names = TRUE)
  
  print("Arquivo combinado criado com sucesso!")
  
} else {
  stop("Os dataframes não têm o mesmo número de linhas.")
}

save.image("AllWithoutSeuratAfterNormalized.RData")

# 1. Verifique se todas as colunas são numéricas
# Para PCA, apenas colunas numéricas devem ser usadas
numeric_data <- combined_matrix[sapply(combined_matrix, is.numeric)]

# 2. Padronizar os dados (opcional, mas comum)
# A função scale() centraliza e escala os dados
scaled_data <- scale(numeric_data)

# 3. Realizar a PCA
pca_result <- prcomp(scaled_data, center = TRUE, scale. = TRUE)

# 4. Resumo dos resultados da PCA
summary(pca_result)

# 5. Visualizar a variação explicada pelas componentes principais
screeplot(pca_result, type = "lines", main = "Scree Plot", npcs = 20)

# Gráfico das duas primeiras componentes principais
biplot(pca_result)

# Alternativamente, se você quiser um gráfico mais bonito, pode usar ggplot2:
library(ggplot2)

# Criar um dataframe com os scores da PCA
pca_scores <- as.data.frame(pca_result$x)

# Adicionar os nomes das amostras se necessário
pca_scores$Sample <- rownames(pca_scores)

# Plotar as duas primeiras componentes principais
ggplot(pca_scores, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point() +
  geom_text(vjust = 1.5) + # Ajustar a posição do texto
  theme_minimal() +
  labs(title = "PCA: Primeiras duas componentes principais", x = "PC1", y = "PC2")










# save.image("AllAfterMeanNotreatedMatrix.RData")
write.csv(matrix_mean_notreated, file = "matrix_mean_notreated_normalized.csv", row.names = FALSE)

# Converta a matriz resultante de volta para um data frame, se necessário
matrix_mean_notreated <- as.data.frame(matrix_mean_notreated)

# Visualizando a nova matriz com médias das não tratadas.

View(matrix_mean_notreated)

### FAZER TRANSCRIPTOGRAMAS DOS SEIS DIAS #####################################

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
dictionary <- getBM(attributes = c("ensembl_peptide_id", "ensembl_gene_id"),
                    mart = ensembl)
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

# Crie um dicionário único para cada ENSG -> ENSP usando as colunas corretas
unique_dictionary <- dictionary[!duplicated(dictionary$ensembl_gene_id), c("ensembl_peptide_id", "ensembl_gene_id")]

# Verifique a interseção entre os `rownames` de `combined_matrix` e `dictionary`
identifiers_expression <- rownames(combined_matrix)
identifiers_dictionary <- unique_dictionary$ensembl_gene_id

# Identificadores que estão em comum
common_identifiers <- intersect(identifiers_expression, identifiers_dictionary)
cat("Número de identificadores em comum:", length(common_identifiers), "\n")

# Filtre `combined_matrix` para manter apenas os identificadores em comum
combined_matrix_gene_common <- combined_matrix[rownames(combined_matrix) %in% common_identifiers, ]

# Filtre `unique_dictionary` para manter apenas os identificadores em comum
unique_dictionary <- unique_dictionary[unique_dictionary$ensembl_gene_id %in% common_identifiers, ]

# Crie um mapeamento ENSG -> ENSP e substitua os `rownames` da matriz
rownames_mapping <- setNames(unique_dictionary$ensembl_peptide_id, unique_dictionary$ensembl_gene_id)

# Substituir os rownames da matriz `combined_matrix` para ENSP
rownames(combined_matrix_gene_common) <- rownames_mapping[rownames(combined_matrix_gene_common)]

# Run transcriptogram preprocess, radius = 30
t_combined <- transcriptogramPreprocess(association = assoc, ordering = ord$Protein, radius = 30)
t_combined <- transcriptogramStep1(object = t_combined, expression = combined_matrix, dictionary = dictionary)
t_combined <- transcriptogramStep2(object = t_combined)

save.image("AllWithT_Combined.RData")

write.csv(t_combined@transcriptogramS2, file = "t_combined.csv", row.names = T)

# Extraindo o dataframe do slot transcriptogramS2
df <- t_combined@transcriptogramS2

# Verificando o tipo de cada coluna
sapply(df, class)  # Para identificar as colunas que não são numéricas

# Mantendo apenas colunas numéricas
df_numeric <- df[, sapply(df, is.numeric)]

# Verificando se sobrou algum dado
if (ncol(df_numeric) == 0) {
  stop("Nenhuma coluna numérica encontrada no dataframe!")
}

# Realizando a padronização dos dados numéricos
df_scaled <- scale(df_numeric)

# Continuando com o PCA
pca_result <- prcomp(df_scaled, center = TRUE, scale. = TRUE)

# Resumo dos resultados
summary(pca_result)

save.image("AfterPCAofT_Combined.RData")

# Scree plot da variância explicada
plot(pca_result, type = "l", main = "Scree Plot - PCA")

# Gráfico da variância acumulada
variancia_explicada <- pca_result$sdev^2 / sum(pca_result$sdev^2)
variancia_acumulada <- cumsum(variancia_explicada)
plot(variancia_acumulada, type = "b", pch = 19, col = "blue",
     xlab = "Número de Componentes Principais", 
     ylab = "Variância Acumulada",
     main = "Gráfico de Variância Acumulada")
abline(h = 0.9, col = "red", lty = 2)  # Linha para 90% da variância expli



# Criando o gráfico combinado
plot(variancia_explicada, type = "b", pch = 19, col = "blue", 
     xlab = "Número de Componentes Principais", 
     ylab = "Variância Explicada / Acumulada",
     ylim = c(0, 1), main = "Scree Plot e Variância Acumulada")
lines(variancia_acumulada, type = "b", pch = 19, col = "red")
legend("bottomright", legend = c("Variância Explicada", "Variância Acumulada"), 
       col = c("blue", "red"), pch = 19, lty = 1)
abline(h = 0.9, col = "darkgreen", lty = 2)  # Linha para 90% da variância acumulada


# Número de componentes principais a exibir (10 ou 20)
num_componentes <- 10  # Altere para 20, se preferir

# Limitando os dados ao número de componentes principais desejados
variancia_explicada_limited <- variancia_explicada[1:num_componentes]
variancia_acumulada_limited <- variancia_acumulada[1:num_componentes]

png("images/PCAeVarianciaAcumuladaTop10TranscriptogramaCombinadoSeisDias.png", width = 1920, height = 1080, res = 200)
# Criando o gráfico combinado
plot(variancia_explicada_limited, type = "b", pch = 19, col = "blue", 
     xlab = "Número de Componentes Principais", 
     ylab = "Variância Acumulada / PC",
     ylim = c(0, 1), main = "PCA e Variância Acumulada - Top 10 \n Transcriptograma Combinado dos 6 Dias",
     xaxt = "n")
axis(1, at = 1:num_componentes, labels = 1:num_componentes)  # Ajusta o eixo X
lines(variancia_acumulada_limited, type = "b", pch = 19, col = "red")
legend("bottomright", legend = c("PC", "Variância Acumulada"), 
       col = c("blue", "red"), pch = 19, lty = 1)
abline(h = 0.9, col = "darkgreen", lty = 2)  # Linha para 90% da variância acumulada
dev.off()

t_combined_r0 <- transcriptogramPreprocess(association = assoc, ordering = ord$Protein, radius = 0)
t_combined_r0 <- transcriptogramStep1(object = t_combined_r0, expression = combined_matrix, dictionary = dictionary)
t_combined_r0 <- transcriptogramStep2(object = t_combined_r0)

write.csv(t_combined_r0@transcriptogramS2, file = "t_combined_r0.csv", row.names = T)

# Extraindo o dataframe do slot transcriptogramS2
df_t_combined_r0 <- t_combined_r0@transcriptogramS2

# Verificando o tipo de cada coluna
sapply(df_t_combined_r0, class)  # Para identificar as colunas que não são numéricas

# Mantendo apenas colunas numéricas
df_t_combined_r0_numeric <- df_t_combined_r0[, sapply(df_t_combined_r0, is.numeric)]

# Verificando se sobrou algum dado
if (ncol(df_t_combined_r0_numeric) == 0) {
  stop("Nenhuma coluna numérica encontrada no dataframe R = 0!")
}

# Realizando a padronização dos dados numéricos
df_t_combined_r0_scaled <- scale(df_t_combined_r0_numeric)

# Continuando com o PCA
pca_df_t_combined_r0 <- prcomp(df_t_combined_r0_scaled, center = TRUE, scale. = TRUE)

# Resumo dos resultados
summary(pca_df_t_combined_r0)

save.image("AfterPCAofT_Combined_R0.RData")

# Scree plot da variância explicada
plot(pca_df_t_combined_r0, type = "l", main = "Scree Plot - PCA")


# Número de componentes principais a exibir (10 ou 20)
num_componentes <- 10  # Altere para 20, se preferir

# Limitando os dados ao número de componentes principais desejados
variancia_explicada_limited <- variancia_explicada[1:num_componentes]
variancia_acumulada_limited <- variancia_acumulada[1:num_componentes]

png("images/PCAeVarianciaAcumuladaTop10TranscriptogramaCombinadoSeisDiasR0.png", width = 1920, height = 1080, res = 200)
# Criando o gráfico combinado
plot(variancia_explicada_limited, type = "b", pch = 19, col = "blue", 
     xlab = "Número de Componentes Principais", 
     ylab = "Variância Acumulada / PC",
     ylim = c(0, 1), main = "PCA e Variância Acumulada - Top 10 \n Transcriptograma Combinado dos 6 Dias - R = 0",
     xaxt = "n")
axis(1, at = 1:num_componentes, labels = 1:num_componentes)  # Ajusta o eixo X
lines(variancia_acumulada_limited, type = "b", pch = 19, col = "red")
legend("bottomright", legend = c("PC", "Variância Acumulada"), 
       col = c("blue", "red"), pch = 19, lty = 1)
abline(h = 0.9, col = "darkgreen", lty = 2)  # Linha para 90% da variância acumulada
dev.off()

     



t_TGFbeta1_1day_batch2 <- transcriptogramPreprocess(association = assoc, ordering = ord$Protein, radius = 30)
t_TGFbeta1_1day_batch2 <- transcriptogramStep1(object = t_TGFbeta1_1day_batch2, expression = matrix_counts_TGFbeta1_1day_batch2_normalized, dictionary = dictionary)
t_TGFbeta1_1day_batch2 <- transcriptogramStep2(object = t_TGFbeta1_1day_batch2)

write.csv(t_TGFbeta1_1day_batch2@transcriptogramS2, file = "t_TGFbeta1_1day_batch2.csv", row.names = FALSE)

t_TGFbeta1_2day_batch2 <- transcriptogramPreprocess(association = assoc, ordering = ord$Protein, radius = 30)
t_TGFbeta1_2day_batch2 <- transcriptogramStep1(object = t_TGFbeta1_2day_batch2, expression = matrix_counts_TGFbeta1_2day_batch2_normalized, dictionary = dictionary)
t_TGFbeta1_2day_batch2 <- transcriptogramStep2(object = t_TGFbeta1_2day_batch2)

write.csv(t_TGFbeta1_2day_batch2@transcriptogramS2, file = "t_TGFbeta1_2day_batch2.csv", row.names = FALSE)

t_TGFbeta1_3day_batch2 <- transcriptogramPreprocess(association = assoc, ordering = ord$Protein, radius = 30)
t_TGFbeta1_3day_batch2 <- transcriptogramStep1(object = t_TGFbeta1_3day_batch2, expression = matrix_counts_TGFbeta1_3day_batch2_normalized, dictionary = dictionary)
t_TGFbeta1_3day_batch2 <- transcriptogramStep2(object = t_TGFbeta1_3day_batch2)

write.csv(t_TGFbeta1_3day_batch2@transcriptogramS2, file = "t_TGFbeta1_3day_batch2.csv", row.names = FALSE)

t_TGFbeta1_4day_batch1 <- transcriptogramPreprocess(association = assoc, ordering = ord$Protein, radius = 30)
t_TGFbeta1_4day_batch1 <- transcriptogramStep1(object = t_TGFbeta1_4day_batch1, expression = matrix_counts_TGFbeta1_4day_batch1_normalized, dictionary = dictionary)
t_TGFbeta1_4day_batch1 <- transcriptogramStep2(object = t_TGFbeta1_4day_batch1)

write.csv(t_TGFbeta1_4day_batch1@transcriptogramS2, file = "t_TGFbeta1_4day_batch1.csv", row.names = FALSE)

t_TGFbeta1_8day_batch1 <- transcriptogramPreprocess(association = assoc, ordering = ord$Protein, radius = 30)
t_TGFbeta1_8day_batch1 <- transcriptogramStep1(object = t_TGFbeta1_8day_batch1, expression = matrix_counts_TGFbeta1_8day_batch1_normalized, dictionary = dictionary)
t_TGFbeta1_8day_batch1 <- transcriptogramStep2(object = t_TGFbeta1_8day_batch1)

write.csv(t_TGFbeta1_8day_batch1@transcriptogramS2, file = "t_TGFbeta1_8day_batch1.csv", row.names = FALSE)

save.image("AllAfterTranscriptogramers.RData")
###############################################################################

####Fazer um arquivo com todas as células (transcriptograma ou matrix de expressão) de todas as 6 classes.

t_combined_all_days <- data.frame(
  t_notreated_mean = t_notreated_mean@transcriptogramS2,
  t_TGFbeta1_1day_batch2 = t_TGFbeta1_1day_batch2@transcriptogramS2,
  t_TGFbeta1_2day_batch2 = t_TGFbeta1_2day_batch2@transcriptogramS2,
  t_TGFbeta1_3day_batch2 = t_TGFbeta1_3day_batch2@transcriptogramS2,
  t_TGFbeta1_4day_batch1 = t_TGFbeta1_4day_batch1@transcriptogramS2,
  t_TGFbeta1_8day_batch1 = t_TGFbeta1_8day_batch1@transcriptogramS2
)

# Salve o arquivo em CSV
write.csv(t_combined_all_days, "t_combined_all_days.csv", row.names = TRUE)

library(irlba)
### Fazendo PCA do t_combined_all_days

# Verifique se há valores ausentes e, se necessário, remova ou substitua-os
t_combined_all_days[is.na(t_combined_all_days)] <- 0  # Opcional: substitui NA por zero

# Verifique se todas as colunas são numéricas
t_combined_all_days <- t_combined_all_days[sapply(t_combined_all_days, is.numeric)]

# Padronize os dados, centralizando e escalando as colunas
# t_combined_all_days_scaled <- scale(numeric_data)

# Realize a PCA

# Selecionando as 50 primeiras componentes principais
pca_result <- prcomp_irlba(t_combined_all_days, n = 50, center = TRUE, scale. = TRUE)

# pca_result <- prcomp(t_combined_all_days_scaled, center = TRUE, scale. = TRUE)
# pca_result <- prcomp(numeric_data, center = TRUE, scale. = FALSE)

# Veja um resumo dos resultados
summary(pca_result)

# Visualize os componentes principais
# Para ver a variação explicada pelas primeiras componentes principais
screeplot(pca_result, type = "lines", main = "Scree Plot")

# Visualizar as duas primeiras componentes principais
library(ggplot2)
pca_data <- as.data.frame(pca_result$x)  # Extrai as coordenadas da PCA
ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point() +
  ggtitle("PCA Plot das duas primeiras componentes") +
  xlab("PC1") +
  ylab("PC2")

# Extrair variância explicada por cada componente
variancia_explicada <- pca_result$sdev^2 / sum(pca_result$sdev^2)
variancia_acumulada <- cumsum(variancia_explicada)

# Criar um dataframe para facilitar o plot
scree_data <- data.frame(
  Component = seq_along(variancia_explicada),
  VarianciaExplicada = variancia_explicada,
  VarianciaAcumulada = variancia_acumulada)

png(file = "images/ComponentesPrincipais+VarianciaCumulativaTop10.png", width = 1920, height = 1080, res = 200)

# Gráfico Scree Plot usando ggplot2
ggplot(scree_data[1:10, ], aes(x = Component)) +
  geom_bar(aes(y = VarianciaExplicada), stat = "identity", fill = "steelblue") +
  geom_line(aes(y = VarianciaAcumulada), color = "red", size = 0.5) +
  geom_point(aes(y = VarianciaAcumulada), color = "red") +
  geom_text(aes(y = VarianciaExplicada, label = scales::percent(VarianciaExplicada)),
            vjust = -0.5, size = 3) +
  labs(title = "Scree Plot - Variância Explicada por Componentes Principais",
       x = "Componente Principal",
       y = "Variância Explicada") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14)) +
  scale_y_continuous(labels = scales::percent)
dev.off()


save.image("AfterPCAsAllDays.RData")






# Identificar genes variáveis
MCF10A_notreated_batch1 <- FindVariableFeatures(MCF10A_notreated_batch1, selection.method = "vst", nfeatures = 2000)
MCF10A_notreated_batch2 <- FindVariableFeatures(MCF10A_notreated_batch2, selection.method = "vst", nfeatures = 2000)
MCF10A_TGFbeta1_1day_batch2 <- FindVariableFeatures(MCF10A_TGFbeta1_1day_batch2, selection.method = "vst", nfeatures = 2000)
MCF10A_TGFbeta1_2day_batch2 <- FindVariableFeatures(MCF10A_TGFbeta1_2day_batch2, selection.method = "vst", nfeatures = 2000)
MCF10A_TGFbeta1_3day_batch2 <- FindVariableFeatures(MCF10A_TGFbeta1_3day_batch2, selection.method = "vst", nfeatures = 2000)
MCF10A_TGFbeta1_4day_batch1 <- FindVariableFeatures(MCF10A_TGFbeta1_4day_batch1, selection.method = "vst", nfeatures = 2000)
MCF10A_TGFbeta1_8day_batch1 <- FindVariableFeatures(MCF10A_TGFbeta1_8day_batch1, selection.method = "vst", nfeatures = 2000)

# Selecionar as 2000 features variáveis mais comuns entre os datasets
features <- SelectIntegrationFeatures(object.list = list(MCF10A_notreated_batch1, 
                                                         MCF10A_notreated_batch2,
                                                         MCF10A_TGFbeta1_1day_batch2,
                                                         MCF10A_TGFbeta1_2day_batch2,
                                                         MCF10A_TGFbeta1_3day_batch2,
                                                         MCF10A_TGFbeta1_4day_batch1,
                                                         MCF10A_TGFbeta1_8day_batch1))


# Preparar os datasets para integração encontrando âncoras
anchors <- FindIntegrationAnchors(object.list = list(MCF10A_notreated_batch1, 
                                                     MCF10A_notreated_batch2,
                                                     MCF10A_TGFbeta1_1day_batch2,
                                                     MCF10A_TGFbeta1_2day_batch2,
                                                     MCF10A_TGFbeta1_3day_batch2,
                                                     MCF10A_TGFbeta1_4day_batch1,
                                                     MCF10A_TGFbeta1_8day_batch1), 
                                  anchor.features = features)
save.image("AntesDeIntegrar.RData")

# Integração dos datasets
integrated_data <- IntegrateData(anchorset = anchors)

# Escalar os dados e realizar PCA
integrated_data <- ScaleData(integrated_data) %>% RunPCA()

# UMAP para visualização
integrated_data <- RunUMAP(integrated_data, dims = 1:30)

# Agrupamento (clustering)
integrated_data <- FindNeighbors(integrated_data, dims = 1:30)
integrated_data <- FindClusters(integrated_data, resolution = 0.5)

# Visualização UMAP com clusters
png(file = "images/ClusterVisualizationUMAPDataIntegrated.png", width = 1920, height = 1080, res = 200)
DimPlot(integrated_data, reduction = "umap", group.by = "seurat_clusters") +
  ggtitle("Cluster Visualization UMAP Data Integrated")
dev.off()

# Visualizar batches no UMAP
png(file = "images/BatchsVisualizationUMAPDataIntegrated.png", width = 1920, height = 1080, res = 200)
DimPlot(integrated_data, reduction = "umap", group.by = "orig.ident") +
  ggtitle("Batchs Visualization UMAP Data Integrated")
dev.off()

# Integrar em dois grupos por batchs.

# Selecionar as 2000 features variáveis mais comuns entre os datasets
features_batch1 <- SelectIntegrationFeatures(object.list = list(MCF10A_notreated_batch1, 
                                                                MCF10A_TGFbeta1_4day_batch1,
                                                                MCF10A_TGFbeta1_8day_batch1))

# Selecionar as 2000 features variáveis mais comuns entre os datasets
features_batch2 <- SelectIntegrationFeatures(object.list = list(MCF10A_notreated_batch2,
                                                                MCF10A_TGFbeta1_1day_batch2,
                                                                MCF10A_TGFbeta1_2day_batch2,
                                                                MCF10A_TGFbeta1_3day_batch2))


# Preparar os datasets para integração encontrando âncoras
anchors_batch1 <- FindIntegrationAnchors(object.list = list(MCF10A_notreated_batch1, 
                                                            MCF10A_TGFbeta1_4day_batch1,
                                                            MCF10A_TGFbeta1_8day_batch1), 
                                         anchor.features = features_batch1)

anchors_batch2 <- FindIntegrationAnchors(object.list = list(MCF10A_notreated_batch2,
                                                            MCF10A_TGFbeta1_1day_batch2,
                                                            MCF10A_TGFbeta1_2day_batch2,
                                                            MCF10A_TGFbeta1_3day_batch2), 
                                         anchor.features = features_batch2)

# Integração dos datasets
integrated_data_batch1 <- IntegrateData(anchorset = anchors_batch1)

integrated_data_batch2 <- IntegrateData(anchorset = anchors_batch2)

# Escalar os dados e realizar PCA
integrated_data_batch1 <- ScaleData(integrated_data_batch1) %>% RunPCA()

integrated_data_batch2 <- ScaleData(integrated_data_batch2) %>% RunPCA()

# UMAP para visualização
integrated_data_batch1 <- RunUMAP(integrated_data_batch1, dims = 1:30)
integrated_data_batch2 <- RunUMAP(integrated_data_batch2, dims = 1:30)

# Agrupamento (clustering)
integrated_data_batch1 <- FindNeighbors(integrated_data_batch1, dims = 1:30)
integrated_data_batch2 <- FindNeighbors(integrated_data_batch2, dims = 1:30)

integrated_data_batch1 <- FindClusters(integrated_data_batch1, resolution = 0.5)
integrated_data_batch2 <- FindClusters(integrated_data_batch2, resolution = 0.5)

# Visualização UMAP com clusters
png(file = "images/ClusterVisualizationUMAPDataIntegrated_batch1.png", width = 1920, height = 1080, res = 200)
DimPlot(integrated_data_batch1, reduction = "umap", group.by = "seurat_clusters") +
  ggtitle("Cluster Visualization UMAP Data Integrated - Batch 1")
dev.off()

png(file = "images/ClusterVisualizationUMAPDataIntegrated_batch2.png", width = 1920, height = 1080, res = 200)
DimPlot(integrated_data_batch2, reduction = "umap", group.by = "seurat_clusters") +
  ggtitle("Cluster Visualization UMAP Data Integrated - Batch 2")
dev.off()

# Visualizar batches no UMAP
png(file = "images/BatchsVisualizationUMAPDataIntegrated_batch1.png", width = 1920, height = 1080, res = 200)
DimPlot(integrated_data_batch1, reduction = "umap", group.by = "orig.ident") +
  ggtitle("Batchs Visualization UMAP Data Integrated - Batch 1")
dev.off()

png(file = "images/BatchsVisualizationUMAPDataIntegrated_batch2.png", width = 1920, height = 1080, res = 200)
DimPlot(integrated_data_batch2, reduction = "umap", group.by = "orig.ident") +
  ggtitle("Batchs Visualization UMAP Data Integrated - Batch 2")
dev.off()

# Escalonar os dados SEM INTEGRAR
MCF10A_notreated_batch1 <- ScaleData(MCF10A_notreated_batch1, features = rownames(MCF10A_notreated_batch1))
MCF10A_notreated_batch2 <- ScaleData(MCF10A_notreated_batch2, features = rownames(MCF10A_notreated_batch2))
MCF10A_TGFbeta1_1day_batch2 <- ScaleData(MCF10A_TGFbeta1_1day_batch2, features = rownames(MCF10A_TGFbeta1_1day_batch2))
MCF10A_TGFbeta1_2day_batch2 <- ScaleData(MCF10A_TGFbeta1_2day_batch2, features = rownames(MCF10A_TGFbeta1_2day_batch2))
MCF10A_TGFbeta1_3day_batch2 <- ScaleData(MCF10A_TGFbeta1_3day_batch2, features = rownames(MCF10A_TGFbeta1_3day_batch2))
MCF10A_TGFbeta1_4day_batch1 <- ScaleData(MCF10A_TGFbeta1_4day_batch1, features = rownames(MCF10A_TGFbeta1_4day_batch1))
MCF10A_TGFbeta1_8day_batch1 <- ScaleData(MCF10A_TGFbeta1_8day_batch1, features = rownames(MCF10A_TGFbeta1_8day_batch1))

# Executar PCA
MCF10A_notreated_batch1 <- RunPCA(MCF10A_notreated_batch1, features = VariableFeatures(object = MCF10A_notreated_batch1))
MCF10A_notreated_batch2 <- RunPCA(MCF10A_notreated_batch2, features = VariableFeatures(object = MCF10A_notreated_batch2))
MCF10A_TGFbeta1_1day_batch2 <- RunPCA(MCF10A_TGFbeta1_1day_batch2, features = VariableFeatures(object = MCF10A_TGFbeta1_1day_batch2))
MCF10A_TGFbeta1_2day_batch2 <- RunPCA(MCF10A_TGFbeta1_2day_batch2, features = VariableFeatures(object = MCF10A_TGFbeta1_2day_batch2))
MCF10A_TGFbeta1_3day_batch2 <- RunPCA(MCF10A_TGFbeta1_3day_batch2, features = VariableFeatures(object = MCF10A_TGFbeta1_3day_batch2))
MCF10A_TGFbeta1_4day_batch1 <- RunPCA(MCF10A_TGFbeta1_4day_batch1, features = VariableFeatures(object = MCF10A_TGFbeta1_4day_batch1))
MCF10A_TGFbeta1_8day_batch1 <- RunPCA(MCF10A_TGFbeta1_8day_batch1, features = VariableFeatures(object = MCF10A_TGFbeta1_8day_batch1))

# Visualizar os resultados do PCA
ElbowPlot(MCF10A_notreated_batch1)
ElbowPlot(MCF10A_notreated_batch2)
ElbowPlot(MCF10A_TGFbeta1_1day_batch2)
ElbowPlot(MCF10A_TGFbeta1_2day_batch2)
ElbowPlot(MCF10A_TGFbeta1_3day_batch2)
ElbowPlot(MCF10A_TGFbeta1_4day_batch1)
ElbowPlot(MCF10A_TGFbeta1_8day_batch1)

# Clustering
MCF10A_notreated_batch1 <- FindNeighbors(MCF10A_notreated_batch1, dims = 1:10)
MCF10A_notreated_batch2 <- FindNeighbors(MCF10A_notreated_batch2, dims = 1:10)
MCF10A_TGFbeta1_1day_batch2 <- FindNeighbors(MCF10A_TGFbeta1_1day_batch2, dims = 1:10)
MCF10A_TGFbeta1_2day_batch2 <- FindNeighbors(MCF10A_TGFbeta1_2day_batch2, dims = 1:10)
MCF10A_TGFbeta1_3day_batch2 <- FindNeighbors(MCF10A_TGFbeta1_3day_batch2, dims = 1:10)
MCF10A_TGFbeta1_4day_batch1 <- FindNeighbors(MCF10A_TGFbeta1_4day_batch1, dims = 1:10)
MCF10A_TGFbeta1_8day_batch1 <- FindNeighbors(MCF10A_TGFbeta1_8day_batch1, dims = 1:10)


MCF10A_notreated_batch1 <- FindClusters(MCF10A_notreated_batch1, resolution = 0.5)
MCF10A_notreated_batch2 <- FindClusters(MCF10A_notreated_batch2, resolution = 0.5)
MCF10A_TGFbeta1_1day_batch2 <- FindClusters(MCF10A_TGFbeta1_1day_batch2, resolution = 0.5)
MCF10A_TGFbeta1_2day_batch2 <- FindClusters(MCF10A_TGFbeta1_2day_batch2, resolution = 0.5)
MCF10A_TGFbeta1_3day_batch2 <- FindClusters(MCF10A_TGFbeta1_3day_batch2, resolution = 0.5)
MCF10A_TGFbeta1_4day_batch1 <- FindClusters(MCF10A_TGFbeta1_4day_batch1, resolution = 0.5)
MCF10A_TGFbeta1_8day_batch1 <- FindClusters(MCF10A_TGFbeta1_8day_batch1, resolution = 0.5)


# UMAP para visualização
MCF10A_notreated_batch1 <- RunUMAP(MCF10A_notreated_batch1, dims = 1:10)
MCF10A_notreated_batch2 <- RunUMAP(MCF10A_notreated_batch2, dims = 1:10)
MCF10A_TGFbeta1_1day_batch2 <- RunUMAP(MCF10A_TGFbeta1_1day_batch2, dims = 1:10)
MCF10A_TGFbeta1_2day_batch2 <- RunUMAP(MCF10A_TGFbeta1_2day_batch2, dims = 1:10)
MCF10A_TGFbeta1_3day_batch2 <- RunUMAP(MCF10A_TGFbeta1_3day_batch2, dims = 1:10)
MCF10A_TGFbeta1_4day_batch1 <- RunUMAP(MCF10A_TGFbeta1_4day_batch1, dims = 1:10)
MCF10A_TGFbeta1_8day_batch1 <- RunUMAP(MCF10A_TGFbeta1_8day_batch1, dims = 1:10)


DimPlot(MCF10A_notreated_batch1, reduction = "umap")
DimPlot(MCF10A_notreated_batch1, reduction = "pca")

DimPlot(MCF10A_notreated_batch2, reduction = "umap")
DimPlot(MCF10A_notreated_batch2, reduction = "pca")

DimPlot(MCF10A_TGFbeta1_1day_batch2, reduction = "umap")
DimPlot(MCF10A_TGFbeta1_1day_batch2, reduction = "pca")

DimPlot(MCF10A_TGFbeta1_2day_batch2, reduction = "umap")
DimPlot(MCF10A_TGFbeta1_2day_batch2, reduction = "pca")

DimPlot(MCF10A_TGFbeta1_3day_batch2, reduction = "umap")
DimPlot(MCF10A_TGFbeta1_3day_batch2, reduction = "pca")

DimPlot(MCF10A_TGFbeta1_4day_batch1, reduction = "umap")
DimPlot(MCF10A_TGFbeta1_4day_batch1, reduction = "pca")

DimPlot(MCF10A_TGFbeta1_8day_batch1, reduction = "umap")
DimPlot(MCF10A_TGFbeta1_8day_batch1, reduction = "pca")

# 2. Extração da Matriz de Expressão

# Extrair matriz de expressão
expression_matrix_MCF10A_notreated_batch1 <- as.data.frame(GetAssayData(MCF10A_notreated_batch1, layer = "data", assay = "RNA"))
expression_matrix_MCF10A_notreated_batch2 <- as.data.frame(GetAssayData(MCF10A_notreated_batch2, layer = "data", assay = "RNA"))
expression_matrix_MCF10A_TGFbeta1_1day_batch2 <- as.data.frame(GetAssayData(MCF10A_TGFbeta1_1day_batch2, layer = "data", assay = "RNA"))
expression_matrix_MCF10A_TGFbeta1_2day_batch2 <- as.data.frame(GetAssayData(MCF10A_TGFbeta1_2day_batch2, layer = "data", assay = "RNA"))
expression_matrix_MCF10A_TGFbeta1_3day_batch2 <- as.data.frame(GetAssayData(MCF10A_TGFbeta1_3day_batch2, layer = "data", assay = "RNA"))
expression_matrix_MCF10A_TGFbeta1_4day_batch1 <- as.data.frame(GetAssayData(MCF10A_TGFbeta1_4day_batch1, layer = "data", assay = "RNA"))
expression_matrix_MCF10A_TGFbeta1_8day_batch1 <- as.data.frame(GetAssayData(MCF10A_TGFbeta1_8day_batch1, layer = "data", assay = "RNA"))



expression_matrix_MCF10A_notreated_batch1$RowNames <- rownames(expression_matrix_MCF10A_notreated_batch1)
expression_matrix_MCF10A_notreated_batch2$RowNames <- rownames(expression_matrix_MCF10A_notreated_batch2)
expression_matrix_MCF10A_TGFbeta1_1day_batch2$RowNames <- rownames(expression_matrix_MCF10A_TGFbeta1_1day_batch2)
expression_matrix_MCF10A_TGFbeta1_2day_batch2$RowNames <- rownames(expression_matrix_MCF10A_TGFbeta1_2day_batch2)
expression_matrix_MCF10A_TGFbeta1_3day_batch2$RowNames <- rownames(expression_matrix_MCF10A_TGFbeta1_3day_batch2)
expression_matrix_MCF10A_TGFbeta1_4day_batch1$RowNames <- rownames(expression_matrix_MCF10A_TGFbeta1_4day_batch1)
expression_matrix_MCF10A_TGFbeta1_8day_batch1$RowNames <- rownames(expression_matrix_MCF10A_TGFbeta1_8day_batch1)

expression_matrix_MCF10A <- full_join(expression_matrix_MCF10A_notreated_batch1, expression_matrix_MCF10A_notreated_batch2, by = "RowNames") %>%
  full_join(expression_matrix_MCF10A_TGFbeta1_1day_batch2, by = "RowNames") %>%
  full_join(expression_matrix_MCF10A_TGFbeta1_2day_batch2, by = "RowNames") %>%
  full_join(expression_matrix_MCF10A_TGFbeta1_3day_batch2, by = "RowNames") %>%
  full_join(expression_matrix_MCF10A_TGFbeta1_4day_batch1, by = "RowNames") %>%
  full_join(expression_matrix_MCF10A_TGFbeta1_8day_batch1, by = "RowNames")

# Preenchendo os NAs com 0
expression_matrix_MCF10A[is.na(expression_matrix_MCF10A)] <- 0



# Movendo os nomes das linhas de volta para rownames e removendo a coluna auxiliar
rownames(expression_matrix_MCF10A) <- expression_matrix_MCF10A$RowNames
expression_matrix_MCF10A$RowNames <- NULL

expression_matrix_MCF10A_notreated_batch1$RowNames <- NULL
expression_matrix_MCF10A_notreated_batch2$RowNames <- NULL
expression_matrix_MCF10A_TGFbeta1_1day_batch2$RowNames <- NULL
expression_matrix_MCF10A_TGFbeta1_2day_batch2$RowNames <- NULL
expression_matrix_MCF10A_TGFbeta1_3day_batch2$RowNames <- NULL
expression_matrix_MCF10A_TGFbeta1_4day_batch1$RowNames <- NULL
expression_matrix_MCF10A_TGFbeta1_8day_batch1$RowNames <- NULL

dim(expression_matrix_MCF10A_notreated_batch1)
dim(expression_matrix_MCF10A_notreated_batch2)
dim(expression_matrix_MCF10A_TGFbeta1_1day_batch2)
dim(expression_matrix_MCF10A_TGFbeta1_2day_batch2)
dim(expression_matrix_MCF10A_TGFbeta1_3day_batch2)
dim(expression_matrix_MCF10A_TGFbeta1_4day_batch1)
dim(expression_matrix_MCF10A_TGFbeta1_8day_batch1)
dim(expression_matrix_MCF10A)

sum(dim(expression_matrix_MCF10A_notreated_batch1)[2],
    dim(expression_matrix_MCF10A_notreated_batch2)[2],
    dim(expression_matrix_MCF10A_TGFbeta1_1day_batch2)[2],
    dim(expression_matrix_MCF10A_TGFbeta1_2day_batch2)[2],
    dim(expression_matrix_MCF10A_TGFbeta1_3day_batch2)[2],
    dim(expression_matrix_MCF10A_TGFbeta1_4day_batch1)[2],
    dim(expression_matrix_MCF10A_TGFbeta1_8day_batch1)[2])








# Import expression matrix
# exp <- read.csv("expressionMatrixDayZero.csv")
# rownames(exp) <- exp$X
# exp$X <- NULL

# See transcriptogramerStep1 documentation
# First column of the dictionary MUST BE THE ENSEMBL PEPTIDE ID
# SECOND COLUMN MUST BE THE SAME ID AS ON THE GENE EXPRESSION MATRIX ROWNAMES
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
dictionary <- getBM(attributes = c("ensembl_peptide_id", "hgnc_symbol"),
                    mart = ensembl)
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



# Run transcriptogram preprocess, radius = 30
t_notreated_batch1 <- transcriptogramPreprocess(association = assoc, ordering = ord$Protein, radius = 30)
t_notreated_batch1 <- transcriptogramStep1(object = t_notreated_batch1, expression = expression_matrix_MCF10A_notreated_batch1, dictionary = dictionary)
t_notreated_batch1 <- transcriptogramStep2(object = t_notreated_batch1)

t_notreated_batch2 <- transcriptogramPreprocess(association = assoc, ordering = ord$Protein, radius = 30)
t_notreated_batch2 <- transcriptogramStep1(object = t_notreated_batch2, expression = expression_matrix_MCF10A_notreated_batch2, dictionary = dictionary)
t_notreated_batch2 <- transcriptogramStep2(object = t_notreated_batch2)

t_TGFbeta1_1day_batch2 <- transcriptogramPreprocess(association = assoc, ordering = ord$Protein, radius = 30)
t_TGFbeta1_1day_batch2 <- transcriptogramStep1(object = t_TGFbeta1_1day_batch2, expression = expression_matrix_MCF10A_TGFbeta1_1day_batch2, dictionary = dictionary)
t_TGFbeta1_1day_batch2 <- transcriptogramStep2(object = t_TGFbeta1_1day_batch2)

t_TGFbeta1_2day_batch2 <- transcriptogramPreprocess(association = assoc, ordering = ord$Protein, radius = 30)
t_TGFbeta1_2day_batch2 <- transcriptogramStep1(object = t_TGFbeta1_2day_batch2, expression = expression_matrix_MCF10A_TGFbeta1_2day_batch2, dictionary = dictionary)
t_TGFbeta1_2day_batch2 <- transcriptogramStep2(object = t_TGFbeta1_2day_batch2)

t_TGFbeta1_3day_batch2 <- transcriptogramPreprocess(association = assoc, ordering = ord$Protein, radius = 30)
t_TGFbeta1_3day_batch2 <- transcriptogramStep1(object = t_TGFbeta1_3day_batch2, expression = expression_matrix_MCF10A_TGFbeta1_3day_batch2, dictionary = dictionary)
t_TGFbeta1_3day_batch2 <- transcriptogramStep2(object = t_TGFbeta1_3day_batch2)

t_TGFbeta1_4day_batch1 <- transcriptogramPreprocess(association = assoc, ordering = ord$Protein, radius = 30)
t_TGFbeta1_4day_batch1 <- transcriptogramStep1(object = t_TGFbeta1_4day_batch1, expression = expression_matrix_MCF10A_TGFbeta1_4day_batch1, dictionary = dictionary)
t_TGFbeta1_4day_batch1 <- transcriptogramStep2(object = t_TGFbeta1_4day_batch1)

t_TGFbeta1_8day_batch1 <- transcriptogramPreprocess(association = assoc, ordering = ord$Protein, radius = 30)
t_TGFbeta1_8day_batch1 <- transcriptogramStep1(object = t_TGFbeta1_8day_batch1, expression = expression_matrix_MCF10A_TGFbeta1_8day_batch1, dictionary = dictionary)
t_TGFbeta1_8day_batch1 <- transcriptogramStep2(object = t_TGFbeta1_8day_batch1)

# save.image("AposTranscriptogramasDiarios.RData")

save(t_notreated_batch1,
     t_notreated_batch2,
     t_TGFbeta1_1day_batch2,
     t_TGFbeta1_2day_batch2,
     t_TGFbeta1_3day_batch2,
     t_TGFbeta1_4day_batch1,
     t_TGFbeta1_8day_batch1, file = 
       "SomenteTranscriptogramas.RData")
# load("AposTranscriptogramasDiarios.RData")
# load("SomenteTranscriptogramas.RData")

# Criando matrizes de média e desvio padrão.

mean_expression_matrix_MCF10A_notreated_batch1 <- 
  rowMeans(expression_matrix_MCF10A_notreated_batch1[sapply(expression_matrix_MCF10A_notreated_batch1, is.numeric)])
sd_expression_matrix_MCF10A_notreated_batch1 <- 
  apply(expression_matrix_MCF10A_notreated_batch1[sapply(expression_matrix_MCF10A_notreated_batch1, is.numeric)], 1, sd)
mean_sd_expression_matrix_MCF10A_notreated_batch1 <- as.data.frame(cbind(MeanNotreatedBatch1 = mean_expression_matrix_MCF10A_notreated_batch1, 
                                                                         SD_NotreatedBatch1= sd_expression_matrix_MCF10A_notreated_batch1))


mean_expression_matrix_MCF10A_notreated_batch2 <- 
  rowMeans(expression_matrix_MCF10A_notreated_batch2[sapply(expression_matrix_MCF10A_notreated_batch2, is.numeric)])
sd_expression_matrix_MCF10A_notreated_batch2 <- 
  apply(expression_matrix_MCF10A_notreated_batch2[sapply(expression_matrix_MCF10A_notreated_batch2, is.numeric)], 1, sd)
mean_sd_expression_matrix_MCF10A_notreated_batch2 <- as.data.frame(cbind(MeanNotreatedBatch2 = mean_expression_matrix_MCF10A_notreated_batch2, 
                                                                         SD_NotreatedBatch2= sd_expression_matrix_MCF10A_notreated_batch2))

mean_expression_matrix_MCF10A_TGFbeta1_1day_batch2 <- 
  rowMeans(expression_matrix_MCF10A_TGFbeta1_1day_batch2[sapply(expression_matrix_MCF10A_TGFbeta1_1day_batch2, is.numeric)])
sd_expression_matrix_MCF10A_TGFbeta1_1day_batch2 <- 
  apply(expression_matrix_MCF10A_TGFbeta1_1day_batch2[sapply(expression_matrix_MCF10A_TGFbeta1_1day_batch2, is.numeric)], 1, sd)
mean_sd_expression_matrix_MCF10A_TGFbeta1_1day_batch2 <- as.data.frame(cbind(MeanTGFbeta1_1dayBatch2 = mean_expression_matrix_MCF10A_TGFbeta1_1day_batch2, 
                                                                             SD_TGFbeta1_1dayBatch2= sd_expression_matrix_MCF10A_TGFbeta1_1day_batch2))

mean_expression_matrix_MCF10A_TGFbeta1_2day_batch2 <- 
  rowMeans(expression_matrix_MCF10A_TGFbeta1_2day_batch2[sapply(expression_matrix_MCF10A_TGFbeta1_2day_batch2, is.numeric)])
sd_expression_matrix_MCF10A_TGFbeta1_2day_batch2 <- 
  apply(expression_matrix_MCF10A_TGFbeta1_2day_batch2[sapply(expression_matrix_MCF10A_TGFbeta1_2day_batch2, is.numeric)], 1, sd)
mean_sd_expression_matrix_MCF10A_TGFbeta1_2day_batch2 <- as.data.frame(cbind(MeanTGFbeta1_2dayBatch2 = mean_expression_matrix_MCF10A_TGFbeta1_2day_batch2, 
                                                                             SD_TGFbeta1_2dayBatch2= sd_expression_matrix_MCF10A_TGFbeta1_2day_batch2))

mean_expression_matrix_MCF10A_TGFbeta1_3day_batch2 <- 
  rowMeans(expression_matrix_MCF10A_TGFbeta1_3day_batch2[sapply(expression_matrix_MCF10A_TGFbeta1_3day_batch2, is.numeric)])
sd_expression_matrix_MCF10A_TGFbeta1_3day_batch2 <- 
  apply(expression_matrix_MCF10A_TGFbeta1_3day_batch2[sapply(expression_matrix_MCF10A_TGFbeta1_3day_batch2, is.numeric)], 1, sd)
mean_sd_expression_matrix_MCF10A_TGFbeta1_3day_batch2 <- as.data.frame(cbind(MeanTGFbeta1_3dayBatch2 = mean_expression_matrix_MCF10A_TGFbeta1_3day_batch2, 
                                                                             SD_TGFbeta1_3dayBatch2= sd_expression_matrix_MCF10A_TGFbeta1_3day_batch2))

mean_expression_matrix_MCF10A_TGFbeta1_4day_batch1 <- 
  rowMeans(expression_matrix_MCF10A_TGFbeta1_4day_batch1[sapply(expression_matrix_MCF10A_TGFbeta1_4day_batch1, is.numeric)])
sd_expression_matrix_MCF10A_TGFbeta1_4day_batch1 <- 
  apply(expression_matrix_MCF10A_TGFbeta1_4day_batch1[sapply(expression_matrix_MCF10A_TGFbeta1_4day_batch1, is.numeric)], 1, sd)
mean_sd_expression_matrix_MCF10A_TGFbeta1_4day_batch1 <- as.data.frame(cbind(MeanTGFbeta1_4dayBatch1 = mean_expression_matrix_MCF10A_TGFbeta1_4day_batch1, 
                                                                             SD_TGFbeta1_4dayBatch1= sd_expression_matrix_MCF10A_TGFbeta1_4day_batch1))

mean_expression_matrix_MCF10A_TGFbeta1_8day_batch1 <- 
  rowMeans(expression_matrix_MCF10A_TGFbeta1_8day_batch1[sapply(expression_matrix_MCF10A_TGFbeta1_8day_batch1, is.numeric)])
sd_expression_matrix_MCF10A_TGFbeta1_8day_batch1 <- 
  apply(expression_matrix_MCF10A_TGFbeta1_8day_batch1[sapply(expression_matrix_MCF10A_TGFbeta1_8day_batch1, is.numeric)], 1, sd)
mean_sd_expression_matrix_MCF10A_TGFbeta1_8day_batch1 <- as.data.frame(cbind(MeanTGFbeta1_8dayBatch1 = mean_expression_matrix_MCF10A_TGFbeta1_8day_batch1, 
                                                                             SD_TGFbeta1_8dayBatch1= sd_expression_matrix_MCF10A_TGFbeta1_8day_batch1))

# Juntando as matrizes de médias e desvios.

mean_sd_expression_matrix_MCF10A_notreated_batch1$RowNames <- rownames(mean_sd_expression_matrix_MCF10A_notreated_batch1)
mean_sd_expression_matrix_MCF10A_notreated_batch2$RowNames <- rownames(mean_sd_expression_matrix_MCF10A_notreated_batch2)
mean_sd_expression_matrix_MCF10A_TGFbeta1_1day_batch2$RowNames <- rownames(mean_sd_expression_matrix_MCF10A_TGFbeta1_1day_batch2)
mean_sd_expression_matrix_MCF10A_TGFbeta1_2day_batch2$RowNames <- rownames(mean_sd_expression_matrix_MCF10A_TGFbeta1_2day_batch2)
mean_sd_expression_matrix_MCF10A_TGFbeta1_3day_batch2$RowNames <- rownames(mean_sd_expression_matrix_MCF10A_TGFbeta1_3day_batch2)
mean_sd_expression_matrix_MCF10A_TGFbeta1_4day_batch1$RowNames <- rownames(mean_sd_expression_matrix_MCF10A_TGFbeta1_4day_batch1)
mean_sd_expression_matrix_MCF10A_TGFbeta1_8day_batch1$RowNames <- rownames(mean_sd_expression_matrix_MCF10A_TGFbeta1_8day_batch1)

mean_sd_expression_matrix_MCF10A <- full_join(mean_sd_expression_matrix_MCF10A_notreated_batch1, 
                                              mean_sd_expression_matrix_MCF10A_notreated_batch2, by = "RowNames") %>%
  full_join(mean_sd_expression_matrix_MCF10A_TGFbeta1_1day_batch2, by = "RowNames") %>%
  full_join(mean_sd_expression_matrix_MCF10A_TGFbeta1_2day_batch2, by = "RowNames") %>%
  full_join(mean_sd_expression_matrix_MCF10A_TGFbeta1_3day_batch2, by = "RowNames") %>%
  full_join(mean_sd_expression_matrix_MCF10A_TGFbeta1_4day_batch1, by = "RowNames") %>%
  full_join(mean_sd_expression_matrix_MCF10A_TGFbeta1_8day_batch1, by = "RowNames")

# Preenchendo os NAs com 0 na matrizes de médias e desvios.

mean_sd_expression_matrix_MCF10A[is.na(mean_sd_expression_matrix_MCF10A)] <- 0

rownames(mean_sd_expression_matrix_MCF10A) <- mean_sd_expression_matrix_MCF10A$RowNames
mean_sd_expression_matrix_MCF10A$RowNames <- NULL

# expression_matrix_MCF10A using colnames from days
expression_matrix_MCF10A_notreated_batch1_forDay <- expression_matrix_MCF10A_notreated_batch1
names_notreated_batch1 <- paste0("notreated_batch1_", seq_len(ncol(expression_matrix_MCF10A_notreated_batch1_forDay)))
colnames(expression_matrix_MCF10A_notreated_batch1_forDay) <- names_notreated_batch1

expression_matrix_MCF10A_notreated_batch2_forDay <- expression_matrix_MCF10A_notreated_batch2
names_notreated_batch2 <- paste0("notreated_batch2_", seq_len(ncol(expression_matrix_MCF10A_notreated_batch2_forDay)))
colnames(expression_matrix_MCF10A_notreated_batch2_forDay) <- names_notreated_batch2

expression_matrix_MCF10A_TGFbeta1_1day_batch2_forDay <- expression_matrix_MCF10A_TGFbeta1_1day_batch2
names_TGFbeta1_1day_batch2 <- paste0("TGFbeta1_1day_batch2_", seq_len(ncol(expression_matrix_MCF10A_TGFbeta1_1day_batch2_forDay)))
colnames(expression_matrix_MCF10A_TGFbeta1_1day_batch2_forDay) <- names_TGFbeta1_1day_batch2

expression_matrix_MCF10A_TGFbeta1_2day_batch2_forDay <- expression_matrix_MCF10A_TGFbeta1_2day_batch2
names_TGFbeta1_2day_batch2 <- paste0("TGFbeta1_2day_batch2_", seq_len(ncol(expression_matrix_MCF10A_TGFbeta1_2day_batch2_forDay)))
colnames(expression_matrix_MCF10A_TGFbeta1_2day_batch2_forDay) <- names_TGFbeta1_2day_batch2

expression_matrix_MCF10A_TGFbeta1_3day_batch2_forDay <- expression_matrix_MCF10A_TGFbeta1_3day_batch2
names_TGFbeta1_3day_batch2 <- paste0("TGFbeta1_3day_batch2_", seq_len(ncol(expression_matrix_MCF10A_TGFbeta1_3day_batch2_forDay)))
colnames(expression_matrix_MCF10A_TGFbeta1_3day_batch2_forDay) <- names_TGFbeta1_3day_batch2

expression_matrix_MCF10A_TGFbeta1_4day_batch1_forDay <- expression_matrix_MCF10A_TGFbeta1_4day_batch1
names_TGFbeta1_4day_batch1 <- paste0("TGFbeta1_4day_batch1_", seq_len(ncol(expression_matrix_MCF10A_TGFbeta1_4day_batch1_forDay)))
colnames(expression_matrix_MCF10A_TGFbeta1_4day_batch1_forDay) <- names_TGFbeta1_4day_batch1

expression_matrix_MCF10A_TGFbeta1_8day_batch1_forDay <- expression_matrix_MCF10A_TGFbeta1_8day_batch1
names_TGFbeta1_8day_batch1 <- paste0("TGFbeta1_8day_batch1_", seq_len(ncol(expression_matrix_MCF10A_TGFbeta1_8day_batch1_forDay)))
colnames(expression_matrix_MCF10A_TGFbeta1_8day_batch1_forDay) <- names_TGFbeta1_8day_batch1

expression_matrix_MCF10A_notreated_batch1_forDay$RowNames <- rownames(expression_matrix_MCF10A_notreated_batch1_forDay)
expression_matrix_MCF10A_notreated_batch2_forDay$RowNames <- rownames(expression_matrix_MCF10A_notreated_batch2_forDay)
expression_matrix_MCF10A_TGFbeta1_1day_batch2_forDay$RowNames <- rownames(expression_matrix_MCF10A_TGFbeta1_1day_batch2_forDay)
expression_matrix_MCF10A_TGFbeta1_2day_batch2_forDay$RowNames <- rownames(expression_matrix_MCF10A_TGFbeta1_2day_batch2_forDay)
expression_matrix_MCF10A_TGFbeta1_3day_batch2_forDay$RowNames <- rownames(expression_matrix_MCF10A_TGFbeta1_3day_batch2_forDay)
expression_matrix_MCF10A_TGFbeta1_4day_batch1_forDay$RowNames <- rownames(expression_matrix_MCF10A_TGFbeta1_4day_batch1_forDay)
expression_matrix_MCF10A_TGFbeta1_8day_batch1_forDay$RowNames <- rownames(expression_matrix_MCF10A_TGFbeta1_8day_batch1_forDay)

expression_matrix_MCF10A_forDay <- 
  full_join(expression_matrix_MCF10A_notreated_batch1_forDay, expression_matrix_MCF10A_notreated_batch2_forDay, by = "RowNames") %>%
  full_join(expression_matrix_MCF10A_TGFbeta1_1day_batch2_forDay, by = "RowNames") %>%
  full_join(expression_matrix_MCF10A_TGFbeta1_2day_batch2_forDay, by = "RowNames") %>%
  full_join(expression_matrix_MCF10A_TGFbeta1_3day_batch2_forDay, by = "RowNames") %>%
  full_join(expression_matrix_MCF10A_TGFbeta1_4day_batch1_forDay, by = "RowNames") %>%
  full_join(expression_matrix_MCF10A_TGFbeta1_8day_batch1_forDay, by = "RowNames")


# Preenchendo os NAs com 0
expression_matrix_MCF10A_forDay[is.na(expression_matrix_MCF10A_forDay)] <- 0


# Movendo os nomes das linhas de volta para rownames e removendo a coluna auxiliar
rownames(expression_matrix_MCF10A_forDay) <- expression_matrix_MCF10A_forDay$RowNames
expression_matrix_MCF10A_forDay$RowNames <- NULL

expression_matrix_MCF10A_notreated_batch1_forDay$RowNames <- NULL
expression_matrix_MCF10A_notreated_batch2_forDay$RowNames <- NULL
expression_matrix_MCF10A_TGFbeta1_1day_batch2_forDay$RowNames <- NULL
expression_matrix_MCF10A_TGFbeta1_2day_batch2_forDay$RowNames <- NULL
expression_matrix_MCF10A_TGFbeta1_3day_batch2_forDay$RowNames <- NULL
expression_matrix_MCF10A_TGFbeta1_4day_batch1_forDay$RowNames <- NULL
expression_matrix_MCF10A_TGFbeta1_8day_batch1_forDay$RowNames <- NULL

View(expression_matrix_MCF10A_forDay)
write.table(expression_matrix_MCF10A_forDay, 
            file = "expression_matrix_MCF10A_forDay_ENSG.txt", 
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
save.image("AposMatrizFullForDay.RData")



mean_t_notreated_batch1 <- rowMeans(t_notreated_batch1@transcriptogramS2[, -c(1,2)])
sd_t_notreated_batch1 <- apply(t_notreated_batch1@transcriptogramS2[, -c(1,2)], 1, sd)
mean_sd_t_notreated_batch1 <- cbind(mean_t_notreated_batch1, sd_t_notreated_batch1)
full_mean_sd_t_notreated_batch1 <- cbind(t_notreated_batch1@transcriptogramS2, mean_sd_t_notreated_batch1)
full_mean_sd_t_notreated_batch1 <- full_mean_sd_t_notreated_batch1[, c(1, 2, ncol(full_mean_sd_t_notreated_batch1)-1, ncol(full_mean_sd_t_notreated_batch1))]
write.table(full_mean_sd_t_notreated_batch1, file = "full_mean_sd_t_notreated_batch1.txt", row.names = F, col.names = T, sep = "\t")

mean_t_notreated_batch2 <- rowMeans(t_notreated_batch2@transcriptogramS2[, -c(1,2)])
sd_t_notreated_batch2 <- apply(t_notreated_batch2@transcriptogramS2[, -c(1,2)], 1, sd)
mean_sd_t_notreated_batch2 <- cbind(mean_t_notreated_batch2, sd_t_notreated_batch2)
full_mean_sd_t_notreated_batch2 <- cbind(t_notreated_batch2@transcriptogramS2, mean_sd_t_notreated_batch2)
full_mean_sd_t_notreated_batch2 <- full_mean_sd_t_notreated_batch2[, c(1, 2, ncol(full_mean_sd_t_notreated_batch2)-1, ncol(full_mean_sd_t_notreated_batch2))]
write.table(full_mean_sd_t_notreated_batch2, file = "full_mean_sd_t_notreated_batch2.txt", row.names = F, col.names = T, sep = "\t")

mean_t_TGFbeta1_1day_batch2 <- rowMeans(t_TGFbeta1_1day_batch2@transcriptogramS2[, -c(1,2)])
sd_t_TGFbeta1_1day_batch2 <- apply(t_TGFbeta1_1day_batch2@transcriptogramS2[, -c(1,2)], 1, sd)
mean_sd_t_TGFbeta1_1day_batch2 <- cbind(mean_t_TGFbeta1_1day_batch2, sd_t_TGFbeta1_1day_batch2)
full_mean_sd_t_TGFbeta1_1day_batch2 <- cbind(t_TGFbeta1_1day_batch2@transcriptogramS2, mean_sd_t_TGFbeta1_1day_batch2)
full_mean_sd_t_TGFbeta1_1day_batch2 <- full_mean_sd_t_TGFbeta1_1day_batch2[, c(1, 2, ncol(full_mean_sd_t_TGFbeta1_1day_batch2)-1, ncol(full_mean_sd_t_TGFbeta1_1day_batch2))]
write.table(full_mean_sd_t_TGFbeta1_1day_batch2, file = "full_mean_sd_t_TGFbeta1_1day_batch2.txt", row.names = F, col.names = T, sep = "\t")

mean_t_TGFbeta1_2day_batch2 <- rowMeans(t_TGFbeta1_2day_batch2@transcriptogramS2[, -c(1,2)])
sd_t_TGFbeta1_2day_batch2 <- apply(t_TGFbeta1_2day_batch2@transcriptogramS2[, -c(1,2)], 1, sd)
mean_sd_t_TGFbeta1_2day_batch2 <- cbind(mean_t_TGFbeta1_2day_batch2, sd_t_TGFbeta1_2day_batch2)
full_mean_sd_t_TGFbeta1_2day_batch2 <- cbind(t_TGFbeta1_2day_batch2@transcriptogramS2, mean_sd_t_TGFbeta1_2day_batch2)
full_mean_sd_t_TGFbeta1_2day_batch2 <- full_mean_sd_t_TGFbeta1_2day_batch2[, c(1, 2, ncol(full_mean_sd_t_TGFbeta1_2day_batch2)-1, ncol(full_mean_sd_t_TGFbeta1_2day_batch2))]
write.table(full_mean_sd_t_TGFbeta1_2day_batch2, file = "full_mean_sd_t_TGFbeta1_2day_batch2.txt", row.names = F, col.names = T, sep = "\t")

mean_t_TGFbeta1_3day_batch2 <- rowMeans(t_TGFbeta1_3day_batch2@transcriptogramS2[, -c(1,2)])
sd_t_TGFbeta1_3day_batch2 <- apply(t_TGFbeta1_3day_batch2@transcriptogramS2[, -c(1,2)], 1, sd)
mean_sd_t_TGFbeta1_3day_batch2 <- cbind(mean_t_TGFbeta1_3day_batch2, sd_t_TGFbeta1_3day_batch2)
full_mean_sd_t_TGFbeta1_3day_batch2 <- cbind(t_TGFbeta1_3day_batch2@transcriptogramS2, mean_sd_t_TGFbeta1_3day_batch2)
full_mean_sd_t_TGFbeta1_3day_batch2 <- full_mean_sd_t_TGFbeta1_3day_batch2[, c(1, 2, ncol(full_mean_sd_t_TGFbeta1_3day_batch2)-1, ncol(full_mean_sd_t_TGFbeta1_3day_batch2))]
write.table(full_mean_sd_t_TGFbeta1_3day_batch2, file = "full_mean_sd_t_TGFbeta1_3day_batch2.txt", row.names = F, col.names = T, sep = "\t")

mean_t_TGFbeta1_4day_batch1 <- rowMeans(t_TGFbeta1_4day_batch1@transcriptogramS2[, -c(1,2)])
sd_t_TGFbeta1_4day_batch1 <- apply(t_TGFbeta1_4day_batch1@transcriptogramS2[, -c(1,2)], 1, sd)
mean_sd_t_TGFbeta1_4day_batch1 <- cbind(mean_t_TGFbeta1_4day_batch1, sd_t_TGFbeta1_4day_batch1)
full_mean_sd_t_TGFbeta1_4day_batch1 <- cbind(t_TGFbeta1_4day_batch1@transcriptogramS2, mean_sd_t_TGFbeta1_4day_batch1)
full_mean_sd_t_TGFbeta1_4day_batch1 <- full_mean_sd_t_TGFbeta1_4day_batch1[, c(1, 2, ncol(full_mean_sd_t_TGFbeta1_4day_batch1)-1, ncol(full_mean_sd_t_TGFbeta1_4day_batch1))]
write.table(full_mean_sd_t_TGFbeta1_4day_batch1, file = "full_mean_sd_t_TGFbeta1_4day_batch1.txt", row.names = F, col.names = T, sep = "\t")

mean_t_TGFbeta1_8day_batch1 <- rowMeans(t_TGFbeta1_8day_batch1@transcriptogramS2[, -c(1,2)])
sd_t_TGFbeta1_8day_batch1 <- apply(t_TGFbeta1_8day_batch1@transcriptogramS2[, -c(1,2)], 1, sd)
mean_sd_t_TGFbeta1_8day_batch1 <- cbind(mean_t_TGFbeta1_8day_batch1, sd_t_TGFbeta1_8day_batch1)
full_mean_sd_t_TGFbeta1_8day_batch1 <- cbind(t_TGFbeta1_8day_batch1@transcriptogramS2, mean_sd_t_TGFbeta1_8day_batch1)
full_mean_sd_t_TGFbeta1_8day_batch1 <- full_mean_sd_t_TGFbeta1_8day_batch1[, c(1, 2, ncol(full_mean_sd_t_TGFbeta1_8day_batch1)-1, ncol(full_mean_sd_t_TGFbeta1_8day_batch1))]
write.table(full_mean_sd_t_TGFbeta1_8day_batch1, file = "full_mean_sd_t_TGFbeta1_8day_batch1.txt", row.names = F, col.names = T, sep = "\t")

# Visualization Of The Transcriptogramas

load("SomenteTranscriptogramas.RData")

map_dfr(ls()[grep("^t", ls())], function(x) {
  
  #normalizar
  transcriptogram <- apply(get(x, envir = .GlobalEnv)@transcriptogramS2[, -c(1,2)], 2, function(x) {x / sum(x)})
  
  data.frame(
    Protein = get(x, envir = .GlobalEnv)@transcriptogramS2[,1],
    Position = get(x, envir = .GlobalEnv)@transcriptogramS2[,2],
    mean = apply(transcriptogram, 1, mean, na.rm = T),
    sd = apply(transcriptogram, 1, sd,  na.rm = T),
    treatment = x
  ) -> res
  return(res)
  
}) -> df


png(file = "images/MeanVariationInFunctionInThePositionProtein.png", width = 1920, height = 1080, res = 150)
ggplot(df, aes(x = Position, y = mean, col = treatment, group = treatment)) +
  geom_line() +
  facet_wrap(~treatment)
dev.off()

# Divide o dataframe por batchs.

df_batch1 <- subset(df,
                    treatment =="t_notreated_batch1" | 
                      treatment =="t_TGFbeta1_4day_batch1" | 
                      treatment =="t_TGFbeta1_8day_batch1")

df_batch2 <- subset(df,
                    treatment =="t_notreated_batch2" | 
                      treatment =="t_TGFbeta1_1day_batch2" |
                      treatment =="t_TGFbeta1_2day_batch2" |
                      treatment =="t_TGFbeta1_3day_batch2")



# Divide CASO por CONTROLE

df_notreated_batch1 <- subset(df,treatment =="t_notreated_batch1")

df_notreated_batch2 <- subset(df,treatment =="t_notreated_batch2")

df_TGFbeta1_1day_batch2 <- subset(df,treatment =="t_TGFbeta1_1day_batch2")

df_TGFbeta1_2day_batch2 <- subset(df,treatment =="t_TGFbeta1_2day_batch2")

df_TGFbeta1_3day_batch2 <- subset(df,treatment =="t_TGFbeta1_3day_batch2")

df_TGFbeta1_4day_batch1 <- subset(df,treatment =="t_TGFbeta1_4day_batch1")

df_TGFbeta1_8day_batch1 <- subset(df,treatment =="t_TGFbeta1_8day_batch1")

# Divide t_TGFbeta1_4day_batch1 por t_nontreated_batch1

df$RatioMean_4day_per_notreated_batch1 <- df_TGFbeta1_4day_batch1$mean / df_notreated_batch1$mean

# Divide t_TGFbeta1_8day_batch1 por t_nontreated_batch1

df$RatioMean_8day_per_notreated_batch1 <- df_TGFbeta1_8day_batch1$mean / df_notreated_batch1$mean

# Divide t_TGFbeta1_1day_batch2 por t_nontreated_batch2

df$RatioMean_1day_per_notreated_batch2 <- df_TGFbeta1_1day_batch2$mean / df_notreated_batch2$mean

# Divide t_TGFbeta1_2day_batch2 por t_nontreated_batch2

df$RatioMean_2day_per_notreated_batch2 <- df_TGFbeta1_2day_batch2$mean / df_notreated_batch2$mean

# Divide t_TGFbeta1_3day_batch2 por t_nontreated_batch2

df$RatioMean_3day_per_notreated_batch2 <- df_TGFbeta1_3day_batch2$mean / df_notreated_batch2$mean



df_long <- df %>%
  pivot_longer(cols = c(RatioMean_4day_per_notreated_batch1, 
                        RatioMean_8day_per_notreated_batch1,
                        RatioMean_1day_per_notreated_batch2,
                        RatioMean_2day_per_notreated_batch2,
                        RatioMean_3day_per_notreated_batch2),
               names_to = "Ratio_Type",
               values_to = "Ratio_Value")

png(file = "images/MeanRatioPerPositionProtein.png", width = 1920, height = 1080, res = 200)
ggplot(df_long, aes(x = Position, y = Ratio_Value, color = Ratio_Type, group = Ratio_Type)) +
  geom_line() +
  labs(title = "Comparação de Ratios ao longo dos dias",
       x = "Position",
       y = "Ratio") +
  theme_minimal()
dev.off()

## COMPARANDO AS MÉDIAS DOS TRANSCRIPTOGRAMAS  NOTREATED DE AMBOS OS BATCHS.

# Adicionar uma coluna identificando cada batch
df_notreated_batch1 <- df_notreated_batch1 %>% mutate(batch = "notreated_batch1")
df_notreated_batch2 <- df_notreated_batch2 %>% mutate(batch = "notreated_batch2")

# Combinar os dois dataframes em um só
df_combined <- bind_rows(df_notreated_batch1, df_notreated_batch2)

# Criar o gráfico com as duas linhas
ggplot(df_combined, aes(x = Position, y = mean, color = batch, group = batch)) +
  geom_line() +
  labs(title = "Comparation of Notreated Batches", x = "Position", y = "Mean") +
  theme_minimal() -> PlotComparationOfNotreatedBatches

# Plotando em qualidade de artigo.

ggsave(filename = "images/ComparationOfNotreatedBatchesMean_GGSAVE.png",
       plot = PlotComparationOfNotreatedBatches, 
       width = 10, height = 8, dpi = 600, units = "in")


# Salvando o gráfico em PNG com qualidade de artigo científico
ggplot(df, aes(x = Position, y = mean, col = treatment, group = treatment)) +
  geom_line() +
  facet_wrap(~treatment) -> p

ggsave(filename = "images/MeanVariationInFunctionInThePositionProtein_GGSAVE.png", plot = p, 
       width = 19, height = 10, dpi = 1200, units = "in")


# # Somente os transcriptogramas.
# 
# png(file = "images/t_notreated_batch1_HeatMap.png", width = 1920, height = 1080, res = 150)
# Heatmap(t_notreated_batch1@transcriptogramS2[, -c(1,2)])
# dev.off()
# 
# 
# 
# png(file = "images/t_notreated_batch2_HeatMap.png", width = 1920, height = 1080, res = 150)
# Heatmap(t_notreated_batch2@transcriptogramS2, parallel = 40)
# dev.off()
# 
# 
# 
# png(file = "images/t_TGFbeta1_1day_batch2_HeatMap.png", width = 1920, height = 1080, res = 150)
# Heatmap(t_TGFbeta1_1day_batch2@transcriptogramS2, parallel = 40)
# dev.off()
# 
# 
# png(file = "images/t_TGFbeta1_2day_batch2_HeatMap.png", width = 1920, height = 1080, res = 150)
# Heatmap(t_TGFbeta1_2day_batch2@transcriptogramS2, parallel = 40)
# dev.off()
# 
# 
# 
# png(file = "images/t_TGFbeta1_3day_batch2_HeatMap.png", width = 1920, height = 1080, res = 150)
# Heatmap(t_TGFbeta1_3day_batch2@transcriptogramS2, parallel = 40)
# dev.off()
# 
# 
# 
# png(file = "images/t_TGFbeta1_4day_batch1_HeatMap.png", width = 1920, height = 1080, res = 150)
# Heatmap(t_TGFbeta1_4day_batch1@transcriptogramS2, parallel = 40)
# dev.off()
# 
# 
# 
# png(file = "images/t_TGFbeta1_8day_batch1_HeatMap.png", width = 1920, height = 1080, res = 150)
# Heatmap(t_TGFbeta1_8day_batch1@transcriptogramS2, parallel = 40)
# dev.off()
# 
# save.image("AposHeatMapTranscriptogramasDiarios.RData")



# t_notreated_batch1 <- differentiallyExpressed(object = t_notreated_batch1, levels = levels, pValue = 0.01,
#                                               trend = FALSE, title = "radius 80")
# 
# t_notreated_batch2 <- differentiallyExpressed(object = t_notreated_batch2, pValue = 0.01,
#                                               trend = FALSE, title = "radius 80")
# 
# t_TGFbeta1_1day_batch2 <- differentiallyExpressed(object = t_TGFbeta1_1day_batch2, pValue = 0.01,
#                                               trend = FALSE, title = "radius 80")
# 
# t_TGFbeta1_2day_batch2 <- differentiallyExpressed(object = t_TGFbeta1_2day_batch2, pValue = 0.01,
#                                               trend = FALSE, title = "radius 80")
# 
# t_TGFbeta1_3day_batch2 <- differentiallyExpressed(object = t_TGFbeta1_3day_batch2, pValue = 0.01,
#                                               trend = FALSE, title = "radius 80")
# 
# t_TGFbeta1_4day_batch1 <- differentiallyExpressed(object = t_TGFbeta1_4day_batch1, pValue = 0.01,
#                                               trend = FALSE, title = "radius 80")
# 
# t_TGFbeta1_8day_batch1 <- differentiallyExpressed(object = t_TGFbeta1_8day_batch1, pValue = 0.01,
#                                               trend = FALSE, title = "radius 80")

# Transcriptogramer de todos os dias.

# t <- transcriptogramPreprocess(association = assoc, ordering = ord$Protein, radius = 80)
# t <- transcriptogramStep1(object = t, expression = expression_matrix_MCF10A, dictionary = dictionary, nCores = 10)
# t <- transcriptogramStep2(object = t, nCores = 10)

t_mean_sd <- transcriptogramPreprocess(association = assoc, ordering = ord$Protein, radius = 80)
t_mean_sd <- transcriptogramStep1(object = t_mean_sd, expression = mean_sd_expression_matrix_MCF10A, dictionary = dictionary, nCores = 10)
t_mean_sd <- transcriptogramStep2(object = t_mean_sd, nCores = 10)


# Transcriptogramer com todos os dias juntos.
# save.image("AposTranscriptogramasDiarios.RData")
# load("AposTranscriptogramasDiarios.RData")


# subset <- expression_matrix[, c(1:10)]
# write.csv(subset, file= "~/mestrado-single-cell/expressionMatrixDayZero.csv", row.names = T)
# Subamostragem ou agregação pode ser necessária devido ao grande volume de dados
# Aqui está um exemplo simples de subamostragem (selecionando um subconjunto de células)
# set.seed(123)
# subsampled_matrix <- expression_matrix[, sample(ncol(expression_matrix), 500)]

# 3. Análise com transcriptogramer

# Carregar dados de genes.
# genes <- read.table("mtx_conversions/notreated_batch1/outs/filtered_feature_bc_matrix/genes.tsv", header = FALSE, sep = "\t")
# barcodes <- read.table("mtx_conversions/notreated_batch1/outs/filtered_feature_bc_matrix/barcodes.tsv", header = FALSE, sep = "\t")

# Conectando ao Ensembl BioMart
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Consultando informações de genes
# dictionary <- getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id"),
#                     filters = "ensembl_gene_id", values = genes, mart = ensembl)


# dictionary <- dictionary[!apply(dictionary, 1, function(row) any(is.na(row) | row == "")), ]

# colnames(dictionary) <- NULL
# 
# head(dictionary)
# 
# class(dictionary)
# 
# dim(dictionary)
# 
# proteinas <- dictionary[,2]
# 
# head(proteinas)
# 
# write.csv(proteinas, file = "proteinasLista.csv", row.names = F, quote = F)
# 
# head(expression_matrix)
# 
# class(expression_matrix)
# 
# head(expression_matrix[, 'AAACCTGCATTGGTAC-1'], digits = 5, 20)
# 
# AAACCTGCATTGGTAC <- expression_matrix[, 'AAACCTGCATTGGTAC-1', drop = FALSE]
# 
# AAACCTGCATTGGTAC <- as.data.frame(AAACCTGCATTGGTAC)
# 
# head(AAACCTGCATTGGTAC, digits = 5, 20)
# 
# class(AAACCTGCATTGGTAC)
# 
# colnames(dictionary) <- c("V1","V2")
# 
# head(genes)
# 
# head(dictionary)
# 
# merged <- merge(dictionary, genes, by = "V1",all = TRUE)
# 
# head(merged)
# 
# dictionaryP <- merged[,2:3]
# 
# colnames(dictionaryP) <- c("ensembl_peptide_id","hgnc_id")
# 
# head(dictionaryP)
# 
# dim(dictionaryP)
# 
# head(AAACCTGCATTGGTAC, nsmall = 5, 20)
# 
# dim(AAACCTGCATTGGTAC)
# 
# head(rownames(AAACCTGCATTGGTAC))
# 
# sum(rownames(AAACCTGCATTGGTAC) %in% dictionaryP[,2])
# 
# # Criar um data frame de exemplo
# 
# # Exibir o data frame original
# 
# 
# # Adicionar o termo "9606." à coluna 'ensembl_peptide_id'
# dictionaryP$ensembl_peptide_id <- paste("9606.", dictionaryP$ensembl_peptide_id, sep = "")
# 
# head(dictionaryP)
# 
# 
# ### TRANSCRIPTOGRAMER
# 
# ## creating the object and setting the radius as 80
# t <- transcriptogramPreprocess(association = association, ordering = Hs900,
#                                radius = 80)
## after the preprocessing

## modifying the radius of an existing Transcriptogram object
# radius(object = t) <- 50
# 
# ## getting the radius of an existing Transcriptogram object
# r <- radius(object = t)
# 
# oPropertiesR50 <- orderingProperties(object = t, nCores = 5)
# 
# ## slight change of radius
# radius(object = t) <- 80

## this output is partially different comparing to oPropertiesR50
# oPropertiesR80 <- orderingProperties(object = t, nCores = 10)
# 
# # sum(oPropertiesR50$windowModularity)
# 
# sum(oPropertiesR80$windowModularity)
# 
# cProperties <- connectivityProperties(object = t)
# 
# t <- transcriptogramStep1(object = t, expression = AAACCTGCATTGGTAC,
#                           dictionary = dictionaryP, nCores = 10)
# # t2 <- t
# 
# head(t@transcriptogramS1)
# 
# t <- transcriptogramStep2(object = t)

# # radius(object = t2) <- 50
# 
# t2 <- transcriptogramStep2(object = t2, nCores = 5)

## trend = FALSE for microarray data or voom log2-counts-per-million
## the default value for trend is FALSE

# controle <- dim(expression_matrix_MCF10A_notreated_batch1)[2] + dim(expression_matrix_MCF10A_notreated_batch2)[2]
# caso <- dim(expression_matrix_MCF10A_TGFbeta1_1day_batch2)[2] + 
#   dim(expression_matrix_MCF10A_TGFbeta1_2day_batch2)[2] +
#   dim(expression_matrix_MCF10A_TGFbeta1_3day_batch2)[2] +
#   dim(expression_matrix_MCF10A_TGFbeta1_4day_batch1)[2] +
#   dim(expression_matrix_MCF10A_TGFbeta1_8day_batch1)[2]







#### A PARTIR DAQUI, ESTÁ DANDO ERRO, AFIRMNADO QUE NÃO HÁ DIFERENÇA DE EXPRESSÃO.
levels <- c(rep(FALSE, 4),
            rep(TRUE, 10))


t_mean_sd <- differentiallyExpressed(object = t_mean_sd, levels = levels, pValue = 0.5,
                                     trend = FALSE, title = "radius 80")

# ## the radius 50 will affect the output significantly
# t2 <- differentiallyExpressed(object = t2, levels = levels, pValue = 0.01,
#                               species = DEsymbols, title = "radius 50")

## using the species argument to translate ENSEMBL Peptide IDs to Symbols
## Internet connection is required for this command
t_mean_sd <- differentiallyExpressed(object = t_mean_sd, levels = levels, pValue = 0.01,
                                     species = "Homo sapiens", title = "radius 80")

## translating ENSEMBL Peptide IDs to Symbols using the DEsymbols dataset
t_mean_sd <- differentiallyExpressed(object = t_mean_sd, levels = levels, pValue = 0.01,
                                     species = DEsymbols, title = "radius 80")

DE <- DE(object = t_mean_sd)
# DE2 <- DE(object = t2)
nrow(DE)

nrow(DE)/length(unique(DE$ClusterNumber))

# nrow(DE2)

# nrow(DE2)/length(unique(DE2$ClusterNumber))

rdp <- clusterVisualization(object = t_mean_sd)

## using the HsBPTerms dataset to create the Protein2GO data.frame
t_mean_sd <- clusterEnrichment(object = t_mean_sd, species = HsBPTerms,
                               pValue = 0.005, nCores = 5)


## using the species argument to create the Protein2GO data.frame
## Internet connection is required for this command
t_mean_sd <- clusterEnrichment(object = t_mean_sd, species = "Homo sapiens",
                               pValue = 0.005, nCores = 5)

