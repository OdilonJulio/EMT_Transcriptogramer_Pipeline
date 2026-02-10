library(Seurat)
library(SeuratObject)
library(dplyr)
library(Matrix)
library(SingleCellExperiment)
library(scDblFinder)
library(tidyverse)
library(biomaRt)
library(transcriptogramer)
library(rsparse)
library(Ropj)
library(vroom)
library(ggplot2)
library(patchwork)
library(parallel)
library(ComplexHeatmap)
library(magrittr)
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
  "notreated_batch1", "notreated_batch2", "TGFbeta1_1day_batch2",
  "TGFbeta1_2day_batch2", "TGFbeta1_3day_batch2", "TGFbeta1_4day_batch1", "TGFbeta1_8day_batch1"
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

# Visualizações antes do controle de qualidade
VlnPlot(
  combined_seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0,
  group.by = "orig.ident"
)

# Aplicar controle de qualidade para cada objeto Seurat individual
filtered_seurat_objects <- lapply(seurat_objects, function(seurat) {
  subset(seurat, subset = nFeature_RNA > 500 & percent.mt < 20)
})

# Combinar os objetos filtrados
filtered_combined_seurat <- merge(
  filtered_seurat_objects[[1]],
  y = filtered_seurat_objects[-1],
  add.cell.ids = sample_names,
  project = "FilteredCombined"
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



### Normalização, integração e extração da matriz com colnames ajustados ###

# Normalizar os dados usando RC (Relative Counts)
normalized_seurat_objects <- lapply(filtered_seurat_objects, function(seurat) {
  NormalizeData(seurat, normalization.method = "RC")
})

# Encontrar genes compartilhados entre os objetos para integração
shared_features <- Reduce(intersect, lapply(normalized_seurat_objects, rownames))

# Encontrar os anchors para integração
integration_anchors <- FindIntegrationAnchors(
  object.list = normalized_seurat_objects,
  anchor.features = shared_features
)

# Integrar os dados
integrated_seurat <- IntegrateData(anchorset = integration_anchors)

# Atualizar os colnames para o formato sample_name_barcode
colnames(integrated_seurat) <- make.unique(
  paste0(
    integrated_seurat$orig.ident, "_", Cells(integrated_seurat)
  )
)

# Extrair a matriz de expressão normalizada
integrated_expression_matrix <- as.matrix(integrated_seurat@assays$integrated@data)

# Visualizar a matriz de expressão
head(integrated_expression_matrix)


save(integrated_expression_matrix, file = "integrated_expression_matrix.RData")

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
  values = rownames(integrated_seurat),  # Usar os nomes dos genes
  mart = ensembl
)

# Limpar o mapeamento removendo valores vazios e NAs
gene_mapping <- gene_mapping %>%
  mutate(ensembl_gene_id = ifelse(ensembl_gene_id == "", NA, ensembl_gene_id)) %>%
  na.omit()  # Remover linhas com NA em 'ensembl_gene_id'

# Alinhar os nomes das linhas com o mapeamento
mapped_genes <- gene_mapping$ensembl_gene_id[match(
  rownames(integrated_expression_matrix),
  gene_mapping$external_gene_name
)]

# Filtrar para manter apenas as linhas com mapeamentos válidos
valid_indices <- !is.na(mapped_genes)
integrated_expression_matrix <- integrated_expression_matrix[valid_indices, , drop = FALSE]  # Filtrar linhas válidas
rownames(integrated_expression_matrix) <- mapped_genes[valid_indices]  # Substituir os nomes das linhas

# Verificar os primeiros nomes trocados
head(rownames(integrated_expression_matrix))  # Verifique os novos nomes de genes (ENSEMBL Gene IDs - ENSG)

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
t_integrated_seurat_R0 <- transcriptogramPreprocess(association = assoc, 
                                                    ordering = ord$Protein, 
                                                    radius = 0)

# Rodar o Transcriptogramer
t_integrated_seurat_R0 <- transcriptogramStep1(object = t_integrated_seurat_R0, 
                                               expression = integrated_expression_matrix, 
                                               dictionary = dictionary)
t_integrated_seurat_R0 <- transcriptogramStep2(object = t_integrated_seurat_R0)

# Salvar os objetos processados
save(t_integrated_seurat_R0, file = "t_integrated_seurat_R0.RData")

# Pré-processamento com raio = 30
t_integrated_seurat_R30 <- transcriptogramPreprocess(association = assoc, 
                                                    ordering = ord$Protein, 
                                                    radius = 30)

# Rodar o Transcriptogramer
t_integrated_seurat_R30 <- transcriptogramStep1(object = t_integrated_seurat_R30, 
                                               expression = integrated_expression_matrix, 
                                               dictionary = dictionary)
t_integrated_seurat_R30 <- transcriptogramStep2(object = t_integrated_seurat_R30)

# Salvar os objetos processados
save(t_integrated_seurat_R30, file = "t_integrated_seurat_R30.RData")

##################################################################
##### SUBSETs com 20 células de cada dia.

df_subset_t_integrated_seurat_R0 <- t_integrated_seurat_R0@transcriptogramS2
df_subset_t_integrated_seurat_R30 <- t_integrated_seurat_R30@transcriptogramS2

# Lista de condições para filtrar
conditions <- c("notreated_batch1", "notreated_batch2", "TGFbeta1_1day_batch2", 
                "TGFbeta1_2day_batch2", "TGFbeta1_3day_batch2", 
                "TGFbeta1_4day_batch1", "TGFbeta1_8day_batch1")

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
df_subset_t_integrated_seurat_R0_filtered <- filter_barcodes(df_subset_t_integrated_seurat_R0, conditions)
df_subset_t_integrated_seurat_R30_filtered <- filter_barcodes(df_subset_t_integrated_seurat_R30, conditions)

# Exportar os resultados para CSV
write.csv(df_subset_t_integrated_seurat_R0_filtered, "df_subset_t_integrated_seurat_R0_filtered.csv", row.names = TRUE)
write.csv(df_subset_t_integrated_seurat_R30_filtered, "df_subset_t_integrated_seurat_R30_filtered.csv", row.names = TRUE)



################# PCA dos Trascriptogramas

# Extraindo o dataframe do slot transcriptogramS2

library(ggplot2)
library(gridExtra)


df_t_integrated_seurat_R0 <- t_integrated_seurat_R0@transcriptogramS2[,-2]

df_t_integrated_seurat_R30 <- t_integrated_seurat_R30@transcriptogramS2[,-2]
# Função para calcular PCA e variância explicada
perform_pca <- function(data) {
  pca <- prcomp(data, scale. = TRUE)  # Escalar os dados para evitar viés
  var_explained <- (pca$sdev^2) / sum(pca$sdev^2)  # Variância explicada por cada PC
  cumulative_variance <- cumsum(var_explained)     # Variância acumulada
  list(pca = pca, var_explained = var_explained, cumulative_variance = cumulative_variance)
}

# Calcular PCA para os dois datasets
pca_r0 <- perform_pca(df_t_integrated_seurat_R0[,-1])  # Ignorar a coluna de ENSP
pca_r30 <- perform_pca(df_t_integrated_seurat_R30[,-1])

# Função para criar gráfico combinado de variância explicada e acumulada
plot_variance_explained <- function(var_explained, cumulative_variance, name, top_n = NULL) {
  n <- ifelse(is.null(top_n), length(var_explained), top_n)
  df <- data.frame(
    PC = 1:n,
    VarianceExplained = var_explained[1:n],
    CumulativeVariance = cumulative_variance[1:n]
  )
  
  ggplot(df, aes(x = PC)) +
    geom_bar(aes(y = VarianceExplained), stat = "identity", fill = "steelblue", alpha = 0.8) +
    geom_line(aes(y = CumulativeVariance), color = "darkred", size = 1) +
    geom_point(aes(y = CumulativeVariance), color = "darkred", size = 2) +
    labs(
      title = paste("Explained Variance and Cumulative Variance -", name),
      x = "Principal Component (PC)",
      y = "Variance Explained / Cumulative Variance"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    ) +
    scale_y_continuous(sec.axis = sec_axis(~ ., name = "Cumulative Variance"))
}

# Criar gráficos de variância explicada e acumulada para todos os PCs e top 10
plot_variance_r0 <- plot_variance_explained(pca_r0$var_explained, pca_r0$cumulative_variance, "R0")
plot_variance_top10_r0 <- plot_variance_explained(pca_r0$var_explained, pca_r0$cumulative_variance, "R0", top_n = 10)

plot_variance_r30 <- plot_variance_explained(pca_r30$var_explained, pca_r30$cumulative_variance, "R30")
plot_variance_top10_r30 <- plot_variance_explained(pca_r30$var_explained, pca_r30$cumulative_variance, "R30", top_n = 10)

# Salvar gráficos de variância explicada e acumulada
ggsave("images/variance_explained_R0.png", plot_variance_r0, dpi = 300, width = 8, height = 6)
ggsave("images/variance_explained_top10_R0.png", plot_variance_top10_r0, dpi = 300, width = 8, height = 6)
ggsave("images/variance_explained_R30.png", plot_variance_r30, dpi = 300, width = 8, height = 6)
ggsave("images/variance_explained_top10_R30.png", plot_variance_top10_r30, dpi = 300, width = 8, height = 6)

# Função para criar gráficos de dispersão entre PCs
plot_pcs <- function(pca, name) {
  scores <- as.data.frame(pca$x)  # Componentes principais
  combinations <- list(
    c("PC1", "PC2"),
    c("PC1", "PC3"),
    c("PC1", "PC4"),
    c("PC2", "PC3"),
    c("PC2", "PC4"),
    c("PC3", "PC4")
  )
  
  plots <- list()
  for (comb in combinations) {
    p <- ggplot(scores, aes_string(x = comb[1], y = comb[2])) +
      geom_point(alpha = 0.7, color = "blue", size = 1.5) +
      labs(
        title = paste(name, "-", comb[1], "vs", comb[2]),
        x = comb[1],
        y = comb[2]
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)
      )
    plots[[paste(comb[1], comb[2], sep = "_")]] <- p
  }
  
  return(plots)
}

# Criar gráficos de dispersão dos PCs
plots_r0 <- plot_pcs(pca_r0$pca, "R0")
plots_r30 <- plot_pcs(pca_r30$pca, "R30")

# Salvar gráficos de dispersão
for (name in names(plots_r0)) {
  ggsave(paste0("images/", name, "_R0.png"), plots_r0[[name]], dpi = 300, width = 8, height = 6)
}

for (name in names(plots_r30)) {
  ggsave(paste0("images/", name, "_R30.png"), plots_r30[[name]], dpi = 300, width = 8, height = 6)
}


##### gráficos de PCAs combinados

# Função para criar gráficos combinados de variância explicada e acumulada
plot_combined_variance <- function(var_explained, cumulative_variance, name, top_n = NULL) {
  n <- ifelse(is.null(top_n), length(var_explained), top_n)
  df <- data.frame(
    PC = 1:n,
    VarianceExplained = var_explained[1:n],
    CumulativeVariance = cumulative_variance[1:n]
  )
  
  ggplot(df, aes(x = PC)) +
    geom_bar(aes(y = VarianceExplained), stat = "identity", fill = "steelblue", alpha = 0.8) +
    geom_line(aes(y = CumulativeVariance), color = "darkred", size = 1) +
    geom_point(aes(y = CumulativeVariance), color = "darkred", size = 2) +
    labs(
      title = paste("Variance Explained and Cumulative Variance -", name),
      x = "Principal Component (PC)",
      y = "Variance Explained / Cumulative Variance"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
}

# Gráficos para todos os PCs e top 10 PCs
plot_r0_all <- plot_combined_variance(pca_r0$var_explained, pca_r0$cumulative_variance, "Radius 0")
plot_r0_top10 <- plot_combined_variance(pca_r0$var_explained, pca_r0$cumulative_variance, "Radius 0", top_n = 10)

plot_r30_all <- plot_combined_variance(pca_r30$var_explained, pca_r30$cumulative_variance, "Radius 30")
plot_r30_top10 <- plot_combined_variance(pca_r30$var_explained, pca_r30$cumulative_variance, "Radius 30", top_n = 10)

# Gráfico comparativo para Radius 0 e Radius 30 no mesmo gráfico
plot_combined_top10 <- function(var_explained_r0, var_explained_r30, cumulative_r0, cumulative_r30) {
  df <- data.frame(
    PC = rep(1:10, 2),
    VarianceExplained = c(var_explained_r0[1:10], var_explained_r30[1:10]),
    CumulativeVariance = c(cumulative_r0[1:10], cumulative_r30[1:10]),
    Group = rep(c("Radius 0", "Radius 30"), each = 10)
  )
  
  ggplot(df, aes(x = PC, group = Group, color = Group)) +
    geom_bar(aes(y = VarianceExplained, fill = Group), position = "dodge", stat = "identity", alpha = 0.8) +
    geom_line(aes(y = CumulativeVariance), size = 1) +
    geom_point(aes(y = CumulativeVariance), size = 2) +
    labs(
      title = "Variance Explained and Cumulative Variance - Top 10 PCs (Radius 0 & 30)",
      x = "Principal Component (PC)",
      y = "Variance Explained / Cumulative Variance",
      color = "Dataset",
      fill = "Dataset"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
}

plot_top10_combined <- plot_combined_top10(
  pca_r0$var_explained, pca_r30$var_explained, 
  pca_r0$cumulative_variance, pca_r30$cumulative_variance
)

# Salvar os gráficos
ggsave("images/variance_explained_R0_all.png", plot_r0_all, dpi = 300, width = 8, height = 6)
ggsave("images/variance_explained_R0_top10.png", plot_r0_top10, dpi = 300, width = 8, height = 6)
ggsave("images/variance_explained_R30_all.png", plot_r30_all, dpi = 300, width = 8, height = 6)
ggsave("images/variance_explained_R30_top10.png", plot_r30_top10, dpi = 300, width = 8, height = 6)
ggsave("images/variance_explained_combined_top10.png", plot_top10_combined, dpi = 300, width = 10, height = 6)

rotation_pca_r0 <-pca_r0[["pca"]][["rotation"]][,1:4]
rotation_pca_r30 <-pca_r30[["pca"]][["rotation"]][,1:4]

component_pca_r0 <-pca_r0[["pca"]][["x"]][,1:4]
component_pca_r30 <-pca_r30[["pca"]][["x"]][,1:4]

write.csv(rotation_pca_r0, "rotation_pca_R0_PC1_PC2_PC3_PC4.csv", row.names = TRUE)
write.csv(rotation_pca_r30, "rotation_pca_R30_PC1_PC2_PC3_PC4.csv", row.names = TRUE)

write.csv(component_pca_r0, "components_pca_R0_PC1_PC2_PC3_PC4.csv", row.names = TRUE)
write.csv(component_pca_r30, "components_pca_R30_PC1_PC2_PC3_PC4.csv", row.names = TRUE)

#### DE



