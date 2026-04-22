################################################################################
# TÍTULO: EMT_21_cell_markers_and_outliers.R
# DESCRIÇÃO: 
#   1. Classifica as células em estágios da EMT utilizando o PC1 como pseudotempo
#      e extrai os marcadores diferenciais de cada estágio.
#   2. Isola os "Outliers de Transição" (células do Dia 8 que resistiram à EMT) 
#      e identifica os marcadores que explicam essa estagnação/resistência.
################################################################################

# ==========================================
# 0. SETUP E CARREGAMENTO DE DADOS
# ==========================================
suppressMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(openxlsx)
  library(Matrix) # Necessário para manipular a matriz esparsa
})

dir.create("markers_analysis", showWarnings = FALSE)

# Carregar o objeto Seurat filtrado e o resultado da PCA Topológica
load("filtered_combined_seurat.RData")
load("pca_result_R30.RData")

# ==========================================
# 1. ALINHAMENTO DA TOPOLOGIA COM O SEURAT E INJEÇÃO DE MATRIZ
# ==========================================
# Extrair scores do PC1 (nosso Proxy de progressão da EMT)
pc1_scores <- pca_result_R30$pca_result$x[, 1]

# Carregar a matriz que já passou pela sua correção de lote (EMT_04)
load("copy_matrix.RData")

# Garantir que estamos usando exatamente as mesmas células em todos os objetos
common_cells <- intersect(intersect(names(pc1_scores), colnames(filtered_combined_seurat)), colnames(copy_matrix))
seurat_obj <- subset(filtered_combined_seurat, cells = common_cells)
matriz_alinhada <- copy_matrix[, common_cells]

# Injetar os scores no metadado do Seurat
seurat_obj$Topological_PC1 <- pc1_scores[common_cells]

# Unificar camadas fragmentadas do Seurat V5 e injetar a matriz corrigida na camada 'data'
seurat_obj <- JoinLayers(seurat_obj)
seurat_obj[["RNA"]]$data <- as(matriz_alinhada, "sparseMatrix")

# ==========================================
# 2. CLASSIFICAÇÃO DAS CÉLULAS (ESTÁGIOS EMT)
# ==========================================
# Dividir o contínuo do PC1 em tercis (3 estágios biológicos)
quantiles_pc1 <- quantile(seurat_obj$Topological_PC1, probs = c(0, 0.33, 0.66, 1))

seurat_obj$EMT_Stage <- case_when(
  seurat_obj$Topological_PC1 <= quantiles_pc1[2] ~ "1_Epithelial",
  seurat_obj$Topological_PC1 > quantiles_pc1[2] & seurat_obj$Topological_PC1 <= quantiles_pc1[3] ~ "2_Hybrid_Transition",
  seurat_obj$Topological_PC1 > quantiles_pc1[3] ~ "3_Mesenchymal"
)

Idents(seurat_obj) <- "EMT_Stage"

# Buscar marcadores positivos para cada estágio
message("Iniciando busca de marcadores para classificação de estágios (isso pode levar alguns minutos)...")
markers_emt_stages <- FindAllMarkers(
  seurat_obj, 
  only.pos = TRUE, 
  min.pct = 0.25, 
  logfc.threshold = 0.25
)

# Filtrar os top 20 marcadores por estágio para relatório
top20_stages <- markers_emt_stages %>% 
  group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

# ==========================================
# 3. IDENTIFICAÇÃO DE MARCADORES DOS OUTLIERS
# ==========================================
# Filtrar apenas as células do Dia 8
day8_cells <- colnames(seurat_obj)[seurat_obj$sample == "TGFbeta1-8day-batch1"]
seurat_day8 <- subset(seurat_obj, cells = day8_cells)

# Definir limite de Outlier: Células do Dia 8 que estão no quartil INFERIOR do PC1 geral
threshold_outlier <- quantile(seurat_day8$Topological_PC1, probs = 0.15)

seurat_day8$Outlier_Status <- ifelse(
  seurat_day8$Topological_PC1 <= threshold_outlier, 
  "Outlier_Resistant", 
  "Fully_Transitioned"
)

Idents(seurat_day8) <- "Outlier_Status"

# Comparar diretamente o Outlier com a célula que completou a transição
message("Iniciando busca de marcadores que explicam a resistência/estagnação dos Outliers do Dia 8...")
outlier_markers <- FindMarkers(
  seurat_day8, 
  ident.1 = "Outlier_Resistant", 
  ident.2 = "Fully_Transitioned",
  logfc.threshold = 0.25
)

# Adicionar a coluna de Gene e classificar por significância
outlier_markers$Gene <- rownames(outlier_markers)
outlier_markers <- outlier_markers %>%
  arrange(p_val_adj, desc(abs(avg_log2FC)))

# ==========================================
# 4. EXPORTAÇÃO E RESULTADOS
# ==========================================
message("Exportando tabelas de marcadores para Excel...")

wb <- createWorkbook()
addWorksheet(wb, "Marcadores_Estagios_EMT")
addWorksheet(wb, "Top20_Estagios")
addWorksheet(wb, "Marcadores_Outliers_Dia8")

writeData(wb, "Marcadores_Estagios_EMT", markers_emt_stages)
writeData(wb, "Top20_Estagios", top20_stages)
writeData(wb, "Marcadores_Outliers_Dia8", outlier_markers)

saveWorkbook(wb, "markers_analysis/EMT_Cell_and_Outlier_Markers.xlsx", overwrite = TRUE)

message("Análise de marcadores concluída com sucesso! Verifique o arquivo Excel na pasta 'markers_analysis'.")