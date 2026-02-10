# reduced_normalization.R
# Normaliza os dados do subset de 70 células

library(Seurat)

# 1. Carregar objeto combinado filtrado
load("filtered_subset_140cells.RData")

# 2. Extrair camadas (atualizado para Seurat v5+)
layers <- Layers(filtered_combined, assay = "RNA")

# 3. Combinar matrizes de expressão (versão eficiente para subset)
count_matrices <- lapply(layers, function(layer) {
  LayerData(filtered_combined, assay = "RNA", layer = layer)
})

combined_expression <- do.call(cbind, count_matrices)

# 4. Verificação rápida
cat("\n=== DIMENSÕES DA MATRIZ ===\n")
print(dim(combined_expression))
cat("Total de células:", ncol(combined_expression), "\n")
cat("Total de genes:", nrow(combined_expression), "\n")

# 5. Normalização Total Count (mesmo método do original)
# Converter para matriz densa (mais eficiente para pequenos datasets)
expr_matrix <- as.matrix(combined_expression)

# Calcular fatores de normalização
size_factors <- colSums(expr_matrix)

# Aplicar normalização
normalized_matrix <- sweep(expr_matrix, 2, size_factors, "/")

# Verificação
cat("\n=== VERIFICAÇÃO DA NORMALIZAÇÃO ===\n")
cat("Somas das colunas após normalização:\n")
print(summary(colSums(normalized_matrix)))

# 6. Salvar resultados
save(normalized_matrix, file = "normalized_subset_matrix.RData")
cat("\nMatriz normalizada salva como 'normalized_subset_matrix.RData'\n")

# 7. Opcional: Criar objeto Seurat normalizado
seurat_normalized <- CreateSeuratObject(
  counts = expr_matrix,
  meta.data = filtered_combined@meta.data
) 
seurat_normalized <- SetAssayData(
  seurat_normalized,
  layer = "data",
  new.data = normalized_matrix
)
save(seurat_normalized, file = "seurat_normalized_subset.RData")
