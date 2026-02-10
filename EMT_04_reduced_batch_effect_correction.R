# reduced_batch_effect_correction.R
# Versão reduzida para 70 células do script de correção de batch effect

# Carregar matriz normalizada do subset
load("normalized_subset_matrix.RData")  # Carrega normalized_matrix

# 1. Identificar colunas para cada batch no subset
cols_batch1 <- grep("notreated_batch1", colnames(normalized_matrix))
cols_batch2 <- grep("notreated_batch2", colnames(normalized_matrix))

# 2. Calcular médias por batch
mean_batch1 <- rowMeans(normalized_matrix[, cols_batch1, drop = FALSE])
mean_batch2 <- rowMeans(normalized_matrix[, cols_batch2, drop = FALSE])

# 3. Calcular razão entre batches (com proteção contra divisão por zero)
ratio <- ifelse(mean_batch1 == 0, 0, mean_batch2 / mean_batch1)

# 4. Aplicar correção multiplicando as células do batch1 pela razão
corrected_matrix <- normalized_matrix
corrected_matrix[, cols_batch1] <- sweep(normalized_matrix[, cols_batch1, drop = FALSE], 
                                         1, ratio, "*")

# 5. Salvar resultados
save(corrected_matrix, file = "corrected_subset_matrix.RData")

# Mensagem de conclusão
cat("Correção de batch effect concluída para o subset de células.\n")
cat("Matriz corrigida salva como 'corrected_subset_matrix.RData'\n")
