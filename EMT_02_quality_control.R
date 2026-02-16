# quality_control.R
# Realiza o controle de qualidade (QC) e filtra os dados com plots em qualidade de artigo.

library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)
library(ggplot2)
library(patchwork)

# Carregar objetos Seurat
load("seurat_objects.RData")

# --- 1. Plots pré-QC com todas as amostras no mesmo gráfico ---
for (i in seq_along(seurat_objects)) {
  seurat_objects[[i]]$sample <- sample_names[i]
}
combined_preQC <- merge(seurat_objects[[1]], y = seurat_objects[-1])

# Calcular percent.mt antes do filtro se necessário
if (!"percent.mt" %in% colnames(combined_preQC@meta.data)) {
  combined_preQC[["percent.mt"]] <- PercentageFeatureSet(combined_preQC, pattern = "^MT-")
}

# --- 1. Plots pré-QC com todas as amostras no mesmo gráfico ---
# Gerar violin plots como uma lista separada
plot_list_pre <- VlnPlot(
  combined_preQC,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "sample",
  pt.size = 0, # Remove pontos originais grandes
  combine = FALSE
)

# Adicionar pontos pretos (poeira) e manter títulos
plot_list_pre <- lapply(plot_list_pre, function(x) {
  x + 
    geom_jitter(height = 0, width = 0.2, size = 0.1, alpha = 0.1, color = "black") + # PONTOS PRETOS
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold") # Mantém e centraliza o título (nFeature, etc)
    )
})

# Juntar tudo, coletando a legenda para ficar única na direita
p_pre <- patchwork::wrap_plots(plot_list_pre, ncol = 3) + 
  patchwork::plot_layout(guides = "collect") + # UNIFICA A LEGENDA COLORIDA
  patchwork::plot_annotation(
    title = "Pre-QC Distributions - All Samples",
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  )

# Salvar PDF de alta qualidade
pdf("images/QC_pre_QC.pdf", width = 14, height = 6)
print(p_pre)
dev.off()


# --- 2. Identificar doublets e porcentagem de genes mitocondriais ---
seurat_objects <- lapply(seurat_objects, function(seurat) {
  sce <- as.SingleCellExperiment(seurat)
  sce <- scDblFinder(sce)
  seurat$doublet <- sce$scDblFinder.class
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
  seurat <- subset(seurat, subset = doublet == "singlet")
  return(seurat)
})

# --- 3. Aplicar filtros de QC ---
filtered_seurat_objects <- lapply(seurat_objects, function(seurat) {
  subset(seurat, subset = nFeature_RNA > 500 & percent.mt < 20)
})

# --- 4. Plots pós-QC para cada amostra ---
for (i in seq_along(filtered_seurat_objects)) {
  filtered_seurat_objects[[i]]$sample <- sample_names[i]
}
combined_postQC <- merge(filtered_seurat_objects[[1]], y = filtered_seurat_objects[-1])

# Garantir que percent.mt está presente
if (!"percent.mt" %in% colnames(combined_postQC@meta.data)) {
  combined_postQC[["percent.mt"]] <- PercentageFeatureSet(combined_postQC, pattern = "^MT-")
}

# --- 4. Plots pós-QC para cada amostra ---
# Gerar violin plots pós-QC como lista
plot_list_post <- VlnPlot(
  combined_postQC,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "sample",
  pt.size = 0, 
  combine = FALSE
)

# Adicionar pontos pretos e manter títulos
plot_list_post <- lapply(plot_list_post, function(x) {
  x + 
    geom_jitter(height = 0, width = 0.2, size = 0.1, alpha = 0.1, color = "black") + # PONTOS PRETOS
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold") # Mantém e centraliza o título
    )
})

# Juntar tudo, coletando a legenda
p_post <- patchwork::wrap_plots(plot_list_post, ncol = 3) + 
  patchwork::plot_layout(guides = "collect") + # UNIFICA A LEGENDA COLORIDA
  patchwork::plot_annotation(
    title = "Post-QC Distributions - All Samples",
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  )

# Salvar PDF de alta qualidade

pdf("images/QC_post_QC.pdf", width = 14, height = 6)
print(p_post)
dev.off()


# --- 5. Comparação antes e depois (todas as amostras no mesmo gráfico) ---


# Unir todos os dados pré e pós-QC em um único data.frame
comparison_df <- do.call(rbind, lapply(seq_along(seurat_objects), function(i) {
  pre <- FetchData(seurat_objects[[i]], vars = c("nFeature_RNA", "percent.mt"))
  pre$status <- "Pre-QC"
  pre$sample <- sample_names[i]
  
  post <- FetchData(filtered_seurat_objects[[i]], vars = c("nFeature_RNA", "percent.mt"))
  post$status <- "Post-QC"
  post$sample <- sample_names[i]
  
  rbind(pre, post)
}))

# Plot combinando todas as amostras
p_comp_nf <- ggplot(comparison_df, aes(x = sample, y = nFeature_RNA, fill = status)) +
  geom_violin(scale = "width", trim = TRUE) +
  scale_y_log10() +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sample", y = "Genes detected (log10)", fill = "Status")

p_comp_mt <- ggplot(comparison_df, aes(x = sample, y = percent.mt, fill = status)) +
  geom_violin(scale = "width", trim = TRUE) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Amostra", y = "% Mitocondrial", fill = "Status")

# Combinar e adicionar título
p_comp <- (p_comp_nf | p_comp_mt) + plot_annotation(
  title = "Pre-QC vs Post-QC Comparison – All Samples"
)

# Salvar PDF único de alta resolução
pdf("images/QC_comparison.pdf", width = 14, height = 6)
print(p_comp)
dev.off()


# --- 6. Combinar os objetos filtrados ---
filtered_combined_seurat <- merge(
  filtered_seurat_objects[[1]],
  y = filtered_seurat_objects[-1],
  add.cell.ids = sample_names,
  project = "FilteredCombined"
)

# Salvar objeto combinado
save(filtered_combined_seurat, file = "filtered_combined_seurat.RData")
