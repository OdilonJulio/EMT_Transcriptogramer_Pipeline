# quality_control_subset.R
# Versão completamente corrigida para Seurat v5

library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)
library(Matrix)

# 1. Carregar dados do script anterior
load("seurat_3500cells_random.RData")

# 2. Definir sample_names conforme script original
sample_names <- c("notreated-batch1", "notreated-batch2", "TGFbeta1-1day-batch2",
                  "TGFbeta1-2day-batch2", "TGFbeta1-3day-batch2", "TGFbeta1-4day-batch1",
                  "TGFbeta1-8day-batch1")

# 3. Processamento individual para cada amostra (atualizado para v5)
seurat_objects <- lapply(sample_names, function(name) {
  subset(combined, subset = orig.ident == name)
})

# 4. Controle de qualidade robusto para cada amostra
filtered_seurat_objects <- lapply(seurat_objects, function(seurat) {
  # Converter para SCE (método compatível com v5)
  counts <- LayerData(seurat, layer = "counts")
  sce <- SingleCellExperiment(list(counts = counts))
  colData(sce) <- DataFrame(seurat@meta.data)
  
  # Detecção de dupletos
  sce <- scDblFinder(sce)
  
  # Adicionar metadados
  seurat$doublet_score <- sce$scDblFinder.score
  seurat$doublet_class <- sce$scDblFinder.class
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
  seurat[["percent.rb"]] <- PercentageFeatureSet(seurat, pattern = "^RP[SL]")
  
  # Filtros rigorosos
  seurat <- subset(seurat,
                   subset = nFeature_RNA > 500 & 
                     nFeature_RNA < 6000 &
                     percent.mt < 15 &
                     percent.rb < 50 &
                     doublet_class == "singlet")
  
  return(seurat)
})

# 5. Padronização dos objetos (versão otimizada para v5)
all_genes <- unique(unlist(lapply(filtered_seurat_objects, function(obj) {
  rownames(LayerData(obj, layer = "counts"))
})))

std_objects <- lapply(filtered_seurat_objects, function(obj) {
  obj_counts <- LayerData(obj, layer = "counts")
  missing_genes <- setdiff(all_genes, rownames(obj_counts))
  
  if(length(missing_genes) > 0) {
    missing_matrix <- Matrix(0,
                             nrow = length(missing_genes),
                             ncol = ncol(obj),
                             sparse = TRUE,
                             dimnames = list(missing_genes, colnames(obj)))
    
    combined_counts <- rbind(obj_counts, missing_matrix)
    combined_counts <- combined_counts[all_genes, ]  # Manter ordem
    
    # Criar novo objeto Seurat mantendo todos os metadados
    new_obj <- CreateSeuratObject(
      counts = combined_counts,
      meta.data = obj@meta.data
    )
    return(new_obj)
  }
  return(obj)
})

# 6. Combinar objetos de forma otimizada (atualizado para v5)
filtered_combined <- merge(
  x = std_objects[[1]],
  y = std_objects[2:7],
  add.cell.ids = sample_names,
  project = "Filtered_3500cells"
)

# 7. Análise de qualidade pós-filtro
cat("\n=== RELATÓRIO DE QUALIDADE ===\n")
cat("Células totais após filtro:", ncol(filtered_combined), "\n")
cat("Genes totais:", nrow(filtered_combined), "\n\n")

cat("Distribuição por amostra:\n")
print(table(filtered_combined$orig.ident))

cat("\nMétricas médias pós-filtro:\n")
print(aggregate(cbind(nFeature_RNA, nCount_RNA, percent.mt) ~ orig.ident, 
                data = filtered_combined@meta.data, 
                FUN = median))

# 8. Salvar objeto filtrado
output_file <- "filtered_3500cells_QC_seuratv5.RData"
save(filtered_combined, file = output_file)
cat("\nObjeto salvo como:", output_file, "\n")

# 9. Visualização rápida (opcional)
if(interactive()) {
  VlnPlot(filtered_combined,
          features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
          group.by = "orig.ident",
          pt.size = 0.1) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}