# quality_control_subset.R
# Versão final compatível com Seurat v5+

library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)

# 1. Carregar dados
load("seurat_140cells.RData")

# 2. Definir sample_names conforme script original
sample_names <- c("notreated-batch1", "notreated-batch2", "TGFbeta1-1day-batch2",
                  "TGFbeta1-2day-batch2", "TGFbeta1-3day-batch2", "TGFbeta1-4day-batch1",
                  "TGFbeta1-8day-batch1")

# 3. Processamento individual IDÊNTICO ao original
seurat_objects <- lapply(sample_names, function(name) {
  subset(combined, subset = orig.ident == name)
})

# 4. QC para cada amostra (ATUALIZADO para Seurat v5+)
filtered_seurat_objects <- lapply(seurat_objects, function(seurat) {
  sce <- as.SingleCellExperiment(seurat)
  sce <- scDblFinder(sce)
  seurat$doublet <- sce$scDblFinder.class
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
  subset(seurat, subset = nFeature_RNA > 500 & percent.mt < 20 & doublet == "singlet")
})

# 5. Solução definitiva atualizada para Seurat v5+
# Passo 1: Extrair todos os genes únicos
all_genes <- unique(unlist(lapply(filtered_seurat_objects, rownames)))

# Passo 2: Padronizar os objetos (versão atualizada)
std_objects <- lapply(filtered_seurat_objects, function(obj) {
  missing_genes <- setdiff(all_genes, rownames(obj))
  if(length(missing_genes) > 0) {
    missing_matrix <- matrix(0, 
                             nrow = length(missing_genes), 
                             ncol = ncol(obj),
                             dimnames = list(missing_genes, colnames(obj)))
    
    # Versão compatível com Seurat v5+
    new_counts <- rbind(LayerData(obj, layer = "counts"), missing_matrix)
    obj <- CreateSeuratObject(
      counts = new_counts,
      meta.data = obj@meta.data
    )
  }
  return(obj)
})

# 6. Combinar objetos (MESMA ABORDAGEM do original)
filtered_combined <- merge(
  x = std_objects[[1]],
  y = std_objects[2:7],
  add.cell.ids = sample_names,
  project = "FilteredSubset"
)

# 7. Verificação final
cat("\n=== RESUMO FINAL ===\n")
cat("Total de células:", ncol(filtered_combined), "\n") 
cat("Total de genes:", nrow(filtered_combined), "\n")
cat("Distribuição por amostra:\n")
print(table(filtered_combined$orig.ident))

# 8. Salvar
save(filtered_combined, file = "filtered_subset_140cells.RData")
