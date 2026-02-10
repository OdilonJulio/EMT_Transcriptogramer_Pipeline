# reduced_transcriptogram_analysis_final.R
# Versão corrigida com tratamento de tipos de dados


library(transcriptogramer)
library(dplyr)
library(Ropj)

# Carregar matriz com efeito de lote corrigido.
load("corrected_subset_matrix.RData")

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
  values = rownames(corrected_matrix),  # Usar os nomes dos genes
  mart = ensembl
)

# Limpar o mapeamento removendo valores vazios e NAs
gene_mapping <- gene_mapping %>%
  mutate(ensembl_gene_id = ifelse(ensembl_gene_id == "", NA, ensembl_gene_id)) %>%
  na.omit()  # Remover linhas com NA em 'ensembl_gene_id'

# Alinhar os nomes das linhas com o mapeamento
mapped_genes <- gene_mapping$ensembl_gene_id[match(
  rownames(corrected_matrix),
  gene_mapping$external_gene_name
)]

# Filtrar para manter apenas as linhas com mapeamentos válidos
valid_indices <- !is.na(mapped_genes)
corrected_matrix <- corrected_matrix[valid_indices, , drop = FALSE]  # Filtrar linhas válidas
rownames(corrected_matrix) <- mapped_genes[valid_indices]  # Substituir os nomes das linhas

# Ordering e Association matrix
ord <- vroom("ordering_HomoSapiensScore800-2024-C.txt")
assoc <- read.opj("Associationmatrix.opj")
assoc <- assoc$associationMa

inner_join(assoc, ord, by = c("A" = "dim1")) %>% 
  dplyr::rename(protein1 = Protein) %>% 
  dplyr::inner_join(ord, by = c("B" = "dim1")) %>% 
  dplyr::rename(protein2 = Protein) %>% 
  dplyr::select(protein1, protein2) -> assoc

# Pré-processamento com raio = 0
t_matrix_R0 <- transcriptogramPreprocess(association = assoc, ordering = ord$Protein, radius = 0)
t_matrix_R0 <- transcriptogramStep1(object = t_matrix_R0, expression = corrected_matrix, dictionary = dictionary)
t_matrix_R0 <- transcriptogramStep2(object = t_matrix_R0)

# Pré-processamento com raio = 30
t_matrix_R30 <- transcriptogramPreprocess(association = assoc, ordering = ord$Protein, radius = 30)
t_matrix_R30 <- transcriptogramStep1(object = t_matrix_R30, expression = corrected_matrix, dictionary = dictionary)
t_matrix_R30 <- transcriptogramStep2(object = t_matrix_R30)

# Salvar objetos processados
save(t_matrix_R0, file = "t_matrix_R0_reduced.RData")
save(t_matrix_R30, file = "t_matrix_R30_reduced.RData")
