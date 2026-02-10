### 1. Pacotes e Ambiente ----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pacotes_bioc <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE", "ReactomePA", "STRINGdb", "biomaRt")
pacotes_cran <- c("ggplot2", "dplyr", "tidyr", "readr")

BiocManager::install(pacotes_bioc, ask = FALSE)
install.packages(pacotes_cran)

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ReactomePA)
library(DOSE)
library(STRINGdb)
library(biomaRt)
library(dplyr)
library(ggplot2)

### 2. Lista de Proteínas ENSP ----
ensp_ids <- genes_in_both_groups$GeneID

### 3. Conversão ENSP → Entrez ID + Gene Symbol ----
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

conversion_df <- getBM(
  attributes = c("ensembl_peptide_id", "entrezgene_id", "external_gene_name"),
  filters = "ensembl_peptide_id",
  values = ensp_ids,
  mart = mart
) %>% filter(!is.na(entrezgene_id))

# Remover duplicatas
conversion_df <- distinct(conversion_df, entrezgene_id, .keep_all = TRUE)

# Vetor de Entrez IDs
entrez_ids <- conversion_df$entrezgene_id

### 4. Enriquecimento Funcional ----

#### 4.1 GO: Processos Biológicos (BP)
ego_bp <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

#### 4.2 KEGG
ekegg <- enrichKEGG(
  gene = entrez_ids,
  organism = "hsa",
  pvalueCutoff = 0.05
)

#### 4.3 Reactome
ereact <- enrichPathway(
  gene = entrez_ids,
  organism = "human",
  pvalueCutoff = 0.05,
  readable = TRUE
)

### 5. Visualização dos Resultados ----

#### Dotplot GO
dotplot(ego_bp, showCategory = 10) +
  ggtitle("GO: Processos Biológicos Enriquecidos")

#### Rede GO
emapplot(ego_bp, showCategory = 15)

#### Cnetplot com genes
cnetplot(ego_bp, categorySize = "geneNum", foldChange = NULL)

### 6. STRINGdb: Rede de Interações ----
string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400)
string_mapped <- string_db$map(conversion_df, "external_gene_name", removeUnmappedRows = TRUE)
string_db$plot_network(string_mapped$STRING_id)

### 7. Exportação ----
write_csv(conversion_df, "genes_emt_convertidos.csv")
write_csv(ego_bp@result, "enriquecimento_GO_BP.csv")
write_csv(ekegg@result, "enriquecimento_KEGG.csv")
write_csv(ereact@result, "enriquecimento_Reactome.csv")
