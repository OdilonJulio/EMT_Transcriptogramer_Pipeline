library(transcriptogramer)
library(dplyr)
library(biomaRt)
library(Ropj)

# Import expression matrix
exp <- read.csv("expressionMatrixDayZero.csv")
rownames(exp) <- exp$X
exp$X <- NULL

# See transcriptogramerStep1 documentation
# First column of the dictionary MUST BE THE ENSEMBL PEPTIDE ID
# SECOND COLUMN MUST BE THE SAME ID AS ON THE GENE EXPRESSION MATRIX ROWNAMES
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
dict <- getBM(attributes = c("ensembl_peptide_id", "ensembl_gene_id"),
              mart = ensembl)
dict %>% 
  mutate(ensembl_peptide_id = ifelse(ensembl_peptide_id == "", NA, ensembl_peptide_id)) %>% 
  na.omit() -> dict

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
t <- transcriptogramPreprocess(association = assoc, ordering = ord$Protein, radius = 30)
t <- transcriptogramStep1(object = t, expression = exp, dictionary = dict, nCores = 10)
t <- transcriptogramStep2(object = t, nCores = 10)

# Transcriptogram values are on the t@transcriptogramS2 slot
t@transcriptogramS2
