a <-  rownames(expression_matrix_MCF10A_forDay)
sum(a %in% dictionary$ensembl_gene_id)
a <- a |> as.data.frame()
colnames(a)[1] <- "ENSEMBL_ID"

dictionary2 <- dictionary[dictionary$ensembl_gene_id %in% a$ENSEMBL_ID,]

sum(duplicated(dictionary2$ensembl_peptide_id))


b <- expression_matrix_MCF10A_forDay 

b <- b |> rownames_to_column("ENSEMBL_ID")

c <-  merge(b, dictionary2, by.x = "ENSEMBL_ID", by.y = "ensembl_gene_id", all.x = T)

# url <- "https://ftp.ensembl.org/pub/release-111/gff3/homo_sapiens/Homo_sapiens.GRCh38.111.gff3.gz"
# destfile <- "Homo_sapiens.GRCh38.111.gff3.gz"
# 
# download.file(url, destfile)
# 
# library("rtracklayer")
# 
# gff <- import.gff3("Homo_sapiens.GRCh38.111.gff3")
# gff<- as.data.frame(gff@elementMetadata) 
# #Filtrar tres colunas do gff
# gff <- gff[,c("type" ,"ID", "Name", "biotype")]
# 
# table(gff$type)
# 
# #Tirar linhas duplicadas
# gff <- dplyr::distinct(gff)
