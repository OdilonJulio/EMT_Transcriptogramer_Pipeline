## ----message = FALSE----------------------------------------------------------
library(transcriptogramer)
t <- transcriptogramPreprocess(association = association, ordering = Hs900)

## ----message = FALSE----------------------------------------------------------
## during the preprocessing

## creating the object and setting the radius as 0
t <- transcriptogramPreprocess(association = association, ordering = Hs900)

## creating the object and setting the radius as 80
t <- transcriptogramPreprocess(association = association, ordering = Hs900,
                               radius = 80)

## ----message = FALSE----------------------------------------------------------
## after the preprocessing

## modifying the radius of an existing Transcriptogram object
radius(object = t) <- 50

## getting the radius of an existing Transcriptogram object
r <- radius(object = t)

## ----message = FALSE, comment = ""--------------------------------------------
oPropertiesR50 <- orderingProperties(object = t, nCores = 1)

## slight change of radius
radius(object = t) <- 80

## this output is partially different comparing to oPropertiesR50
oPropertiesR80 <- orderingProperties(object = t, nCores = 1)

sum(oPropertiesR50$windowModularity)
sum(oPropertiesR80$windowModularity)

## ----message = FALSE----------------------------------------------------------
cProperties <- connectivityProperties(object = t)

## ----message = FALSE----------------------------------------------------------
t <- transcriptogramStep1(object = t, expression = GSE9988,
                          dictionary = GPL570, nCores = 1)
t2 <- t

## ----message = FALSE----------------------------------------------------------
t <- transcriptogramStep2(object = t, nCores = 1)

## ----message = FALSE----------------------------------------------------------
radius(object = t2) <- 50
t2 <- transcriptogramStep2(object = t2, nCores = 1)

## ----message = FALSE, fig.show = "hide"---------------------------------------
## trend = FALSE for microarray data or voom log2-counts-per-million
## the default value for trend is FALSE
levels <- c(rep(FALSE, 3), rep(TRUE, 3))
t <- differentiallyExpressed(object = t, levels = levels, pValue = 0.01,
                             trend = FALSE, title = "radius 80")

## the radius 50 will affect the output significantly
t2 <- differentiallyExpressed(object = t2, levels = levels, pValue = 0.01,
                              species = DEsymbols, title = "radius 50")

## ----eval = FALSE-------------------------------------------------------------
#  ## using the species argument to translate ENSEMBL Peptide IDs to Symbols
#  ## Internet connection is required for this command
#  t <- differentiallyExpressed(object = t, levels = levels, pValue = 0.01,
#                               species = "Homo sapiens", title = "radius 80")
#  
#  ## translating ENSEMBL Peptide IDs to Symbols using the DEsymbols dataset
#  t <- differentiallyExpressed(object = t, levels = levels, pValue = 0.01,
#                               species = DEsymbols, title = "radius 80")

## ----message = FALSE, comment = ""--------------------------------------------
DE <- DE(object = t)
DE2 <- DE(object = t2)
nrow(DE)
nrow(DE)/length(unique(DE$ClusterNumber))
nrow(DE2)
nrow(DE2)/length(unique(DE2$ClusterNumber))

## ----eval = FALSE-------------------------------------------------------------
#  rdp <- clusterVisualization(object = t)

## ----message = FALSE----------------------------------------------------------
## using the HsBPTerms dataset to create the Protein2GO data.frame
t <- clusterEnrichment(object = t, species = HsBPTerms,
                       pValue = 0.005, nCores = 1)

## ----eval = FALSE-------------------------------------------------------------
#  ## using the species argument to create the Protein2GO data.frame
#  ## Internet connection is required for this command
#  t <- clusterEnrichment(object = t, species = "Homo sapiens",
#                             pValue = 0.005, nCores = 1)

## ----eval = FALSE-------------------------------------------------------------
#  head(Terms(t))

## ----echo = FALSE, comment = ""-----------------------------------------------
load("terms.RData")
head(terms)

## ----eval = FALSE-------------------------------------------------------------
#  enrichmentPlot(t)

## -----------------------------------------------------------------------------
sessionInfo()

## -----------------------------------------------------------------------------
warnings()
