# main.R
# Script principal que chama todos os outros scripts.

# Carregar pacotes
library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)
library(transcriptogramer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(ggplot2)
library(patchwork)

# Executar scripts
source("EMT_01_load_data.R")
source("EMT_02_quality_control.R")
source("EMT_03_normalization.R")
source("EMT_04_batch_effect_correction.R")
source("EMT_05_transcriptogram_analysis.R")
source("EMT_06_pca_analysis.R")
source("EMT_07_reconstruction_analysis.R")
source("EMT_08_r2_analysis.R")
source("EMT_09_plot_transcriptogram.R")
source("EMT_10_line_plots.R")

# Mensagem de conclusão
cat("Análise concluída com sucesso!\n")