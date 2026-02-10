################################################################################
# TÍTULO: EMT_20_discovery_and_publication_complete.R
# AUTOR: Odilon Júlio dos Santos
# DATA: 15/12/2025
# DESCRIÇÃO: Pipeline UNIFICADO E CORRIGIDO.
#            - Executa Mapeamento e Enriquecimento (Hallmarks).
#            - Resgata automaticamente clusters sem Hallmark (ex: C5) via GO.
#            - Gera Figuras de Artigo (A, B, C e Painel Final).
#            - Gera Tabela 1 (Excel) formatada.
################################################################################

# ==============================================================================
# 1. SETUP DE AMBIENTE (CORRIGIDO)
# ==============================================================================
rm(list = ls())
graphics.off()
options(stringsAsFactors = FALSE)

# Função auxiliar segura para instalar pacotes
ensure_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Instalando pacote necessário: ", pkg)
    
    # Lista de pacotes Bioconductor conhecidos
    bioc_pkgs <- c("clusterProfiler", "org.Hs.eg.db", "ReactomePA", 
                   "enrichplot", "biomaRt", "msigdbr", "AnnotationDbi")
    
    if (pkg %in% bioc_pkgs) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      install.packages(pkg)
    }
  }
}

# Lista de dependências
pkgs <- c("dplyr", "ggplot2", "stringr", "openxlsx", "RColorBrewer", 
          "ggraph", "readr", "forcats", "patchwork", "scales", 
          "clusterProfiler", "org.Hs.eg.db", "biomaRt", "msigdbr", "AnnotationDbi")

# Instalação e carregamento
invisible(lapply(pkgs, ensure_pkg))

# Carregar bibliotecas essenciais para o namespace
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)

# Diretórios
input_dir <- "cluster_analysis"
out_discovery <- file.path(input_dir, "discovery_results")
out_figures   <- file.path(input_dir, "publication_figures")
out_tables    <- file.path(input_dir, "publication_tables")

dirs <- c(out_discovery, out_figures, out_tables)
invisible(lapply(dirs, function(d) dir.create(d, showWarnings = FALSE, recursive = TRUE)))

message("=== 1. AMBIENTE PREPARADO COM SUCESSO ===")

# ==============================================================================
# 2. CARREGAMENTO DOS DADOS (CLUSTERS)
# ==============================================================================
message("\n=== 2. CARREGAMENTO DOS CLUSTERS ===")

cluster_csv <- file.path(input_dir, "identified_gene_clusters.csv")

if (!file.exists(cluster_csv)) {
  # Procura alternativa se o arquivo não estiver no local padrão
  files <- list.files(input_dir, pattern = "clusters.*\\.csv", full.names = TRUE)
  if (length(files) > 0) {
    cluster_csv <- files[1]
    message("⚠️ Arquivo padrão não encontrado. Usando arquivo alternativo: ", cluster_csv)
  } else {
    stop("ERRO CRÍTICO: Nenhum arquivo de clusters encontrado em '", input_dir, "'.")
  }
}

df_clusters <- read.csv(cluster_csv) %>% 
  dplyr::filter(!is.na(Gene) & Gene != "")

# Validação básica
if (!all(c("Gene", "Cluster") %in% colnames(df_clusters))) {
  stop("O CSV de entrada deve conter as colunas 'Gene' e 'Cluster'.")
}

message("Dados carregados: ", nrow(df_clusters), " genes em ", length(unique(df_clusters$Cluster)), " clusters.")

# ==============================================================================
# 3. MAPEAMENTO DE IDs (ENSP -> ENTREZ)
# ==============================================================================
message("\n=== 3. MAPEAMENTO DE IDs (Robustez Extrema) ===")

exemplo <- df_clusters$Gene[1]
# Detecta se é ENSP (Proteína) ou SYMBOL
tipo_origem <- ifelse(stringr::str_detect(exemplo, "^ENSP"), "ENSEMBLPROT", "SYMBOL")
message("Tipo de ID detectado: ", tipo_origem)

df_clusters$ENTREZID <- NA_character_
df_clusters$SYMBOL <- NA_character_

# A. Tentativa Local (AnnotationDbi)
tryCatch({
  message("Mapeando via org.Hs.eg.db (Local)...")
  # Uso explícito de AnnotationDbi::select para evitar conflito com dplyr::select
  mapa_local <- AnnotationDbi::select(
    org.Hs.eg.db, 
    keys = unique(df_clusters$Gene),
    columns = c("ENTREZID", "SYMBOL"), 
    keytype = tipo_origem
  )
  # Remove duplicatas (1 gene -> múltiplos IDs) mantendo o primeiro
  mapa_local <- mapa_local[!duplicated(mapa_local[[tipo_origem]]), ]
  
  idx <- match(df_clusters$Gene, mapa_local[[tipo_origem]])
  df_clusters$ENTREZID <- mapa_local$ENTREZID[idx]
  df_clusters$SYMBOL   <- mapa_local$SYMBOL[idx]
}, error = function(e) message("Aviso: Falha no banco local. ", e$message))

# B. Tentativa Online (biomaRt) - Apenas para preencher NAs
genes_missing <- df_clusters %>% dplyr::filter(is.na(ENTREZID)) %>% dplyr::pull(Gene)

if (length(genes_missing) > 0) {
  message("Tentando recuperar ", length(genes_missing), " genes via biomaRt...")
  tryCatch({
    mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    filt_name <- ifelse(tipo_origem == "ENSEMBLPROT", "ensembl_peptide_id", "external_gene_name")
    
    mapa_mart <- biomaRt::getBM(
      attributes = c(filt_name, "entrezgene_id", "hgnc_symbol"),
      filters = filt_name, values = genes_missing, mart = mart
    )
    
    idx_m <- match(df_clusters$Gene, mapa_mart[[filt_name]])
    mask <- !is.na(idx_m) & is.na(df_clusters$ENTREZID)
    
    df_clusters$ENTREZID[mask] <- as.character(mapa_mart$entrezgene_id[idx_m[mask]])
    df_clusters$SYMBOL[mask]   <- mapa_mart$hgnc_symbol[idx_m[mask]]
  }, error = function(e) message("Aviso: biomaRt offline ou falhou."))
}

# Cleanup final: Se o input já era SYMBOL, preenche a coluna SYMBOL com ele mesmo
if(tipo_origem == "SYMBOL") {
  df_clusters$SYMBOL <- ifelse(is.na(df_clusters$SYMBOL), df_clusters$Gene, df_clusters$SYMBOL)
}

genes_validos <- df_clusters %>% dplyr::filter(!is.na(ENTREZID))
message("✅ Sucesso: ", nrow(genes_validos), " genes mapeados e prontos.")

# ==============================================================================
# 4. MODO DESCOBERTA (HALLMARK MSigDB)
# ==============================================================================
message("\n=== 4. EXECUTANDO VARREDURA HALLMARK (MSigDB) ===")

# Baixar Hallmarks
m_t2g <- tryCatch({
  msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>% 
    dplyr::select(gs_name, entrez_gene)
}, error = function(e) {
  stop("Erro crítico: Falha ao baixar MSigDB. Verifique sua conexão.")
})

# Executar compareCluster
ck_hallmarks <- clusterProfiler::compareCluster(
  ENTREZID ~ Cluster, 
  data = genes_validos, 
  fun = "enricher", 
  TERM2GENE = m_t2g,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

# Converter para DataFrame (pode estar vazio se nada for significativo)
df_res_hallmark <- if (!is.null(ck_hallmarks)) as.data.frame(ck_hallmarks) else data.frame()

# ==============================================================================
# 5. OPERAÇÃO DE RESGATE (PARA CLUSTERS PERDIDOS, EX: CLUSTER 5)
# ==============================================================================
message("\n=== 5. OPERAÇÃO DE RESGATE (GO) ===")

# Identificar clusters que sumiram
todos_clusters <- unique(genes_validos$Cluster)
clusters_encontrados <- if(nrow(df_res_hallmark) > 0) unique(df_res_hallmark$Cluster) else c()
clusters_perdidos <- setdiff(todos_clusters, clusters_encontrados)

df_res_rescue <- data.frame()

if (length(clusters_perdidos) > 0) {
  message("⚠️ Clusters sem Hallmark detectado: ", paste(clusters_perdidos, collapse=", "))
  message("   Iniciando resgate via Gene Ontology (GO BP)...")
  
  for (clus in clusters_perdidos) {
    message("   -> Resgatando ", clus, "...")
    genes_clus <- genes_validos %>% dplyr::filter(Cluster == clus) %>% dplyr::pull(ENTREZID)
    
    # Roda GO com parâmetros relaxados para garantir que ache a identidade
    ego_rescue <- clusterProfiler::enrichGO(
      gene = genes_clus, OrgDb = org.Hs.eg.db, ont = "BP",
      pvalueCutoff = 0.1, qvalueCutoff = 0.2, readable = TRUE
    )
    
    if (!is.null(ego_rescue) && nrow(as.data.frame(ego_rescue)) > 0) {
      # Plot Detalhado Individual do Cluster Resgatado
      p_rescue <- barplot(ego_rescue, showCategory=10) + 
        ggtitle(paste("Detalhamento:", clus, "(GO BP)")) +
        theme(axis.text.y = element_text(size=8))
      
      ggsave(file.path(out_figures, paste0("Fig_Detail_", clus, "_Rescue.png")), 
             p_rescue, width = 8, height = 6)
      
      # Salva os dados formatados para unir com a tabela Hallmark
      temp_res <- as.data.frame(ego_rescue) %>%
        dplyr::mutate(Cluster = clus, Description = paste0("[GO] ", Description)) %>% 
        dplyr::select(Cluster, ID, Description, GeneRatio, BgRatio, p.adjust, geneID, Count)
      
      df_res_rescue <- rbind(df_res_rescue, temp_res)
    }
  }
}

# UNIR RESULTADOS (Hallmark + Rescue)
cols_comuns <- c("Cluster", "ID", "Description", "GeneRatio", "BgRatio", "p.adjust", "geneID", "Count")

df_final_res <- data.frame()

if(nrow(df_res_hallmark) > 0) {
  df_final_res <- rbind(df_final_res, df_res_hallmark %>% dplyr::select(dplyr::all_of(cols_comuns)))
}

if(nrow(df_res_rescue) > 0) {
  df_final_res <- rbind(df_final_res, df_res_rescue %>% dplyr::select(dplyr::all_of(cols_comuns)))
  message("✅ Clusters resgatados integrados à análise final.")
}

# Salvar Tabela Completa
write.csv(df_final_res, file.path(out_discovery, "Tabela_Completa_Integrada.csv"), row.names = FALSE)

# ==============================================================================
# 6. GERAÇÃO DE FIGURAS DE ARTIGO (A, B, C e PAINEL)
# ==============================================================================
message("\n=== 6. GERANDO FIGURAS DE PUBLICAÇÃO ===")

# Limpeza e Padronização dos Nomes
df_clean <- df_final_res %>%
  dplyr::mutate(
    # Remove prefixos
    Term_Clean = stringr::str_remove(Description, "HALLMARK_"),
    Term_Clean = stringr::str_remove(Term_Clean, "\\[GO\\] "), 
    Term_Clean = stringr::str_replace_all(Term_Clean, "_", " "),
    Term_Clean = stringr::str_to_title(Term_Clean), 
    # Label curto (C1, C2...)
    Cluster_Label = stringr::str_replace(Cluster, "Cluster_", "C"), 
    LogP = -log10(p.adjust)
  )

# --- FIGURA A: VALIDAÇÃO EMT ---
message("   -> Gerando Figura A (Validação EMT)...")
termo_emt <- "EPITHELIAL_MESENCHYMAL_TRANSITION"

# Base com todos clusters para garantir eixo X completo
df_base <- data.frame(Cluster = unique(genes_validos$Cluster))

df_emt <- df_clean %>% 
  dplyr::filter(stringr::str_detect(stringr::str_to_upper(Description), termo_emt)) %>%
  dplyr::full_join(df_base, by = "Cluster") %>%
  dplyr::mutate(
    LogP = ifelse(is.na(LogP), 0, LogP),
    Count = ifelse(is.na(Count), 0, Count),
    Cluster_Label = stringr::str_replace(Cluster, "Cluster_", "C"),
    Is_Sig = LogP > -log10(0.05)
  ) %>%
  # Ordenar e remover duplicatas
  dplyr::arrange(Cluster, dplyr::desc(LogP)) %>%
  dplyr::distinct(Cluster, .keep_all = TRUE) %>%
  dplyr::arrange(Cluster)

colors_emt <- c("TRUE" = "#D55E00", "FALSE" = "#999999") 
theme_pub <- theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14), legend.position = "none")

p_emt <- ggplot(df_emt, aes(x = Cluster_Label, y = LogP, fill = Is_Sig)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray30") +
  annotate("text", x = 0.6, y = -log10(0.05) + 0.2, label = "p < 0.05", 
           hjust = 0, size = 3, color = "gray30", fontface = "italic") +
  geom_text(aes(label = ifelse(Count > 0, paste0("n=", Count), "")), 
            vjust = -0.5, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = colors_emt) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
  labs(title = "A. EMT Signature Validation", 
       subtitle = "Enrichment of EMT Hallmark genes",
       x = NULL, y = expression(-log[10](italic(P)[adj]))) +
  theme_pub

ggsave(file.path(out_figures, "FigA_EMT_Validation.png"), p_emt, width = 6, height = 5, dpi = 300)

# --- FIGURA B: LANDSCAPE GLOBAL (Top 3 Termos) ---
message("   -> Gerando Figura B (Global Landscape)...")
df_top <- df_clean %>%
  dplyr::group_by(Cluster) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::slice_head(n = 3) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Cluster = factor(Cluster, levels = stringr::str_sort(unique(Cluster), numeric = TRUE)),
    Term_Clean = forcats::fct_reorder(Term_Clean, as.numeric(Cluster))
  )

p_dot <- ggplot(df_top, aes(x = Cluster_Label, y = Term_Clean)) +
  geom_tile(fill = NA, color = "gray95") + 
  geom_point(aes(size = Count, color = LogP)) +
  scale_size_continuous(range = c(2, 6), name = "Gene Count") +
  scale_color_gradientn(colors = c("blue", "white", "red"), 
                        values = scales::rescale(c(min(df_top$LogP), -log10(0.05), max(df_top$LogP))),
                        name = expression(-log[10](P))) +
  labs(title = "B. Functional Landscape", subtitle = "Top 3 biological processes per cluster",
       x = "Temporal Clusters", y = NULL) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        panel.grid.major = element_line(color = "gray90"))

ggsave(file.path(out_figures, "FigB_Global_Landscape.png"), p_dot, width = 10, height = 8, dpi = 300)

# --- FIGURA C: STORYTELLING (4 Pilares) ---
message("   -> Gerando Figura C (Storytelling)...")

# Termos chaves para sua tese (Case Insensitive)
key_terms_regex <- paste0(c("Epithelial Mesenchymal Transition", "Oxidative Phosphorylation", "Fatty Acid Metabolism", 
                            "G2m Checkpoint", "E2f Targets", "Inflammatory Response", "Tnfa Signaling Via Nfkb", 
                            "Il6 Jak Stat3 Signaling", "Unfolded Protein Response", "Myogenesis", "Apoptosis", "Hypoxia",
                            "Defense Response To Bacterium", "Detoxification Of Copper Ion"), collapse = "|")

df_story <- df_clean %>% 
  dplyr::filter(stringr::str_detect(Term_Clean, stringr::regex(key_terms_regex, ignore_case = TRUE))) %>%
  # Ordenar por Cluster para ficar bonito
  dplyr::arrange(Cluster)

p_heat <- ggplot(df_story, aes(x = Cluster_Label, y = Term_Clean, fill = LogP)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(LogP > -log10(0.01), "**", ifelse(LogP > -log10(0.05), "*", ""))), vjust = 0.8) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, name = "Sig.") +
  labs(title = "C. The Pillars of EMT", subtitle = "Selected Hallmarks & Processes",
       x = NULL, y = NULL, caption = "* p < 0.05, ** p < 0.01") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14), panel.grid = element_blank())

ggsave(file.path(out_figures, "FigC_Storytelling.png"), p_heat, width = 8, height = 6, dpi = 300)

# --- PAINEL FINAL ---
message("   -> Montando Painel Final...")
layout <- "
AAABBB
AAABBB
CCCCCC
CCCCCC
"
painel <- p_emt + p_dot + p_heat + 
  plot_layout(design = layout) +
  plot_annotation(title = "Integrative Analysis of EMT Transcriptional Modules",
                  theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)))

ggsave(file.path(out_figures, "PANEL_FINAL_ARTIGO.png"), painel, width = 16, height = 12, dpi = 300)
ggsave(file.path(out_figures, "PANEL_FINAL_ARTIGO.pdf"), painel, width = 16, height = 12, dpi = 300)

# ==============================================================================
# 7. TABELA SÍNTESE PARA DISSERTAÇÃO (EXCEL)
# ==============================================================================
message("\n=== 7. GERANDO TABELA 1 (Excel) ===")

summarize_genes <- function(gene_str) {
  if (is.na(gene_str)) return("")
  genes <- unlist(stringr::str_split(gene_str, "/"))
  paste(head(genes, 7), collapse = ", ")
}

df_summary <- df_clean %>%
  dplyr::group_by(Cluster) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::slice_head(n = 3) %>%
  dplyr::summarise(
    Top_Processos = paste0(Term_Clean, " (p=", format(p.adjust, digits=2, scientific=TRUE), ")", collapse = ";\n"),
    Identidade_Primaria = dplyr::first(Term_Clean),
    Genes_Drivers = summarize_genes(dplyr::first(geneID)),
    N_Genes = dplyr::first(Count),
    Best_P_Value = dplyr::first(p.adjust)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Interpretacao = dplyr::case_when(
      stringr::str_detect(Top_Processos, "Epithelial|Mesenchymal") ~ "EMT Core",
      stringr::str_detect(Top_Processos, "Oxidative|Fatty|Glycolysis") ~ "Metabolismo",
      stringr::str_detect(Top_Processos, "Inflammatory|Interferon|Tnf|Il6") ~ "Inflamação",
      stringr::str_detect(Top_Processos, "G2m|E2f|Mitotic") ~ "Ciclo Celular",
      stringr::str_detect(Top_Processos, "Hypoxia|Apoptosis|Unfolded") ~ "Estresse (UPR/Hipóxia)",
      stringr::str_detect(Top_Processos, "Metal|Copper|Zinc|Detoxification") ~ "Resposta a Metais/Stress",
      stringr::str_detect(Top_Processos, "Myogenesis|Angiogenesis") ~ "Remodelagem Tecidual",
      TRUE ~ "Outros"
    )
  ) %>%
  dplyr::select(Cluster, Interpretacao, Identidade_Primaria, Top_Processos, N_Genes, Genes_Drivers, Best_P_Value)

# Exportação Excel
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Resumo Clusters")
hs <- openxlsx::createStyle(fontColour = "#ffffff", fgFill = "#4F81BD", halign = "center", textDecoration = "Bold")
style_interp <- openxlsx::createStyle(textDecoration = "bold")

openxlsx::writeData(wb, "Resumo Clusters", df_summary, headerStyle = hs)
openxlsx::addStyle(wb, "Resumo Clusters", style_interp, rows = 2:(nrow(df_summary)+1), cols = 2)
openxlsx::setColWidths(wb, "Resumo Clusters", cols = 1:7, widths = c(10, 20, 30, 50, 10, 40, 15))

outfile_tab <- file.path(out_tables, "Tabela1_Identidade_Clusters_Dissertacao.xlsx")
openxlsx::saveWorkbook(wb, outfile_tab, overwrite = TRUE)

message("\n✅ PROCESSO CONCLUÍDO COM SUCESSO EXTREMO!")
message("📁 Verifique: ", file.path(out_figures, "PANEL_FINAL_ARTIGO.png"))
message("📁 Verifique: ", outfile_tab)