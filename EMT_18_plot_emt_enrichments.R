################################################################################
# plot_emt_enrichments_top30.R
# Versão melhorada e corrigida: seleciona top 30 termos EMT-related e gera múltiplos plots
# Saídas PNG em images/emt_plots/
################################################################################

# --- Configurações ------------------------------------------------------------
base_dirs <- c("enrichment_results", "tables", ".")
out_dir <- file.path("images", "emt_plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# parâmetros
PLOT_TOP_N_TERMS <- 30   # apenas os 30 primeiros termos por score
MAX_GENES_FOR_MATRIX <- 40  # limitar número de genes mostrados na matriz / rede
emt_keywords <- c("epithel", "mesench", "emt", "epithelial", "mesenchymal",
                  "extracell", "extracellular", "matrix", "adhesion",
                  "migration", "wnt", "tgf", "transforming growth", "smad",
                  "focal adhesion", "cell migration", "collagen", "integrin")

# pacotes necessários
required <- c("dplyr","readr","ggplot2","stringr","tibble","forcats","pheatmap",
              "tidyr","scales","igraph","ggraph","ggforce","RColorBrewer")
missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if(length(missing_pkgs) > 0){
  message("Instalando pacotes em falta: ", paste(missing_pkgs, collapse = ", "))
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(stringr);
  library(tibble); library(forcats); library(pheatmap); library(tidyr);
  library(scales); library(igraph); library(ggraph); library(ggforce);
  library(RColorBrewer)
})

# --- utilitários robustos -----------------------------------------------------
find_file <- function(patterns, dirs = base_dirs) {
  for(d in dirs) {
    if(!dir.exists(d)) next
    files <- list.files(d, pattern = patterns, full.names = TRUE, ignore.case = TRUE)
    if(length(files) > 0) return(files[1])
  }
  return(NULL)
}

safe_read_csv <- function(path) {
  if(is.null(path) || !file.exists(path)) return(NULL)
  df <- tryCatch(readr::read_csv(path, show_col_types = FALSE),
                 error = function(e) tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
                                              error = function(e2) NULL))
  return(df)
}

get_pvalue_col <- function(df) {
  if(is.null(df)) return(NULL)
  names_lower <- tolower(names(df))
  i <- grep("^p.adjust$|p_adjust|padj|adj.p|adj_p|qvalue", names_lower)
  if(length(i) > 0) return(names(df)[i[1]])
  i2 <- grep("pvalue|p.value|p_val|p.value", names_lower)
  if(length(i2) > 0) return(names(df)[i2[1]])
  return(NULL)
}

get_desc_col <- function(df) {
  if(is.null(df)) return(NULL)
  names_lower <- tolower(names(df))
  i <- grep("description|term_name|term|name|title", names_lower)
  if(length(i) > 0) return(names(df)[i[1]])
  return(NULL)
}

get_genes_col <- function(df) {
  if(is.null(df)) return(NULL)
  names_lower <- tolower(names(df))
  i <- grep("geneid|gene_id|genes|leadingedge|leadingEdge|intersection|geneids|mapped_genes", names_lower)
  if(length(i) > 0) return(names(df)[i[1]])
  return(NULL)
}

is_emt_term <- function(term_label, keywords = emt_keywords) {
  if(is.null(term_label) || length(term_label) == 0) return(FALSE)
  lab <- tolower(term_label)
  any(vapply(keywords, function(k) grepl(k, lab, perl = TRUE), logical(1)))
}

parse_gene_list <- function(s) {
  if(is.null(s)) return(character(0))
  s <- as.character(s)
  s[s==""] <- NA
  s <- s[!is.na(s)]
  if(length(s) == 0) return(character(0))
  # separadores comuns: "/", ",", ";", "|" ou espaços
  out <- unlist(strsplit(s, "[/;,\\|]+|\\s+"))
  out <- out[nzchar(out)]
  unique(out)
}

# --- localizar arquivos -------------------------------------------------------
f_reactome <- find_file("ReactomePA_enrichPathway\\.csv$")
f_go_bp   <- find_file("GO_BP_clusterProfiler\\.csv$|GO_BP_.*\\.csv$")
f_kegg    <- find_file("KEGG_clusterProfiler_readable\\.csv$|KEGG_.*\\.csv$")
f_gprof   <- find_file("gprofiler2_results\\.csv$|gprofiler_.*\\.csv$")

message("Arquivos detectados (primeiros encontrados):")
message(" Reactome: ", ifelse(is.null(f_reactome), "<none>", f_reactome))
message(" GO BP  : ", ifelse(is.null(f_go_bp), "<none>", f_go_bp))
message(" KEGG   : ", ifelse(is.null(f_kegg), "<none>", f_kegg))
message(" g:Prof : ", ifelse(is.null(f_gprof), "<none>", f_gprof))

df_reactome <- safe_read_csv(f_reactome)
df_go_bp    <- safe_read_csv(f_go_bp)
df_kegg     <- safe_read_csv(f_kegg)
df_gprof    <- safe_read_csv(f_gprof)

# --- extrair termos EMT de cada fonte ----------------------------------------
process_enrich_df <- function(df, src_name) {
  if(is.null(df)) return(NULL)
  desc_col <- get_desc_col(df)
  pcol <- get_pvalue_col(df)
  gcol <- get_genes_col(df)
  if(is.null(desc_col)) return(NULL)
  tdf <- tibble::tibble(
    source = src_name,
    term = as.character(df[[desc_col]]),
    pval_raw = if(!is.null(pcol)) as.numeric(df[[pcol]]) else NA_real_,
    genes_raw = if(!is.null(gcol)) as.character(df[[gcol]]) else NA_character_
  )
  tdf <- tdf %>% mutate(n_genes = ifelse(!is.na(genes_raw),
                                         vapply(genes_raw, function(x) length(parse_gene_list(x)), integer(1)),
                                         NA_integer_))
  tdf_emt <- tdf %>% dplyr::filter(vapply(term, is_emt_term, logical(1)))
  if(nrow(tdf_emt) == 0) return(NULL)
  tdf_emt <- tdf_emt %>% mutate(score = ifelse(!is.na(pval_raw) & pval_raw>0, -log10(pval_raw), NA_real_))
  return(tdf_emt %>% dplyr::select(source, term, pval_raw, score, n_genes, genes_raw))
}

r_reactome_emt <- process_enrich_df(df_reactome, "Reactome")
r_go_emt      <- process_enrich_df(df_go_bp, "GO_BP")
r_kegg_emt    <- process_enrich_df(df_kegg, "KEGG")
r_gprof_emt   <- process_enrich_df(df_gprof, "gProfiler")

emt_terms_all <- dplyr::bind_rows(r_reactome_emt, r_go_emt, r_kegg_emt, r_gprof_emt)
if(is.null(emt_terms_all) || nrow(emt_terms_all) == 0) stop("Nenhum termo EMT detectado nas fontes encontradas.")

# fallback score
emt_terms_all <- emt_terms_all %>%
  mutate(score = ifelse(is.na(score) & !is.na(pval_raw) & pval_raw>0, -log10(pval_raw), score))

if(all(is.na(emt_terms_all$score))) {
  emt_terms_all$score <- emt_terms_all$n_genes
}

# agrupar termos idênticos vindos de fontes diferentes (manter maior score)
emt_terms_collapsed <- emt_terms_all %>%
  group_by(term) %>%
  arrange(desc(score)) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(desc(score))

top_terms <- head(emt_terms_collapsed, PLOT_TOP_N_TERMS)
message("Top terms selecionados: ", nrow(top_terms))

# salvar consolidação
readr::write_csv(emt_terms_all, file.path(out_dir, "emt_terms_all_sources_raw.csv"))
readr::write_csv(top_terms, file.path(out_dir, "emt_terms_topN.csv"))

# --- preparação para plots ---------------------------------------------------
# construir dataframe de pares termo-gene (apenas para top_terms)
term_gene_pairs <- tibble(term = character(0), gene = character(0), source = character(0), score = numeric(0))
for(i in seq_len(nrow(emt_terms_all))) {
  row <- emt_terms_all[i, ]
  if(!(row$term %in% top_terms$term)) next
  genes_vec <- parse_gene_list(row$genes_raw)
  if(length(genes_vec) == 0) next
  tmp <- tibble(term = rep(row$term, length(genes_vec)),
                gene = genes_vec,
                source = rep(as.character(row$source), length(genes_vec)),
                score = rep(row$score, length(genes_vec)))
  term_gene_pairs <- bind_rows(term_gene_pairs, tmp)
}
if(nrow(term_gene_pairs) == 0) stop("Sem pares termo-gene para os termos selecionados.")

# contar genes mais frequentes
gene_freq <- term_gene_pairs %>% count(gene, sort = TRUE)
top_genes <- head(gene_freq$gene, MAX_GENES_FOR_MATRIX)

# filtrar pares para top_genes (evita poluição)
tgp_small <- term_gene_pairs %>% filter(gene %in% top_genes)

# criar matriz term x gene (0/1) para tile plot
mat_df <- tgp_small %>% distinct(term, gene) %>%
  mutate(pres = 1) %>%
  pivot_wider(names_from = gene, values_from = pres, values_fill = 0) %>%
  arrange(desc(term))
if(nrow(mat_df) == 0) stop("Matriz termo-gene vazia após filtragem.")

rownames_mat <- mat_df$term
mat <- as.matrix(mat_df[,-1, drop = FALSE])
rownames(mat) <- rownames_mat

# --- Plot A: Horizontal bar plot (Top terms) --------------------------------
bar_df <- top_terms %>% mutate(term_short = stringr::str_trunc(term, 80)) %>%
  arrange(score) %>%
  mutate(term_short = factor(term_short, levels = term_short))

p_bar <- ggplot(bar_df, aes(x = term_short, y = score)) +
  geom_col(fill = RColorBrewer::brewer.pal(8, "Set2")[1], width = 0.7) +
  coord_flip() +
  labs(title = paste0("Top ", nrow(bar_df), " EMT-related enriched terms"),
       subtitle = "Ordem por intensidade (-log10 p) / score. Fonte: Reactome / GO / KEGG / g:Profiler",
       x = NULL,
       y = expression("Intensidade: " * -log[10](p)),
       caption = "Interpretação: termos com maior barra têm evidência estatística mais forte.") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), plot.caption = element_text(size = 9))

ggsave(file.path(out_dir, "emt_top_terms_bar_horizontal.png"), p_bar, width = 12, height = max(4, 0.3 * nrow(bar_df) + 2), dpi = 600)
message("Saved: emt_top_terms_bar_horizontal.png")

# --- Plot B: Dotplot (# genes vs score) ---------------------------------------
dot_df <- top_terms %>% mutate(term_short = stringr::str_trunc(term, 80)) %>%
  left_join(emt_terms_all %>% group_by(term) %>% summarise(n_genes = max(n_genes, na.rm = TRUE)), by = "term") %>%
  mutate(score_plot = score)

p_dot <- ggplot(dot_df, aes(x = reorder(term_short, score_plot), y = score_plot)) +
  geom_point(aes(size = n_genes, color = score_plot)) +
  coord_flip() +
  scale_color_gradient(low = "lightblue", high = "red", na.value = "grey50") +
  scale_size_continuous(range = c(3,8)) +
  labs(title = paste0("Dotplot: Top ", nrow(dot_df), " EMT-associated terms"),
       subtitle = "Tamanho do ponto = número de genes no termo; cor = intensidade (-log10 p).",
       x = NULL, y = expression(-log[10](p)), caption = "Comparação simultânea de p-value e cobertura do termo (n genes).") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

ggsave(file.path(out_dir, "emt_top_terms_dotplot.png"), p_dot, width = 12, height = max(4, 0.3 * nrow(dot_df) + 2), dpi = 600)
message("Saved: emt_top_terms_dotplot.png")

# --- Plot C: Lollipop plot (ordered) ------------------------------------------
lop_df <- bar_df %>% arrange(score)
p_lolli <- ggplot(lop_df, aes(x = term_short, y = score)) +
  geom_segment(aes(x = term_short, xend = term_short, y = 0, yend = score), color = "grey70", size = 0.4) +
  geom_point(size = 3, color = RColorBrewer::brewer.pal(8, "Set1")[2]) +
  coord_flip() +
  labs(title = paste0("Lollipop: Top ", nrow(lop_df), " EMT-related terms (ordered)"),
       subtitle = "Comprimentos = intensidade (-log10 p).",
       x = NULL, y = expression(-log[10](p)),
       caption = "Lollipop enfatiza magnitude relativa.") +
  theme_minimal(base_size = 12)

ggsave(file.path(out_dir, "emt_top_terms_lollipop.png"), p_lolli, width = 12, height = max(4, 0.3 * nrow(lop_df) + 2), dpi = 600)
message("Saved: emt_top_terms_lollipop.png")

# --- Plot D: term × gene matrix (tile) --------------------------------------
mat_long <- as.data.frame(mat) %>% tibble::rownames_to_column(var = "term")
mat_long <- tidyr::pivot_longer(mat_long, cols = -term, names_to = "gene", values_to = "present")
term_order <- top_terms$term
gene_order <- gene_freq$gene
mat_long$term <- factor(mat_long$term, levels = rev(term_order))
mat_long$gene <- factor(mat_long$gene, levels = gene_order)

p_tile <- ggplot(mat_long, aes(x = gene, y = term, fill = factor(present))) +
  geom_tile(color = "grey70") +
  scale_fill_manual(values = c("0" = "white", "1" = "#2166AC"), labels = c("absent","present"), name = "Presença") +
  labs(title = "Term × Gene (EMT) — Presença de genes nos termos (top terms × top genes)",
       subtitle = paste0("Somente top ", nrow(top_terms), " termos e até ", length(top_genes), " genes mostrados."),
       x = "Gene", y = "Termo (enriquecimento)",
       caption = "Tile preenchido indica que o gene faz parte do termo/enriquecimento. Ordenação por frequência/score.") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 9),
        legend.position = "right")

ggsave(file.path(out_dir, "emt_term_gene_tile_matrix.png"), p_tile, width = max(10, 0.25 * ncol(mat) + 6), height = max(6, 0.25 * nrow(mat) + 4), dpi = 600)
message("Saved: emt_term_gene_tile_matrix.png")

# --- Plot E: bipartite network (corrigido para bipartite) ---------------------
# Build igraph bipartite graph from edges (term <-> gene)
edges <- tgp_small %>% distinct(term, gene) %>% rename(from = term, to = gene)
g <- igraph::graph_from_data_frame(edges, directed = FALSE)

# Explicit bipartite type: logical vector required by layout = 'bipartite' (TRUE/FALSE)
term_nodes <- unique(edges$from)
gene_nodes <- unique(edges$to)
all_nodes <- V(g)$name

# create logical type: TRUE for term nodes, FALSE for gene nodes
vtype_logical <- ifelse(all_nodes %in% term_nodes, TRUE,
                        ifelse(all_nodes %in% gene_nodes, FALSE, NA))

# if any NA remain, remove those vertices (defensive) and warn
if(any(is.na(vtype_logical))) {
  warning("Alguns nós não puderam ser classificados como term/gene. Eles serão removidos para garantir layout bipartite.")
  to_remove <- all_nodes[is.na(vtype_logical)]
  g <- igraph::delete_vertices(g, to_remove)
  # recompute node lists and logical types
  all_nodes <- V(g)$name
  vtype_logical <- ifelse(all_nodes %in% term_nodes, TRUE, FALSE)
}

# assign required logical type attribute
V(g)$type <- as.logical(vtype_logical)

# create a human-readable factor for plotting colors
V(g)$vtype_char <- ifelse(V(g)$type, "term", "gene")

# compute degree-based sizes
deg <- igraph::degree(g)
V(g)$size_attr <- scales::rescale(deg, to = c(3,9))

# Plot using bipartite layout (now safe)
set.seed(42)
p_net <- ggraph::ggraph(g, layout = 'bipartite') +
  geom_edge_link(aes(), edge_colour = "grey80", edge_alpha = 0.6) +
  geom_node_point(aes(color = vtype_char, size = size_attr)) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_manual(values = c("term" = "#D95F02", "gene" = "#1B9E77"), name = "") +
  guides(size = FALSE) +
  labs(title = "Rede bipartida: termos EMT ↔ genes (top terms & frequent genes)",
       subtitle = "Nós laranja = termos; nós verdes = genes. Tamanho do nó ~ grau (n conexões).",
       caption = "Interpretação: genes conectados a muitos termos podem ser hubs funcionais.") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "emt_term_gene_bipartite_network.png"), p_net, width = 14, height = 10, dpi = 600)
message("Saved: emt_term_gene_bipartite_network.png")

# --- Plot F: Heatmap de expressão dos genes EMT (se tabela percent disponível) -
# (Aqui deixo a lógica original apenas se df_percent existir; adaptável)
f_percent <- find_file("region_genes_percentile_.*_0_100.*\\.csv$|genes_percentile_matrix_0_100_interpolated.*\\.csv", base_dirs)
df_percent <- safe_read_csv(f_percent)
genes_for_heatmap <- unique(top_genes)  # por padrão usamos os top genes detectados
if(!is.null(df_percent) && nrow(df_percent) > 0 && length(genes_for_heatmap) >= 2) {
  pct <- df_percent
  # assume first column = gene id
  colnames(pct) <- make.names(colnames(pct))
  gene_colname <- names(pct)[1]
  rownames(pct) <- as.character(pct[[gene_colname]])
  numeric_cols <- setdiff(names(pct), gene_colname)
  expr_num <- as.data.frame(lapply(pct[numeric_cols], function(x) as.numeric(as.character(x))), stringsAsFactors = FALSE)
  rownames(expr_num) <- rownames(pct)
  genes_present <- intersect(genes_for_heatmap, rownames(expr_num))
  if(length(genes_present) >= 2) {
    mat_e <- expr_num[genes_present, , drop = FALSE]
    mat_ez <- t(scale(t(as.matrix(mat_e))))
    mat_ez[is.na(mat_ez)] <- 0
    fname <- file.path(out_dir, "emt_genes_expression_heatmap.png")
    pheatmap(mat_ez, cluster_rows = TRUE, cluster_cols = FALSE, fontsize_row = 8, fontsize_col = 9,
             main = "Heatmap: genes EMT (z-scored) ao longo dos percent steps", filename = fname,
             width = 10, height = max(6, 0.25*nrow(mat_ez)+3))
    message("Saved: emt_genes_expression_heatmap.png")
  } else {
    message("Heatmap não gerado: apenas ", length(genes_present), " genes presentes na tabela percent (precisa >=2).")
  }
} else {
  message("Heatmap não gerado: faltam dados percent/região ou genes para heatmap.")
}

# --- salvar resumo e instruções nas figuras -----------------------------------
readme_txt <- file.path(out_dir, "README_EMT_Plots.txt")
if(file.exists(readme_txt)) file.remove(readme_txt)
cat("EMT Plots Summary\n", file = readme_txt)
cat("=================\n\n", file = readme_txt, append = TRUE)
cat("Arquivos gerados (PNG):\n", file = readme_txt, append = TRUE)
cat(" - emt_top_terms_bar_horizontal.png\n", file = readme_txt, append = TRUE)
cat(" - emt_top_terms_dotplot.png\n", file = readme_txt, append = TRUE)
cat(" - emt_top_terms_lollipop.png\n", file = readme_txt, append = TRUE)
cat(" - emt_term_gene_tile_matrix.png\n", file = readme_txt, append = TRUE)
cat(" - emt_term_gene_bipartite_network.png\n", file = readme_txt, append = TRUE)
cat(" - emt_genes_expression_heatmap.png (opcional)\n\n", file = readme_txt, append = TRUE)
cat("Legendas / interpretação rápida:\n", file = readme_txt, append = TRUE)
cat(" - Barras / pontos: intensidade corresponde a -log10(p) (quanto maior, mais significativo)\n", file = readme_txt, append = TRUE)
cat(" - Tamanho do ponto no dotplot: número de genes envolvidos no termo\n", file = readme_txt, append = TRUE)
cat(" - Tile preenchido na matriz: gene presente no termo; colunas genes ordenadas por frequência.\n", file = readme_txt, append = TRUE)
cat(" - Rede bipartida: nós laranja = termos; nós verdes = genes; tamanho do nó ~ grau (n conexões).\n\n", file = readme_txt, append = TRUE)
cat("Recomendações:\n1) Verifique p.adjust/qvalue nos CSVs originais para confirmar robustez estatística.\n2) Use o heatmap de expressão para confirmar direção (up/down) dos genes EMT ao longo dos percent steps.\n", file = readme_txt, append = TRUE)
message("README salvo em: ", readme_txt)

message("Todos os plots gerados em: ", normalizePath(out_dir))
