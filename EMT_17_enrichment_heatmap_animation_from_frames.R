################################################################################
# pipeline_regions_and_enrichment.R
# Robust pipeline to:
#  1) define arbitrary high-variation regions on transcriptogram PCA loadings,
#  2) export region CSVs and plots,
#  3) map IDs (ENSP/ENSG/Entrez/SYMBOL) to Entrez & gene symbol,
#  4) run enrichment (clusterProfiler / ReactomePA / g:Profiler) when sensible,
#  5) save results and sessionInfo for reproducibility.
#
# Run in R (R >= 4.0 recommended). This script will try to install missing
# packages (CRAN vs Bioconductor handled separately) but prefer to prepare
# the environment manually for production runs.
################################################################################

# 0. Settings -------------------------------------------------------------------
options(stringsAsFactors = FALSE)
# Paths (adjust if necessary)
pca_rdata_path     <- "~/mestrado/pca_result_R30.RData"
t_matrix_rdata_path <- "~/mestrado/t_matrix_R30.RData"

out_images   <- "images"
out_regions  <- "regions_data"
out_tables   <- "enrichment_results"
dir.create(out_images, showWarnings = FALSE, recursive = TRUE)
dir.create(out_regions, showWarnings = FALSE, recursive = TRUE)
dir.create(out_tables, showWarnings = FALSE, recursive = TRUE)

# 1. Packages ------------------------------------------------------------------
# Separate CRAN and Bioconductor packages for more robust installation
cran_pkgs <- c("dplyr","ggplot2","stringr","forcats","RColorBrewer","zoo","tibble","gprofiler2")
bioc_pkgs <- c("biomaRt","org.Hs.eg.db","AnnotationDbi",
               "clusterProfiler","ReactomePA","enrichplot","BiocManager")

install_if_missing <- function(pkgs_cran = NULL, pkgs_bioc = NULL) {
  if (!is.null(pkgs_cran)) {
    for (p in pkgs_cran) {
      if (!requireNamespace(p, quietly = TRUE)) {
        message("Installing CRAN pkg: ", p)
        install.packages(p, repos = "https://cloud.r-project.org")
      }
      suppressPackageStartupMessages(library(p, character.only = TRUE))
    }
  }
  if (!is.null(pkgs_bioc)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    for (p in pkgs_bioc) {
      if (!requireNamespace(p, quietly = TRUE)) {
        message("Installing Bioconductor pkg: ", p)
        BiocManager::install(p, ask = FALSE, update = FALSE)
      }
      suppressPackageStartupMessages(library(p, character.only = TRUE))
    }
  }
}

install_if_missing(cran_pkgs, bioc_pkgs)

# 2. Load PCA and transcriptogram objects --------------------------------------
if (!file.exists(pca_rdata_path)) stop("pca_result_R30.RData not found at: ", pca_rdata_path)
if (!file.exists(t_matrix_rdata_path)) stop("t_matrix_R30.RData not found at: ", t_matrix_rdata_path)
load(pca_rdata_path)       # should load pca_result_R30
load(t_matrix_rdata_path)  # should load t_matrix_R30

# Basic checks
if (!exists("pca_result_R30")) stop("Object 'pca_result_R30' not found in loaded RData.")
if (!exists("t_matrix_R30"))  stop("Object 't_matrix_R30' not found in loaded RData.")
if (is.null(pca_result_R30$pca_result$rotation) && is.null(pca_result_R30[["pca_result"]][["rotation"]])) {
  stop("Cannot find PCA rotation matrix in 'pca_result_R30'. Expecting pca_result_R30$pca_result$rotation")
}

# 3. Prepare data (genes, positions, rotation) --------------------------------
gene_names     <- t_matrix_R30@transcriptogramS2[, 1]
gene_positions <- t_matrix_R30@transcriptogramS2[, 2]
pca_rotation   <- pca_result_R30[["pca_result"]][["rotation"]]
rownames(pca_rotation) <- gene_names

n_pcs <- 6L
if (ncol(pca_rotation) < n_pcs) {
  warning("PCA rotation has fewer than ", n_pcs, " PCs. Using available: ", ncol(pca_rotation))
  n_pcs <- ncol(pca_rotation)
}

variation_signal <- rowSums(abs(pca_rotation[, seq_len(n_pcs), drop = FALSE]))
var_df <- data.frame(
  Gene = gene_names,
  Position = gene_positions,
  Variation = variation_signal,
  stringsAsFactors = FALSE
)

# 4. Smooth signal -------------------------------------------------------------
window_size <- 15L
var_df$VariationSmooth <- zoo::rollmean(var_df$Variation, k = window_size, fill = NA, align = "center")

# 5. Define regions automatically based on VariationSmooth --------------------
# Automatic detection of contiguous above-average regions
threshold <- mean(var_df$VariationSmooth, na.rm = TRUE)
var_df$AboveMean <- var_df$VariationSmooth > threshold
rle_above <- rle(var_df$AboveMean)
ends <- cumsum(rle_above$lengths)
starts <- c(1, head(ends, -1) + 1)
region_tbl <- data.frame(
  start_idx = starts[rle_above$values],
  end_idx   = ends[rle_above$values]
)
region_tbl$start_pos <- var_df$Position[region_tbl$start_idx]
region_tbl$end_pos   <- var_df$Position[region_tbl$end_idx]
min_size <- 20  # minimum region size (adjust if needed)
region_tbl <- region_tbl[region_tbl$end_idx - region_tbl$start_idx + 1 >= min_size, ]

regions_limits <- setNames(
  lapply(seq_len(nrow(region_tbl)), function(i) {
    as.numeric(region_tbl[i, c("start_pos", "end_pos")])
  }),
  paste0("AutoRegion_", seq_len(nrow(region_tbl)))
)

# Validação
if (!is.list(regions_limits) || length(regions_limits) == 0)
  stop("regions_limits must be a non-empty named list of numeric vectors of length 2.")

message("✅ Regiões automáticas detectadas: ", length(regions_limits))
print(head(regions_limits, 5))


message("Regiões automáticas detectadas: ", length(regions_limits))
print(regions_limits)

# You can replace the list above with your custom regions. Ensure names are unique.

# Validate regions_limits format
if (!is.list(regions_limits) || length(regions_limits) == 0) stop("regions_limits must be a non-empty named list of numeric vectors of length 2.")

# Build region data.frames
regions <- list()
for (region_name in names(regions_limits)) {
  limits <- regions_limits[[region_name]]
  if (!is.numeric(limits) || length(limits) != 2) stop("Each region must be numeric vector length 2: ", region_name)
  region_df <- var_df %>%
    dplyr::filter(Position >= limits[1], Position <= limits[2]) %>%
    dplyr::arrange(Position) %>%
    dplyr::mutate(Region = region_name)
  regions[[region_name]] <- region_df
}

# 6. Plot overview with highlighted regions (MODIFIED: highlight by colored curve segments) ----
n_regions <- length(regions)
# choose palette with minimum 3 colors for aesthetics; expand if necessary
pal_base <- RColorBrewer::brewer.pal(min(9, max(3, n_regions)), "Set1")
if (n_regions > length(pal_base)) {
  region_colors <- colorRampPalette(pal_base)(n_regions)
} else {
  region_colors <- pal_base[seq_len(n_regions)]
}

p_main <- ggplot(var_df, aes(x = Position, y = VariationSmooth)) +
  geom_line(color = "gray40", size = 0.6) +                                      # base curve (thin gray)
  geom_hline(yintercept = mean(var_df$VariationSmooth, na.rm = TRUE),
             linetype = "dashed", color = "red") +
  labs(
    title = paste0("High Variation Regions"),
    x = "Gene Position (transcriptogram order)",
    y = "Smoothed Variation (Σ|loading|)"
  ) +
  theme_minimal(base_size = 13)

# Instead of geom_rect, draw thicker colored curve segments for each region.
i <- 1L
for (region_name in names(regions)) {
  region <- regions[[region_name]]
  if (nrow(region) == 0) { i <- i + 1L; next } # skip empty regions
  # overlay a thicker colored line only for the region interval
  p_main <- p_main +
    geom_line(data = region, aes(x = Position, y = VariationSmooth),
              color = region_colors[i], size = 1.2, inherit.aes = FALSE)
    # geom_vline(xintercept = mean(region$Position),
    #            color = region_colors[i], linetype = "dotted", size = 0.7)
  i <- i + 1L
}

fn_overview <- file.path(out_images, "high_variation_regions_overview.png")
ggsave(fn_overview, p_main, width = 14, height = 6, dpi = 600, units = "in")
message("Overview plot saved to: ", fn_overview)

# 7. Export region CSVs (sanitizing names) -------------------------------------
safe_name <- function(x) gsub("[^A-Za-z0-9_\\-]", "_", x)
region_names <- names(regions)
for (i in seq_along(region_names)) {
  rn <- region_names[i]
  region_df <- regions[[rn]]
  fname <- sprintf("high_variation_region_%02d_%s.csv", i, safe_name(rn))
  fpath <- file.path(out_regions, fname)
  write.csv(region_df, fpath, row.names = FALSE)
  message(sprintf("Saved region CSV: %s  (genes: %d)", fpath, nrow(region_df)))
  if (nrow(region_df) == 0) message("  -> region empty.")
}

# 8. Detailed plots per region -------------------------------------------------
plot_region_detail <- function(region_df, color, label) {
  if (nrow(region_df) == 0) {
    p <- ggplot() + geom_blank() +
      labs(title = paste0("Detailed Variation - ", label, " (empty)")) + theme_minimal()
    return(p)
  }
  ggplot(region_df, aes(x = Position, y = Variation)) +
    geom_line(color = "gray40", size = 0.6) +
    geom_line(aes(y = VariationSmooth), color = color, size = 0.9, na.rm = TRUE) +
    labs(title = paste0("Detailed Variation - ", label), x = "Gene Position", y = "Variation (Σ|loading|)") +
    theme_minimal(base_size = 12)
}

region_plots <- Map(function(df, col, nm) plot_region_detail(df, col, nm),
                    regions, region_colors, names(regions))

k <- 1L
for (nm in names(region_plots)) {
  outpng <- file.path(out_images, sprintf("region_detail_%02d_%s.png", k, safe_name(nm)))
  ggsave(outpng, region_plots[[nm]], width = 10, height = 5, dpi = 600, units = "in")
  message("Saved region detail plot: ", outpng)
  k <- k + 1
}

# 9. Enrichment pipeline: detect region files automatically ---------------------
region_files <- list.files(out_regions, pattern = "^high_variation_region_.*\\.csv$", full.names = TRUE)
if (length(region_files) == 0) {
  warning("No region CSV files found in '", out_regions, "'. Skipping enrichment.")
} else {
  # read all region CSVs and build union
  region_dfs <- lapply(region_files, function(f) {
    df <- read.csv(f, stringsAsFactors = FALSE)
    # detect gene column: prefer 'Gene' else first column
    colnames(df) <- trimws(colnames(df))
    gene_col <- if ("Gene" %in% colnames(df)) "Gene" else colnames(df)[1]
    df$Gene <- trimws(as.character(df[[gene_col]]))
    df
  })
  genes_union <- unique(unlist(lapply(region_dfs, function(df) df$Gene)))
  message("Found ", length(region_files), " region files. Union unique IDs: ", length(genes_union))
  
  # basic sanity
  if (length(genes_union) < 3) {
    warning("Very small union (<3 IDs). Enrichment unlikely to work; script will attempt mapping but skip some steps.")
  }
  
  # 10. Detect ID type (heuristic) -----------------------------------------------
  detect_id_type <- function(ids) {
    ids_non_na <- ids[!is.na(ids) & ids != ""]
    if (length(ids_non_na) == 0) return("unknown")
    sample_id <- ids_non_na[1]
    if (grepl("^ENSP", sample_id, ignore.case = TRUE)) return("ensembl_peptide_id")
    if (grepl("^ENSG", sample_id, ignore.case = TRUE)) return("ensembl_gene_id")
    if (grepl("^ENS", sample_id, ignore.case = TRUE))  return("ensembl")
    if (all(grepl("^[0-9]+$", ids_non_na))) return("entrez")
    return("symbol")
  }
  
  id_type <- detect_id_type(genes_union)
  message("Detected ID type (heuristic): ", id_type)
  
  # 11. Map IDs -> Entrez & Symbol via biomaRt (with tryCatch) --------------------
  map_to_entrez_safe <- function(query_ids, id_type) {
    res <- NULL
    try({
      mart <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
      if (id_type == "ensembl_peptide_id") {
        attrs <- c("ensembl_peptide_id","ensembl_gene_id","entrezgene_id","external_gene_name")
        res <- biomaRt::getBM(attributes = attrs, filters = "ensembl_peptide_id",
                              values = unique(query_ids), mart = mart)
        if (nrow(res) > 0) colnames(res)[colnames(res) == "ensembl_peptide_id"] <- "query_id"
      } else if (id_type == "ensembl_gene_id") {
        attrs <- c("ensembl_gene_id","ensembl_peptide_id","entrezgene_id","external_gene_name")
        res <- biomaRt::getBM(attributes = attrs, filters = "ensembl_gene_id",
                              values = unique(query_ids), mart = mart)
        if (nrow(res) > 0) colnames(res)[colnames(res) == "ensembl_gene_id"] <- "query_id"
      } else if (id_type == "entrez") {
        res <- data.frame(query_id = query_ids, entrezgene_id = as.integer(query_ids), external_gene_name = NA_character_, stringsAsFactors = FALSE)
      } else if (id_type == "symbol") {
        attrs <- c("external_gene_name","ensembl_gene_id","ensembl_peptide_id","entrezgene_id")
        res <- biomaRt::getBM(attributes = attrs, filters = "external_gene_name",
                              values = unique(query_ids), mart = mart)
        if (nrow(res) > 0) colnames(res)[colnames(res) == "external_gene_name"] <- "query_id"
      } else {
        # fallback attempt: try mapping as symbol
        attrs <- c("external_gene_name","ensembl_gene_id","ensembl_peptide_id","entrezgene_id")
        res <- biomaRt::getBM(attributes = attrs, filters = "external_gene_name",
                              values = unique(query_ids), mart = mart)
        if (nrow(res) > 0) colnames(res)[colnames(res) == "external_gene_name"] <- "query_id"
      }
    }, silent = TRUE)
    
    # If nothing returned and we have SYMBOLS fallback to org.Hs.eg.db
    if ((is.null(res) || nrow(res) == 0) && id_type %in% c("symbol","unknown")) {
      message("biomaRt returned no results. Trying AnnotationDbi mapIds (org.Hs.eg.db) for SYMBOLs...")
      require(org.Hs.eg.db)
      mapped_entrez <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = unique(query_ids),
                                             column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
      res <- data.frame(query_id = names(mapped_entrez),
                        entrezgene_id = as.integer(mapped_entrez),
                        external_gene_name = names(mapped_entrez),
                        stringsAsFactors = FALSE)
    }
    
    if (is.null(res) || nrow(res) == 0) return(NULL)
    
    # ensure columns exist and normalized names
    if (!"entrezgene_id" %in% colnames(res)) res$entrezgene_id <- NA_integer_
    if (!"external_gene_name" %in% colnames(res)) res$external_gene_name <- NA_character_
    if (!"query_id" %in% colnames(res)) colnames(res)[1] <- "query_id"
    res <- res %>% dplyr::distinct(query_id, .keep_all = TRUE)
    return(res)
  }
  
  mapping <- map_to_entrez_safe(genes_union, id_type)
  if (is.null(mapping) || nrow(mapping) == 0) {
    warning("Mapping returned zero rows. No Entrez/symbols found. Will not run clusterProfiler/Reactome; will attempt g:Profiler only if plausible symbols exist.")
    # Build fallback mapping with query_id only, but DO NOT treat these as gene symbols for enrichment automatically.
    mapping <- data.frame(query_id = genes_union, entrezgene_id = NA_integer_, external_gene_name = NA_character_, stringsAsFactors = FALSE)
  }
  
  # report mapping quality
  mapped_count <- sum(!is.na(mapping$entrezgene_id))
  mapped_symbols_count <- sum(!is.na(mapping$external_gene_name) & mapping$external_gene_name != "")
  message("Mapping summary: total queries = ", length(genes_union),
          "; mapped entrez = ", mapped_count,
          "; mapped symbols = ", mapped_symbols_count)
  
  # if some queries map to multiple rows, warn
  # note: earlier we used distinct(query_id, .keep_all = TRUE) which keeps first occurance.
  # to detect duplicates before distinct we'd have to re-query; warn generically:
  if (nrow(mapping) < length(unique(genes_union))) {
    message("Note: some queries did not map or were dropped due to duplicate handling. Inspect mapping object for details.")
  }
  
  # 12. Prepare vectors for enrichment ------------------------------------------------
  mapping <- mapping %>%
    dplyr::mutate(entrez = as.character(entrezgene_id),
                  entrez = ifelse(is.na(entrez) | entrez == "NA", NA_character_, entrez),
                  symbol = ifelse(is.na(external_gene_name) | external_gene_name == "", NA_character_, external_gene_name))
  
  entrez_ids <- unique(na.omit(mapping$entrez))
  symbol_ids <- unique(na.omit(mapping$symbol))
  
  message("Prepared for enrichment: Entrez IDs = ", length(entrez_ids), "; Symbols = ", length(symbol_ids))
  
  # Save mapping for inspection
  write.csv(mapping, file.path(out_tables, "gene_id_mapping.csv"), row.names = FALSE)
  
  # 13. clusterProfiler / ReactomePA enrichment (prefer Entrez) -------------------------
  save_enrich_to_csv <- function(obj, fname) {
    if (is.null(obj)) return(NULL)
    df <- as.data.frame(obj)
    if (nrow(df) == 0) return(NULL)
    write.csv(df, file = file.path(out_tables, fname), row.names = FALSE)
    return(df)
  }
  
  ego_bp <- ego_mf <- ego_cc <- kegg_res <- reactome_res <- NULL
  if (length(entrez_ids) >= 3) {
    # GO
    try({
      ego_bp <- clusterProfiler::enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                                          ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
      ego_mf <- clusterProfiler::enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                                          ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
      ego_cc <- clusterProfiler::enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                                          ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
      save_enrich_to_csv(ego_bp, "GO_BP_clusterProfiler.csv")
      save_enrich_to_csv(ego_mf, "GO_MF_clusterProfiler.csv")
      save_enrich_to_csv(ego_cc, "GO_CC_clusterProfiler.csv")
    }, silent = TRUE)
    
    # KEGG
    try({
      kegg_res <- clusterProfiler::enrichKEGG(gene = entrez_ids, organism = "hsa", pAdjustMethod = "BH", qvalueCutoff = 0.05)
      if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
        kegg_readable <- clusterProfiler::setReadable(kegg_res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
        save_enrich_to_csv(kegg_readable, "KEGG_clusterProfiler_readable.csv")
      }
    }, silent = TRUE)
    
    # Reactome
    try({
      reactome_res <- ReactomePA::enrichPathway(gene = entrez_ids, organism = "human", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
      save_enrich_to_csv(reactome_res, "ReactomePA_enrichPathway.csv")
    }, silent = TRUE)
  } else {
    message("Too few Entrez IDs (", length(entrez_ids), ") for clusterProfiler/Reactome analyses. Skipping those steps.")
  }
  
  # 14. g:Profiler (symbols) as complementary resource --------------------------------
  gprof_res <- NULL
  if (length(symbol_ids) >= 2) {
    message("Running g:Profiler on mapped gene symbols (n = ", length(symbol_ids), ") ...")
    try({
      gprof_res <- gprofiler2::gost(query = symbol_ids, organism = "hsapiens",
                                    sources = c("GO:BP","GO:MF","GO:CC","KEGG","REAC"), user_threshold = 0.05, correction_method = "fdr")
      if (!is.null(gprof_res) && !is.null(gprof_res$result) && nrow(gprof_res$result) > 0) {
        write.csv(gprof_res$result, file.path(out_tables, "gprofiler2_results.csv"), row.names = FALSE)
        message("g:Profiler results saved.")
      } else {
        message("g:Profiler returned no significant results.")
      }
    }, silent = TRUE)
  } else {
    message("Not enough mapped gene symbols (", length(symbol_ids), ") for g:Profiler. Skipping.")
  }
  
  # 15. Visualization of enrichment results (if any) --------------------------------
  plot_top_terms_bar <- function(enrich_df, title, out_png, topn = 20, label_wrap = 60) {
    if (is.null(enrich_df) || nrow(enrich_df) == 0) { message("No terms to plot for ", title); return(NULL) }
    df <- as.data.frame(enrich_df) %>% dplyr::arrange(p.adjust) %>% head(topn)
    if ("Description" %in% names(df)) df$term_label <- stringr::str_wrap(df$Description, label_wrap)
    else if ("term_name" %in% names(df)) df$term_label <- stringr::str_wrap(df$term_name, label_wrap)
    else df$term_label <- stringr::str_wrap(as.character(df[,1]), label_wrap)
    if (!("p.adjust" %in% names(df))) {
      pcol <- names(df)[grep("p.adjust|p_value|pval|p.value|pvalue|p_val", names(df), ignore.case = TRUE)][1]
      if (!is.na(pcol)) df$p.adjust <- as.numeric(df[[pcol]]) else df$p.adjust <- NA_real_
    }
    df <- df %>% dplyr::mutate(logp = -log10(p.adjust))
    p <- ggplot(df, aes(x = forcats::fct_reorder(term_label, logp), y = logp)) +
      geom_col(fill = "#2C7FB8") + coord_flip() +
      labs(title = title, x = NULL, y = expression(-log[10]~"(adjusted p-value)")) +
      theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold", hjust = 0.5))
    ggsave(out_png, p, width = 9, height = min(7, 0.4 * nrow(df) + 2), dpi = 600)
    return(p)
  }
  
  # Create plots where results exist
  if (!is.null(ego_bp) && nrow(as.data.frame(ego_bp)) > 0) {
    plot_top_terms_bar(as.data.frame(ego_bp), "GO Biological Process (clusterProfiler)", file.path(out_images, "GO_BP_clusterProfiler_bar.png"), topn = 25)
    message("GO BP plot saved.")
  }
  if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
    plot_top_terms_bar(as.data.frame(kegg_res), "KEGG (clusterProfiler)", file.path(out_images, "KEGG_clusterProfiler_bar.png"), topn = 20)
    message("KEGG plot saved.")
  }
  if (!is.null(reactome_res) && nrow(as.data.frame(reactome_res)) > 0) {
    plot_top_terms_bar(as.data.frame(reactome_res), "Reactome (ReactomePA)", file.path(out_images, "Reactome_bar.png"), topn = 20)
    message("Reactome plot saved.")
  }
  if (!is.null(gprof_res) && !is.null(gprof_res$result) && nrow(gprof_res$result) > 0) {
    # prepare a gprofiler-friendly df
    gdf <- gprof_res$result %>% dplyr::arrange(p_value) %>% dplyr::slice_head(n = 30)
    gdf$term_label <- stringr::str_wrap(gdf$term_name, 70)
    gdf <- gdf %>% dplyr::mutate(logp = -log10(p_value))
    p_g <- ggplot(gdf, aes(x = forcats::fct_reorder(term_label, logp), y = logp, fill = source)) +
      geom_col(width = 0.7, alpha = 0.9, color = "gray20") + coord_flip() +
      labs(title = "Functional Enrichment (g:Profiler) — Top Terms", x = NULL, y = expression(-log[10](p~value))) +
      theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "bottom")
    ggsave(file.path(out_images, "gProfiler_top_terms.png"), p_g, width = 10, height = 7, dpi = 600)
    message("g:Profiler top terms plot saved.")
  } else {
    message("No g:Profiler plotting (no results).")
  }
  
  # 16. Hallmark EMT check (msigdbr) - optional but useful
  # Try to fetch Hallmark EMT list via msigdbr if package present
  if (requireNamespace("msigdbr", quietly = TRUE)) {
    msig <- tryCatch({
      msigdbr::msigdbr(species = "Homo sapiens", category = "H")
    }, error = function(e) NULL)
    if (!is.null(msig)) {
      emt_set <- msig %>% dplyr::filter(gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") %>% dplyr::pull(gene_symbol) %>% unique()
      emt_overlap <- intersect(symbol_ids, emt_set)
      message("Hallmark EMT genes overlapped (n = ", length(emt_overlap), "): ", if (length(emt_overlap) > 0) paste(emt_overlap, collapse = ", ") else "(none)")
      write.csv(data.frame(emt_overlap = emt_overlap), file.path(out_tables, "hallmark_emt_overlap.csv"), row.names = FALSE)
    } else {
      message("msigdbr call failed; skipping Hallmark EMT overlap check.")
    }
  } else {
    message("msigdbr not installed; skipping Hallmark EMT overlap check.")
  }
  
  # 17. Save final gene lists and session info ----------------------------------
  write.csv(data.frame(entrez = entrez_ids), file.path(out_tables, "entrez_gene_list.csv"), row.names = FALSE)
  write.csv(data.frame(symbol = symbol_ids), file.path(out_tables, "symbol_gene_list.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(out_tables, "sessionInfo.txt"))
  message("Enrichment pipeline finished. Tables/images saved under: ", normalizePath(out_tables))
}

# End of script
message("All done. Check outputs in directories:\n - ", normalizePath(out_images),
        "\n - ", normalizePath(out_regions),
        "\n - ", normalizePath(out_tables))


# combine_auto_regions.R  — versão corrigida e robusta
# Une todos os arquivos high_variation_region_*_AutoRegion_*.csv com coerção segura de tipos

library(dplyr)
library(readr)
library(stringr)
library(lubridate)
library(purrr)

regions_dir <- "regions_data"
out_fname_prefix <- "high_variation_autoRegions_combined"

# localizar apenas arquivos que contenham "AutoRegion"
files_all <- list.files(regions_dir, pattern = "^high_variation_region_.*AutoRegion.*\\.csv$", full.names = TRUE)
if (length(files_all) == 0) stop("Nenhum arquivo com 'AutoRegion' encontrado.")

# função robusta de conversão numérica
safe_to_numeric <- function(x) {
  x0 <- as.character(x)
  x0[is.na(x0) | x0 == ""] <- NA_character_
  x0 <- str_trim(x0)
  comma_no_dot <- grepl(",", x0) & !grepl("\\.", x0)
  x0[comma_no_dot] <- gsub(",", ".", x0[comma_no_dot])
  x0 <- gsub("\\s+", "", x0)
  x1 <- gsub("[^0-9eE+\\-\\.]", "", x0)
  suppressWarnings(as.numeric(x1))
}

# extrair nome da região
extract_region_name <- function(fname) {
  bn <- basename(fname)
  m <- regmatches(bn, regexpr("AutoRegion[_A-Za-z0-9\\-]*", bn, ignore.case = TRUE))
  if (length(m) == 1 && nzchar(m)) return(m)
  return(sub("\\.csv$", "", bn, ignore.case = TRUE))
}

# função principal de leitura e normalização
process_file <- function(f) {
  df <- tryCatch(read_csv(f, show_col_types = FALSE),
                 error = function(e) { 
                   warning("read_csv falhou em ", basename(f), " - tentando read.csv"); 
                   read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
                 })
  
  colnames(df) <- trimws(colnames(df))
  
  # detectar coluna de gene
  gene_col <- if ("Gene" %in% colnames(df)) "Gene" else if ("gene" %in% colnames(df)) "gene" else colnames(df)[1]
  df$Gene_raw <- as.character(df[[gene_col]])
  
  # garantir colunas essenciais
  if (!"Position" %in% colnames(df)) df$Position <- NA
  if (!"Variation" %in% colnames(df)) df$Variation <- NA
  if (!"VariationSmooth" %in% colnames(df)) df$VariationSmooth <- NA
  
  # conversão numérica segura
  df$Position_num <- safe_to_numeric(df$Position)
  df$Variation_num <- safe_to_numeric(df$Variation)
  df$VariationSmooth_num <- safe_to_numeric(df$VariationSmooth)
  
  region_name <- extract_region_name(f)
  
  df <- df %>%
    mutate(
      Gene = str_trim(as.character(Gene_raw)),
      Position = Position_num,
      Variation = Variation_num,
      VariationSmooth = VariationSmooth_num,
      region_source_file = basename(f),
      region_name = region_name
    ) %>%
    select(Gene, Position, Variation, VariationSmooth, region_source_file, region_name, everything())
  
  # Forçar todas as colunas não numéricas para character
  df[] <- lapply(df, function(x) {
    if (is.logical(x)) as.character(x)
    else if (is.factor(x)) as.character(x)
    else x
  })
  
  return(df)
}

# processar arquivos
processed_list <- lapply(files_all, process_file)

# padronizar colunas antes de unir
all_cols <- unique(unlist(lapply(processed_list, names)))
processed_list <- lapply(processed_list, function(df) {
  missing <- setdiff(all_cols, names(df))
  if (length(missing) > 0) df[missing] <- NA
  df <- df[, all_cols]
  return(df)
})

# unir de forma segura
combined_df <- bind_rows(processed_list)

# salvar com timestamp
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_file <- file.path(regions_dir, paste0(out_fname_prefix, "_", ts, ".csv"))
write_csv(combined_df, out_file)

# gerar resumo
summary_df <- combined_df %>%
  group_by(region_name, region_source_file) %>%
  summarise(rows = n(), unique_genes = n_distinct(Gene), .groups = "drop")

summary_file <- file.path(regions_dir, paste0("combine_summary_", ts, ".csv"))
write_csv(summary_df, summary_file)

message("\n✅ Combinação concluída com sucesso.")
message("Arquivo combinado: ", out_file)
message("Resumo salvo em: ", summary_file)

# ------------------------------------------------------------------
# Salvar agrupamento de genes (clusters) para enriquecimento futuro
# ------------------------------------------------------------------

if (exists("gene_clusters") && nrow(gene_clusters) > 0) {
  message("\nSalvando agrupamento de genes para análise de enriquecimento...")
  
  # Garantir colunas esperadas
  expected_cols <- c("Cluster", "GeneID", "GeneName")
  missing_cols <- setdiff(expected_cols, colnames(gene_clusters))
  if (length(missing_cols) > 0) {
    warning("As seguintes colunas estavam ausentes e serão criadas vazias: ", paste(missing_cols, collapse = ", "))
    for (col in missing_cols) gene_clusters[[col]] <- NA
  }
  
  # Ordenar para consistência
  gene_clusters <- gene_clusters %>%
    dplyr::arrange(Cluster, GeneID)
  
  # Caminho de saída
  cluster_csv <- file.path(output_dir, "clustered_genes_for_enrichment.csv")
  
  # Exportar CSV
  readr::write_csv(gene_clusters, cluster_csv)
  message("Arquivo salvo com sucesso: ", cluster_csv)
  
} else {
  warning("Objeto 'gene_clusters' não encontrado ou vazio. Nenhum arquivo de agrupamento foi salvo.")
}


