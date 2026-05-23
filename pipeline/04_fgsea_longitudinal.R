# 03_fgsea_longitudinal.R
# Pipeline step 3: Gene set enrichment analysis on longitudinal DEG results.
# Genes ranked by avg_log2FC from MAST output. Gene sets from Hallmark,
# KEGG, and Reactome collections. Run per patient per cell type.
# Also identifies pathways showing opposite directionality between
# stable and progressor patients across cell types.
# Input:  MAST results from consolidation/mast_az_l2/
# Output: Per-patient per-celltype fgsea results in consolidation/fgsea_longitudinal/

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(msigdbr)
})

mast_dir <- "~/Dhruvi/consolidation/mast_az_l2/"
out_dir  <- "~/Dhruvi/consolidation/fgsea_longitudinal/"
dir.create(path.expand(out_dir), recursive = TRUE,
           showWarnings = FALSE)

patients <- paste0("P", 1:7)

stable_pts    <- c("P1", "P2")
progressor_pts <- c("P3", "P4", "P5", "P6", "P7")

# Same noise gene filter as DEG step for consistency
NOISE <- paste0(
  "^IG[HKL]|^IGLC|^IGKC|^IGKV|^IGLV|^IGHV|",
  "^RPS|^RPL|^MT-|^MTRNR|MALAT1|NEAT1|",
  "^AC[0-9]|^AL[0-9]|^LINC|^SNHG|",
  "HBA1|HBA2|HBB|HBD|HBE1|HBG1|HBG2|ALAS2|AHSP|GYPA"
)

# Load gene sets from MSigDB: Hallmark, KEGG, and Reactome
cat("Loading gene sets from MSigDB\n")

get_sets <- function(category, subcategory = NULL) {
  args <- list(species = "Homo sapiens", category = category)
  if (!is.null(subcategory)) args$subcategory <- subcategory
  dt <- do.call(msigdbr, args)
  split(dt$gene_symbol, dt$gs_name)
}

all_sets <- c(
  get_sets("H"),
  get_sets("C2", "CP:KEGG_LEGACY"),
  get_sets("C2", "CP:REACTOME")
)

cat(sprintf("Gene sets loaded: %d\n", length(all_sets)))

# Load all MAST results and remove noise genes
cat("Loading MAST results\n")
mast_files <- list.files(path.expand(mast_dir),
                          pattern = "_MAST\\.csv$",
                          full.names = TRUE)

all_mast <- rbindlist(lapply(mast_files, fread), fill = TRUE)
all_mast <- all_mast[
  !grepl(NOISE, gene) &
    !is.na(avg_log2FC) &
    !is.na(pval)
]

cat(sprintf("MAST rows loaded: %d\n", nrow(all_mast)))

# fgsea function: ranks genes by avg_log2FC, deduplicates,
# then runs fgsea with 1000 permutations
run_fgsea <- function(mast_sub, gene_sets,
                       min_size = 10, max_size = 500) {

  # Deduplicate keeping highest absolute logFC per gene
  mast_sub <- mast_sub[order(-abs(avg_log2FC))]
  mast_sub <- mast_sub[!duplicated(gene)]

  ranks <- setNames(mast_sub$avg_log2FC, mast_sub$gene)
  ranks <- sort(ranks, decreasing = TRUE)

  if (length(ranks) < 50) return(NULL)

  res <- tryCatch(
    fgsea(
      pathways    = gene_sets,
      stats       = ranks,
      minSize     = min_size,
      maxSize     = max_size,
      nPermSimple = 1000,
      eps         = 0
    ),
    error = function(e) NULL
  )

  if (is.null(res)) return(NULL)

  as.data.table(res)[order(pval)]
}

all_results  <- list()
summary_rows <- list()

for (pt in patients) {

  cat("\nPatient:", pt, "\n")
  pt_mast    <- all_mast[patient == pt]
  cell_types <- unique(pt_mast$az_l2)

  for (ct in cell_types) {

    ct_mast <- pt_mast[az_l2 == ct]
    if (nrow(ct_mast) < 50) next

    cat(sprintf("  %s %d genes\n", ct, nrow(ct_mast)))

    out_file <- file.path(
      path.expand(out_dir),
      sprintf("%s_%s_fgsea.csv", pt,
              gsub("[^A-Za-z0-9]", "_", ct))
    )

    if (file.exists(out_file)) {
      cat("  Already exists, skipping\n")
      next
    }

    res <- run_fgsea(ct_mast, all_sets)
    if (is.null(res) || nrow(res) == 0) next

    res[, patient     := pt]
    res[, az_l2       := ct]
    res[, leadingEdge := sapply(leadingEdge,
                                 paste, collapse = ",")]

    fwrite(res, out_file)

    n_up <- sum(res$padj < 0.05 & res$NES > 0,
                na.rm = TRUE)
    n_dn <- sum(res$padj < 0.05 & res$NES < 0,
                na.rm = TRUE)

    top_up <- res[NES > 0][order(pval)][
      1:min(3, .N), pathway]
    top_dn <- res[NES < 0][order(pval)][
      1:min(3, .N), pathway]

    cat(sprintf("  sig up=%d dn=%d\n", n_up, n_dn))

    all_results[[length(all_results) + 1]]  <- res
    summary_rows[[length(summary_rows) + 1]] <- data.table(
      patient  = pt,
      az_l2    = ct,
      n_genes  = nrow(ct_mast),
      n_sig_up = n_up,
      n_sig_dn = n_dn,
      top3_up  = paste(top_up, collapse = " | "),
      top3_dn  = paste(top_dn, collapse = " | ")
    )
  }
}

# Save combined results and per-combination summary
if (length(all_results) > 0) {

  combined <- rbindlist(all_results, fill = TRUE)

  fwrite(combined,
         file.path(path.expand(out_dir),
                   "fgsea_longitudinal_all.csv"))

  summary_dt <- rbindlist(summary_rows)
  fwrite(summary_dt,
         file.path(path.expand(out_dir),
                   "fgsea_longitudinal_summary.csv"))

  cat(sprintf("\nTotal fgsea results: %d rows across %d combinations\n",
              nrow(combined), nrow(summary_dt)))
}

# Cross-patient directionality analysis:
# identifies pathways showing opposite NES direction
# between stable (P1 P2) and progressor (P3-P7) patients
# within the same cell type, using relaxed padj threshold 0.25
cat("\nRunning cross-patient directionality analysis\n")

if (length(all_results) > 0) {

  combined    <- rbindlist(all_results, fill = TRUE)
  sig         <- combined[padj < 0.25]
  dir_results <- list()

  for (ct in unique(sig$az_l2)) {

    ct_sig       <- sig[az_l2 == ct]
    pathways_ct  <- unique(ct_sig$pathway)

    for (pw in pathways_ct) {

      pw_data  <- combined[az_l2 == ct & pathway == pw]
      stab_NES <- pw_data[patient %in% stable_pts, NES]
      prog_NES <- pw_data[patient %in% progressor_pts, NES]

      if (length(stab_NES) == 0 ||
          length(prog_NES) < 2) next

      n_stab_up <- sum(stab_NES > 0, na.rm = TRUE)
      n_stab_dn <- sum(stab_NES < 0, na.rm = TRUE)
      n_prog_up <- sum(prog_NES > 0, na.rm = TRUE)
      n_prog_dn <- sum(prog_NES < 0, na.rm = TRUE)

      # Retain only pathways where stable and progressor
      # patients show clearly opposite directions
      opp_A <- (n_stab_up >= 1 && n_prog_dn >= 3)
      opp_B <- (n_stab_dn >= 1 && n_prog_up >= 3)
      if (!opp_A && !opp_B) next

      dir_results[[length(dir_results) + 1]] <- data.table(
        az_l2      = ct,
        pathway    = pw,
        n_stab     = length(stab_NES),
        n_prog     = length(prog_NES),
        n_stab_up  = n_stab_up,
        n_stab_dn  = n_stab_dn,
        n_prog_up  = n_prog_up,
        n_prog_dn  = n_prog_dn,
        mean_stab  = round(mean(stab_NES, na.rm = TRUE), 3),
        mean_prog  = round(mean(prog_NES, na.rm = TRUE), 3),
        NES_diff   = round(
          mean(prog_NES, na.rm = TRUE) -
            mean(stab_NES, na.rm = TRUE), 3),
        pattern    = ifelse(opp_A,
                             "StabUP_ProgDN",
                             "StabDN_ProgUP"),
        stab_NES   = paste(round(stab_NES, 2),
                            collapse = ","),
        prog_NES   = paste(round(prog_NES, 2),
                            collapse = ",")
      )
    }
  }

  if (length(dir_results) > 0) {

    dir_dt <- rbindlist(dir_results)

    fwrite(dir_dt,
           file.path(path.expand(out_dir),
                     "fgsea_pathway_directionality.csv"))

    cat(sprintf("Pathways with opposite directionality: %d\n",
                nrow(dir_dt)))

    # Strictest filter: all stable patients vs 4+ of 5 progressors
    strict_A <- dir_dt[
      pattern == "StabUP_ProgDN" &
        n_stab_up == 2 & n_prog_dn >= 4
    ][order(-abs(NES_diff))]

    strict_B <- dir_dt[
      pattern == "StabDN_ProgUP" &
        n_stab_dn == 2 & n_prog_up >= 4
    ][order(-abs(NES_diff))]

    # Count how many cell types each pathway appears in
    all_strict <- rbind(strict_A, strict_B, fill = TRUE)

    pw_consistency <- all_strict[, .(
      n_cell_types = .N,
      cell_types   = paste(sort(az_l2), collapse = ", "),
      pattern      = paste(sort(unique(pattern)),
                            collapse = ", ")
    ), by = pathway][order(-n_cell_types)]

    fwrite(pw_consistency,
           file.path(path.expand(out_dir),
                     "fgsea_crosscelltype_pathways.csv"))

    cat("Pathways consistent in 3+ cell types:\n")
    print(pw_consistency[n_cell_types >= 3])
  }
}

cat("\nfgsea longitudinal complete\n")
cat(format(Sys.time()), "\n")