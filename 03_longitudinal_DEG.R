# 02_longitudinal_DEG.R
# Pipeline step 2: Differential gene expression between T1 and T2
# for each patient and each cell type using MAST hurdle model.
# Cells are grouped by az_l2 annotation. Cell types with fewer
# than 15 cells at either timepoint are skipped.
# Input:  Per-patient scored Seurat objects from objects/Px/Px_scored.rds
# Output: Per-patient per-celltype MAST results in consolidation/mast_az_l2/

suppressPackageStartupMessages({
  library(Seurat)
  library(MAST)
  library(data.table)
  library(Matrix)
})

obj_dir <- "~/Dhruvi/consolidation/module_scoring/"
out_dir <- "~/Dhruvi/consolidation/mast_az_l2/"
dir.create(path.expand(out_dir), recursive = TRUE,
           showWarnings = FALSE)

patients  <- paste0("P", 1:7)
MIN_CELLS <- 15

# Noise genes excluded from testing: immunoglobulins, ribosomes,
# mitochondrial, erythroid, and non-coding RNAs
NOISE <- paste0(
  "^IG[HKL]|^IGLC|^IGKC|^IGKV|^IGLV|^IGHV|",
  "^RPS|^RPL|^MT-|^MTRNR|",
  "MALAT1|NEAT1|ALAS2|AHSP|GYPA|",
  "HBA1|HBA2|HBB|HBD|HBE1|HBG1|HBG2|",
  "^AC[0-9]|^AL[0-9]|^LINC|^SNHG"
)

# Cell types to pool across multiple az_l2 labels for sufficient cell numbers
pooled_groups <- list(
  HSPC_pooled      = c("HSC","LMPP","GMP","EMP","CLP","Prog Mk"),
  DC_pooled        = c("cDC1","cDC2","pDC","pre-mDC","pre-pDC"),
  Erythroid_pooled = c("Early Eryth","Late Eryth")
)

# Core MAST function: runs hurdle model comparing T2 vs T1
# Returns data frame with log2FC, p-value, adjusted p-value per gene
run_mast <- function(counts_t1, counts_t2, label, pt) {

  n1 <- ncol(counts_t1)
  n2 <- ncol(counts_t2)
  cat(sprintf("  T1=%d T2=%d cells\n", n1, n2))

  if (n1 < MIN_CELLS || n2 < MIN_CELLS) {
    cat(sprintf("  Skipping: need >= %d cells per timepoint\n",
                MIN_CELLS))
    return(NULL)
  }

  counts_all <- cbind(counts_t1, counts_t2)

  # Keep genes expressed in at least 5% of cells in either timepoint
  pct1 <- rowMeans(counts_t1 > 0)
  pct2 <- rowMeans(counts_t2 > 0)
  keep <- ((pct1 >= 0.05) | (pct2 >= 0.05)) &
    !grepl(NOISE, rownames(counts_all))
  mat  <- counts_all[keep, ]

  if (sum(keep) < 10) {
    cat("  Skipping: too few genes pass filter\n")
    return(NULL)
  }
  cat(sprintf("  Genes tested: %d\n", sum(keep)))

  # Library size normalisation to 10k counts then log1p transform
  lib  <- Matrix::colSums(mat)
  lib[lib == 0] <- 1
  norm <- Matrix::t(Matrix::t(mat) / lib) * 1e4
  lmat <- as.matrix(log1p(norm))

  # Cellular detection rate used as covariate to control for
  # differences in library complexity between cells
  cdr  <- scale(colMeans(lmat > 0))[, 1]

  cdat <- data.frame(
    wellKey   = colnames(lmat),
    timepoint = factor(
      c(rep("T1", n1), rep("T2", n2)),
      levels = c("T1", "T2")
    ),
    ngeneson  = cdr,
    row.names = colnames(lmat),
    stringsAsFactors = FALSE
  )

  fdat <- data.frame(
    primerid  = rownames(lmat),
    row.names = rownames(lmat)
  )

  sca <- tryCatch(
    FromMatrix(exprsArray = lmat, cData = cdat,
               fData = fdat, check_sanity = FALSE),
    error = function(e) NULL
  )
  if (is.null(sca)) return(NULL)

  zlm_fit <- tryCatch(
    suppressMessages(
      zlm(~timepoint + ngeneson, sca, silent = TRUE)
    ),
    error = function(e) NULL
  )
  if (is.null(zlm_fit)) return(NULL)

  summ <- tryCatch(
    suppressMessages(
      summary(zlm_fit,
              doLRT    = "timepointT2",
              logFC    = TRUE)$datatable
    ),
    error = function(e) NULL
  )
  if (is.null(summ)) return(NULL)

  hurdle <- summ[summ$component == "H",
                  c("primerid", "Pr(>Chisq)")]
  lfc_dt <- summ[summ$component == "logFC",
                  c("primerid", "coef")]
  colnames(hurdle) <- c("gene", "pval")
  colnames(lfc_dt) <- c("gene", "avg_log2FC")

  res          <- merge(hurdle, lfc_dt, by = "gene")
  res$padj     <- p.adjust(res$pval, method = "BH")
  res$pct_T1   <- rowMeans(counts_t1[res$gene, , drop = FALSE] > 0)
  res$pct_T2   <- rowMeans(counts_t2[res$gene, , drop = FALSE] > 0)
  res$n_T1     <- n1
  res$n_T2     <- n2
  res$patient  <- pt
  res$az_l2    <- label

  res[order(res$pval), ]
}

summary_all <- list()

cat("Starting longitudinal DEG analysis\n")
cat(format(Sys.time()), "\n\n")

for (pt in patients) {

  cat("Patient:", pt, "\n")

  obj_path <- file.path(
    path.expand(obj_dir),
    sprintf("%s_scored.rds", pt)
  )

  obj <- tryCatch(readRDS(obj_path),
                  error = function(e) {
                    cat("  Cannot load object for", pt, "\n")
                    NULL
                  })

  if (is.null(obj)) next

  DefaultAssay(obj) <- "RNA"
  obj <- JoinLayers(obj)

  counts_full <- tryCatch(
    GetAssayData(obj, assay = "RNA", layer = "counts"),
    error = function(e) NULL
  )

  if (is.null(counts_full) || sum(counts_full) == 0) {
    cat("  Count matrix empty for", pt, "\n")
    next
  }

  meta <- obj@meta.data

  # Get all individual az_l2 cell types present in this patient
  individual_types <- setdiff(
    unique(meta$az_l2),
    unlist(pooled_groups)
  )

  # Build full list of tasks: individual types plus pooled groups
  tasks <- lapply(individual_types, function(ct) {
    list(label = ct, az = ct)
  })

  for (group_name in names(pooled_groups)) {
    tasks[[length(tasks) + 1]] <- list(
      label = group_name,
      az    = pooled_groups[[group_name]]
    )
  }

  for (task in tasks) {

    label <- task$label
    az    <- task$az

    cat(sprintf("  [%s]\n", label))

    out_file <- file.path(
      path.expand(out_dir),
      sprintf("%s_%s_MAST.csv", pt, label)
    )

    if (file.exists(out_file)) {
      cat("  Already exists, skipping\n")
      next
    }

    cells_t1 <- rownames(meta)[
      meta$az_l2 %in% az &
        meta$timepoint == "TimePoint_1"
    ]
    cells_t2 <- rownames(meta)[
      meta$az_l2 %in% az &
        meta$timepoint == "TimePoint_2"
    ]

    ct1 <- counts_full[,
                        intersect(cells_t1,
                                  colnames(counts_full)),
                        drop = FALSE]
    ct2 <- counts_full[,
                        intersect(cells_t2,
                                  colnames(counts_full)),
                        drop = FALSE]

    result <- tryCatch(
      run_mast(ct1, ct2, label, pt),
      error = function(e) {
        cat(sprintf("  Error: %s\n", e$message))
        NULL
      }
    )

    if (!is.null(result) && nrow(result) > 0) {

      fwrite(as.data.table(result), out_file)

      n_up <- sum(result$padj < 0.05 &
                    result$avg_log2FC > 0, na.rm = TRUE)
      n_dn <- sum(result$padj < 0.05 &
                    result$avg_log2FC < 0, na.rm = TRUE)

      cat(sprintf("  Saved: %d up %d down (padj<0.05)\n",
                  n_up, n_dn))

      summary_all[[length(summary_all) + 1]] <- data.table(
        patient    = pt,
        cell_type  = label,
        n_cells_T1 = ncol(ct1),
        n_cells_T2 = ncol(ct2),
        n_genes    = nrow(result),
        n_up_sig   = n_up,
        n_dn_sig   = n_dn
      )
    }
  }

  rm(obj, counts_full, meta)
  gc()

  cat(pt, "complete\n\n")
}

# Save summary of all MAST runs across patients and cell types
if (length(summary_all) > 0) {
  summary_dt <- rbindlist(summary_all)
  fwrite(
    summary_dt,
    file.path(path.expand(out_dir), "MAST_summary_all.csv")
  )
  cat("Summary saved\n")
  print(summary_dt)
}

cat("\nLongitudinal DEG complete\n")
cat(format(Sys.time()), "\n")