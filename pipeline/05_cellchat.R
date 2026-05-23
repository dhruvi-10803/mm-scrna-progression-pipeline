# 05_cellchat.R
# Pipeline step 5: Cell-cell communication inference using CellChat.
# Runs CellChat separately for T1 and T2 per patient using az_l2 cell type labels.
# Computes pathway-level communication probabilities and interaction counts.
# Saves T1 and T2 CellChat objects, pathway diff CSV, all interactions CSV,
# and cell type interaction change summary across all patients.
# Input:  Per-patient scored objects from consolidation/module_scoring/Px_scored.rds
# Output: CellChat objects in TME/objects/Px/
#         Results in TME/results/Px_cellchat_results/
#         Cross-patient summary in TME/results/celltype_interaction_changes.csv

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(ggplot2)
  library(dplyr)
  library(data.table)
})

options(future.globals.maxSize = 2000 * 1024^2)
options(CellChat.verbose = FALSE)

obj_dir     <- "~/Dhruvi/consolidation/module_scoring/"
out_obj_dir <- "~/Dhruvi/TME/objects/"
out_res_dir <- "~/Dhruvi/TME/results/"

dir.create(path.expand(out_obj_dir), recursive = TRUE,
           showWarnings = FALSE)
dir.create(path.expand(out_res_dir), recursive = TRUE,
           showWarnings = FALSE)

patients  <- paste0("P", 1:7)
MIN_CELLS <- 10

# Run CellChat on a Seurat object subset.
# Groups cells by az_l2 label, removes groups below MIN_CELLS threshold,
# computes communication probabilities using truncated mean approach.
run_cellchat <- function(seurat_obj, min_cells = MIN_CELLS) {

  cell_counts <- table(seurat_obj@meta.data$az_l2)
  keep_groups <- names(cell_counts)[cell_counts >= min_cells]

  if (length(keep_groups) < 2) {
    cat("  Too few cell groups, skipping\n")
    return(NULL)
  }

  seurat_obj <- subset(
    seurat_obj,
    cells = colnames(seurat_obj)[
      seurat_obj@meta.data$az_l2 %in% keep_groups
    ]
  )

  cat(sprintf("  Cell groups: %d | Cells: %d\n",
              length(keep_groups), ncol(seurat_obj)))

  data_input <- GetAssayData(seurat_obj, assay = "RNA",
                              layer = "data")
  meta_input <- seurat_obj@meta.data[, "az_l2", drop = FALSE]
  colnames(meta_input) <- "labels"

  cc <- createCellChat(object  = data_input,
                        meta    = meta_input,
                        group.by = "labels")

  cc@DB <- CellChatDB.human

  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)

  # Truncated mean reduces influence of outlier-expressing cells
  cc <- computeCommunProb(cc, type = "triMean",
                           nboot = 100,
                           population.size = TRUE)
  cc <- filterCommunication(cc, min.cells = min_cells)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)

  return(cc)
}

# Extract all significant interactions from a CellChat object
# into a flat data frame with source, target, ligand, receptor,
# probability, pval, pathway, annotation, and timepoint columns
extract_interactions <- function(cc, timepoint_label) {

  df <- subsetCommunication(cc)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df$time <- timepoint_label
  df
}

# Compute pathway-level communication probability at T1 and T2
# and log2FC between timepoints for each pathway
compute_pathway_diff <- function(cc_t1, cc_t2) {

  pw_t1 <- cc_t1@netP$prob
  pw_t2 <- cc_t2@netP$prob

  all_pathways <- union(names(pw_t1), names(pw_t2))

  results <- lapply(all_pathways, function(pw) {
    s1 <- if (pw %in% names(pw_t1)) sum(pw_t1[[pw]]) else 0
    s2 <- if (pw %in% names(pw_t2)) sum(pw_t2[[pw]]) else 0
    lfc <- log2((s2 + 1e-10) / (s1 + 1e-10))
    data.frame(pathway = pw, T1 = s1, T2 = s2,
               log2FC = lfc,
               stringsAsFactors = FALSE)
  })

  df <- do.call(rbind, results)
  df[order(-abs(df$log2FC)), ]
}

# Compute total interactions per cell type at T1 and T2
# and classify as Emerging, Disappearing, Increasing, or Decreasing
compute_celltype_changes <- function(cc_t1, cc_t2, patient) {

  mat_t1 <- cc_t1@net$count
  mat_t2 <- cc_t2@net$count

  all_cells <- union(rownames(mat_t1), rownames(mat_t2))

  results <- lapply(all_cells, function(ct) {

    t1_send <- if (ct %in% rownames(mat_t1))
      sum(mat_t1[ct, ]) else 0
    t1_recv <- if (ct %in% colnames(mat_t1))
      sum(mat_t1[, ct]) else 0
    t1_total <- t1_send + t1_recv

    t2_send <- if (ct %in% rownames(mat_t2))
      sum(mat_t2[ct, ]) else 0
    t2_recv <- if (ct %in% colnames(mat_t2))
      sum(mat_t2[, ct]) else 0
    t2_total <- t2_send + t2_recv

    delta <- t2_total - t1_total

    status <- case_when(
      t1_total == 0 & t2_total > 0 ~ "Emerging",
      t1_total > 0 & t2_total == 0 ~ "Disappearing",
      delta > 0                     ~ "Increasing",
      delta < 0                     ~ "Decreasing",
      TRUE                          ~ "Stable"
    )

    data.frame(
      patient   = patient,
      celltype  = ct,
      T1_total  = t1_total,
      T2_total  = t2_total,
      delta     = delta,
      abs_delta = abs(delta),
      status    = status,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, results)
}

all_celltype_changes <- list()

for (pt in patients) {

  cat("\nPatient:", pt, "\n")

  obj_path <- file.path(
    path.expand(obj_dir),
    paste0(pt, "_scored.rds")
  )

  if (!file.exists(obj_path)) {
    cat("  Scored object not found, skipping\n")
    next
  }

  px <- readRDS(obj_path)
  cat(sprintf("  Loaded: %d cells\n", ncol(px)))

  DefaultAssay(px) <- "RNA"
  tryCatch(px <- JoinLayers(px), error = function(e) NULL)
  px <- NormalizeData(px, verbose = FALSE)

  # Create output directories for this patient
  pt_obj_dir <- file.path(path.expand(out_obj_dir), pt)
  pt_res_dir <- file.path(
    path.expand(out_res_dir),
    paste0(pt, "_cellchat_results")
  )
  dir.create(pt_obj_dir, recursive = TRUE,
             showWarnings = FALSE)
  dir.create(pt_res_dir, recursive = TRUE,
             showWarnings = FALSE)

  # Subset T1 and T2
  px_t1 <- subset(px, subset = timepoint == "TimePoint_1")
  px_t2 <- subset(px, subset = timepoint == "TimePoint_2")

  cat("  T1 cells:", ncol(px_t1),
      "| T2 cells:", ncol(px_t2), "\n")

  # Run CellChat on T1
  cat("  Running CellChat T1\n")
  cc_t1 <- tryCatch(
    run_cellchat(px_t1),
    error = function(e) {
      cat("  T1 error:", e$message, "\n")
      NULL
    }
  )

  # Run CellChat on T2
  cat("  Running CellChat T2\n")
  cc_t2 <- tryCatch(
    run_cellchat(px_t2),
    error = function(e) {
      cat("  T2 error:", e$message, "\n")
      NULL
    }
  )

  if (is.null(cc_t1) || is.null(cc_t2)) {
    cat("  CellChat failed for", pt, "\n")
    rm(px, px_t1, px_t2)
    gc()
    next
  }

  # Save T1 and T2 CellChat objects
  saveRDS(cc_t1,
          file.path(pt_obj_dir,
                    paste0(pt, "_cellchat_T1.rds")))
  saveRDS(cc_t2,
          file.path(pt_obj_dir,
                    paste0(pt, "_cellchat_T2.rds")))
  cat("  CellChat objects saved\n")

  # Compute and save pathway-level diff
  pathway_diff <- tryCatch(
    compute_pathway_diff(cc_t1, cc_t2),
    error = function(e) {
      cat("  Pathway diff error:", e$message, "\n")
      NULL
    }
  )

  if (!is.null(pathway_diff)) {
    write.csv(
      pathway_diff,
      file.path(pt_res_dir,
                paste0(pt, "_pathway_diff.csv")),
      row.names = FALSE
    )
    cat(sprintf("  Pathway diff saved: %d pathways\n",
                nrow(pathway_diff)))
  }

  # Extract and save all significant interactions
  int_t1 <- tryCatch(
    extract_interactions(cc_t1, "T1"),
    error = function(e) NULL
  )
  int_t2 <- tryCatch(
    extract_interactions(cc_t2, "T2"),
    error = function(e) NULL
  )

  all_int <- rbind(int_t1, int_t2)

  if (!is.null(all_int) && nrow(all_int) > 0) {
    write.csv(
      all_int,
      file.path(pt_res_dir,
                paste0(pt, "_all_interactions.csv")),
      row.names = FALSE
    )
    cat(sprintf("  All interactions saved: %d rows\n",
                nrow(all_int)))
  }

  # Compute cell type interaction changes T1 vs T2
  ct_changes <- tryCatch(
    compute_celltype_changes(cc_t1, cc_t2, pt),
    error = function(e) {
      cat("  Celltype changes error:", e$message, "\n")
      NULL
    }
  )

  if (!is.null(ct_changes)) {
    all_celltype_changes[[pt]] <- ct_changes
    cat(sprintf("  Cell type changes computed: %d types\n",
                nrow(ct_changes)))
  }

  rm(px, px_t1, px_t2, cc_t1, cc_t2)
  gc()

  cat(pt, "complete\n")
}

# Save cross-patient cell type interaction change summary
if (length(all_celltype_changes) > 0) {

  combined_changes <- do.call(rbind, all_celltype_changes)

  write.csv(
    combined_changes,
    file.path(path.expand(out_res_dir),
              "celltype_interaction_changes.csv"),
    row.names = FALSE
  )

  cat("\nCross-patient cell type changes saved\n")
  print(
    combined_changes %>%
      group_by(patient, status) %>%
      summarise(n = n(), .groups = "drop") %>%
      tidyr::pivot_wider(names_from  = status,
                         values_from = n,
                         values_fill = 0)
  )
}

cat("\nCellChat complete\n")
cat(format(Sys.time()), "\n")