# 06_infercnv_clonal_evolution.R
# Pipeline step 6: Copy number inference and clonal evolution analysis.
# This script has four sequential sections that must be run in order.
#
# Section A: QC and doublet removal for all samples including HD donors
# Section B: Per-patient clustering with Harmony batch correction
#            and manual plasma cell cluster identification
# Section C: HD reference annotation and inferCNV run per patient
# Section D: Clonal evolution classification and Sankey diagram generation
#
# Input:  Raw H5 files in CNV_analysis/data/
#         Gene order file in CNV_analysis/references/gene_order.txt
# Output: Clonal evolution calls in CNV_analysis/results/clonal_evolution/
#         Sankey diagrams in CNV_analysis/results/CNV_plots/
#         inferCNV objects in CNV_analysis/results/inferCNV/

suppressPackageStartupMessages({
  library(Seurat)
  library(infercnv)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(harmony)
  library(glmGamPoi)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
  library(ggrepel)
  library(ape)
  library(phangorn)
  library(ComplexHeatmap)
  library(circlize)
  library(ggalluvial)
  library(future)
})

setwd("~/Dhruvi/CNV_analysis/scripts")
options(future.globals.maxSize = 10 * 1024^3)

patients  <- paste0("P", 1:7)
stage_map <- c(
  P1 = "MGUS-MGUS", P2 = "MGUS-MGUS",
  P3 = "MGUS-SMM",  P4 = "MGUS-MM",
  P5 = "SMM-MM",    P6 = "SMM-MM",   P7 = "SMM-MM"
)
progressor <- c(
  P1 = FALSE, P2 = FALSE,
  P3 = TRUE,  P4 = TRUE,
  P5 = TRUE,  P6 = TRUE, P7 = TRUE
)

gene_order_file <- "../references/gene_order.txt"
if (!file.exists(gene_order_file))
  stop("Gene order file missing. Run EnsDb extraction first.")

chr_order <- paste0("chr", c(1:22, "X"))

# Sample manifest for all 19 samples
sample_map <- data.frame(
  gsm = c(
    "GSM8369863_HD1","GSM8369864_HD2","GSM8369865_HD3",
    "GSM8369866_HD4","GSM8369867_HD5",
    "GSM8369868_MGUS_1","GSM8369869_MGUS_2",
    "GSM8369870_MGUS_3","GSM8369871_MGUS_4",
    "GSM8369872_MGUS_5","GSM8369873_MGUS_6",
    "GSM8369874_SMM_1","GSM8369875_SMM_2",
    "GSM8369876_SMM_3","GSM8369877_SMM_4",
    "GSM8369878_MM_1","GSM8369879_MM_2",
    "GSM8369880_MM_3","GSM8369881_MM_4"
  ),
  stage = c(
    rep("Healthy", 5), rep("MGUS", 6),
    rep("SMM", 4), rep("MM", 4)
  ),
  patient = c(
    "HD1","HD2","HD3","HD4","HD5",
    "P1","P1","P2","P2","P3","P4",
    "P3","P5","P6","P7",
    "P4","P5","P6","P7"
  ),
  timepoint = c(
    rep("HD", 5),
    "T1","T2","T1","T2","T1","T1",
    "T2","T1","T1","T1",
    "T2","T2","T2","T2"
  ),
  stringsAsFactors = FALSE
)

# Plasma cell cluster IDs identified manually per patient after
# inspecting dotplots and marker CSVs from Section B
# These must be updated if Section B is rerun with different parameters
pc_clusters_final <- list(
  P1 = c("23","24"),
  P2 = c("15","16","22"),
  P3 = c("9","18","20","21","23"),
  P4 = c("10","13"),
  P5 = c("3","12","23"),
  P6 = c("13","24"),
  P7 = c("7","21","24")
)

# ============================================================
# SECTION A: QC and doublet removal for all samples
# ============================================================

cat("Section A: QC and doublet removal\n")

dirs <- c(
  "../objects/clean_samples",
  "../results/QC/per_sample_plots",
  "../results/QC/tables"
)
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Adaptive QC thresholds based on sample-specific distribution
# Lower bound: median minus 2 MAD, floored at minimum
# Upper bound: median plus 3 MAD
# No mitochondrial cutoff applied — plasma cells have elevated MT by biology
adaptive_lower <- function(x, floor_val) {
  max(median(x, na.rm = TRUE) - 2 * mad(x, na.rm = TRUE), floor_val)
}
adaptive_upper <- function(x) {
  median(x, na.rm = TRUE) + 3 * mad(x, na.rm = TRUE)
}

summary_df <- data.frame()
data_dir   <- "../data"
h5_files   <- list.files(data_dir, pattern = "\\.h5$",
                           full.names = TRUE)

for (f in h5_files) {

  sample_name <- gsub("\\.h5$", "", basename(f))
  clean_path  <- paste0("../objects/clean_samples/",
                         sample_name, "_clean.rds")

  if (file.exists(clean_path)) {
    cat("Already clean:", sample_name, "\n")
    next
  }

  cat("QC:", sample_name, "\n")

  mat <- Read10X_h5(f)
  obj <- CreateSeuratObject(counts = mat, project = sample_name,
                              min.cells = 3, min.features = 200)

  idx             <- match(sample_name, sample_map$gsm)
  obj$sample      <- sample_name
  obj$stage       <- sample_map$stage[idx]
  obj$patient     <- sample_map$patient[idx]
  obj$timepoint   <- sample_map$timepoint[idx]
  cells_raw       <- ncol(obj)

  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj$complexity      <- log10(obj$nFeature_RNA) /
                           log10(obj$nCount_RNA)

  min_genes   <- adaptive_lower(obj$nFeature_RNA, 200)
  max_genes   <- adaptive_upper(obj$nFeature_RNA)
  min_counts  <- adaptive_lower(obj$nCount_RNA,   500)
  max_counts  <- adaptive_upper(obj$nCount_RNA)

  obj <- subset(obj,
                subset = nFeature_RNA >= min_genes  &
                           nFeature_RNA <= max_genes  &
                           nCount_RNA   >= min_counts &
                           nCount_RNA   <= max_counts &
                           complexity   >= 0.60)

  cells_after_qc <- ncol(obj)

  sce <- as.SingleCellExperiment(obj)
  set.seed(42)
  sce        <- scDblFinder(sce)
  obj$dbl_class <- sce$scDblFinder.class
  doublets   <- sum(obj$dbl_class == "doublet")
  obj        <- subset(obj, dbl_class == "singlet")

  cat(sprintf("  Raw=%d QC=%d Final=%d Doublets=%d\n",
              cells_raw, cells_after_qc, ncol(obj), doublets))

  saveRDS(obj, clean_path)

  summary_df <- rbind(summary_df, data.frame(
    sample           = sample_name,
    stage            = sample_map$stage[idx],
    patient          = sample_map$patient[idx],
    timepoint        = sample_map$timepoint[idx],
    cells_raw        = cells_raw,
    cells_after_qc   = cells_after_qc,
    doublets_removed = doublets,
    cells_final      = ncol(obj)
  ))

  rm(obj, sce, mat); gc()
}

if (nrow(summary_df) > 0) {
  write.csv(summary_df,
            "../results/QC/tables/QC_summary_per_sample.csv",
            row.names = FALSE)
}

cat("Section A complete\n\n")

# ============================================================
# SECTION B: Per-patient clustering for plasma cell identification
# Uses Harmony to correct T1 vs T2 batch within each patient
# After running, inspect dotplots and update pc_clusters_final above
# ============================================================

cat("Section B: Per-patient clustering\n")

dir.create("../objects/patient_objects",    recursive = TRUE,
           showWarnings = FALSE)
dir.create("../results/patient_clustering", recursive = TRUE,
           showWarnings = FALSE)

patient_map <- list(
  P1 = list(T1 = "GSM8369868_MGUS_1", T2 = "GSM8369869_MGUS_2",
             stage_T1 = "MGUS", stage_T2 = "MGUS",
             progression = "Stable"),
  P2 = list(T1 = "GSM8369870_MGUS_3", T2 = "GSM8369871_MGUS_4",
             stage_T1 = "MGUS", stage_T2 = "MGUS",
             progression = "Stable"),
  P3 = list(T1 = "GSM8369872_MGUS_5", T2 = "GSM8369874_SMM_1",
             stage_T1 = "MGUS", stage_T2 = "SMM",
             progression = "Progressor"),
  P4 = list(T1 = "GSM8369873_MGUS_6", T2 = "GSM8369878_MM_1",
             stage_T1 = "MGUS", stage_T2 = "MM",
             progression = "Progressor"),
  P5 = list(T1 = "GSM8369875_SMM_2",  T2 = "GSM8369879_MM_2",
             stage_T1 = "SMM",  stage_T2 = "MM",
             progression = "Progressor"),
  P6 = list(T1 = "GSM8369876_SMM_3",  T2 = "GSM8369880_MM_3",
             stage_T1 = "SMM",  stage_T2 = "MM",
             progression = "Progressor"),
  P7 = list(T1 = "GSM8369877_SMM_4",  T2 = "GSM8369881_MM_4",
             stage_T1 = "SMM",  stage_T2 = "MM",
             progression = "Progressor")
)

pc_markers <- c(
  "MZB1","TNFRSF17","SDC1","JCHAIN","TXNDC5","DERL3","FKBP11",
  "IGHG1","IGHG3","IGHA1","IGHA2","IGKC","IGLC2","IRF4","CCND1",
  "MKI67","MS4A1","CD79A","PAX5","CD19","IGHD",
  "CD3D","CD3E","CD4","CD8A","NKG7","GNLY",
  "LYZ","CD14","S100A8","HBB","PPBP"
)

for (pid in names(patient_map)) {

  harmony_path <- paste0("../objects/patient_objects/",
                          pid, "_processed_harmony.rds")
  if (file.exists(harmony_path)) {
    cat("Already clustered:", pid, "\n")
    next
  }

  cat("Clustering:", pid, "\n")
  pm      <- patient_map[[pid]]
  out_dir <- paste0("../results/patient_clustering/", pid, "/")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  obj1 <- readRDS(paste0("../objects/clean_samples/",
                           pm$T1, "_clean.rds"))
  obj1$timepoint   <- "T1"
  obj1$stage       <- pm$stage_T1
  obj1$patient     <- pid
  obj1$progression <- pm$progression
  obj1$sample_id   <- pm$T1

  obj2 <- readRDS(paste0("../objects/clean_samples/",
                           pm$T2, "_clean.rds"))
  obj2$timepoint   <- "T2"
  obj2$stage       <- pm$stage_T2
  obj2$patient     <- pid
  obj2$progression <- pm$progression
  obj2$sample_id   <- pm$T2

  obj <- merge(obj1, obj2, project = pid)
  obj <- JoinLayers(obj)
  rm(obj1, obj2); gc()

  # SCTransform without regressing MT to preserve plasma cell signal
  obj <- SCTransform(obj, method = "glmGamPoi",
                      variable.features.n = 3000,
                      verbose = FALSE)

  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)

  # Harmony corrects T1 vs T2 batch within patient while
  # preserving biological cell type differences
  obj <- RunHarmony(obj, group.by.vars  = "sample_id",
                     reduction      = "pca",
                     reduction.save = "harmony",
                     verbose        = FALSE)

  obj <- RunUMAP(obj,       reduction = "harmony",
                  dims = 1:20, verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "harmony",
                        dims = 1:20, verbose = FALSE)
  obj <- FindClusters(obj,  resolution = 1.2, verbose = FALSE)

  n_clust <- length(unique(obj$seurat_clusters))
  cat("  Clusters:", n_clust, "\n")

  DefaultAssay(obj) <- "SCT"
  markers_present   <- pc_markers[pc_markers %in% rownames(obj)]

  p_dot <- DotPlot(obj, features = markers_present,
                    group.by = "seurat_clusters") +
    theme(axis.text.x = element_text(angle = 45,
                                       hjust = 1, size = 8)) +
    labs(title = paste(pid, "identify plasma cell clusters"))

  ggsave(paste0(out_dir, pid, "_dotplot.png"),
         p_dot,
         width  = 22,
         height = max(6, n_clust * 0.45 + 2),
         dpi    = 250)

  top20 <- FindAllMarkers(obj, only.pos = TRUE,
                            min.pct = 0.2,
                            logfc.threshold = 0.3,
                            verbose = FALSE) %>%
    group_by(cluster) %>%
    slice_max(avg_log2FC, n = 20)

  write.csv(top20,
            paste0(out_dir, pid, "_top20_markers.csv"),
            row.names = FALSE)

  saveRDS(obj, harmony_path)
  cat("  Saved:", harmony_path, "\n")

  rm(obj, top20); gc()
}

cat("Section B complete\n")
cat("Inspect dotplots in results/patient_clustering/ and update\n")
cat("pc_clusters_final at the top of this script before running\n")
cat("Sections C and D\n\n")

# ============================================================
# SECTION C: Annotate patients and HD reference, run inferCNV
# ============================================================

cat("Section C: Annotation and inferCNV\n")

dir.create("../objects/patient_objects", recursive = TRUE,
           showWarnings = FALSE)

# Annotate each patient object using manually identified
# plasma cell clusters from pc_clusters_final
ref_types <- c("CD4_T","CD8_T","NK","Monocyte")

for (pid in patients) {

  ann_path <- paste0("../objects/patient_objects/",
                      pid, "_annotated.rds")
  if (file.exists(ann_path)) {
    cat("Already annotated:", pid, "\n")
    next
  }

  cat("Annotating:", pid, "\n")

  obj <- readRDS(paste0("../objects/patient_objects/",
                          pid, "_processed_harmony.rds"))
  DefaultAssay(obj) <- "RNA"
  obj <- JoinLayers(obj)
  counts <- LayerData(obj, assay = "RNA", layer = "counts")

  all_cl <- levels(obj$seurat_clusters)

  pct_fn <- function(g, cells) {
    if (g %in% rownames(counts))
      round(mean(counts[g, cells] > 0) * 100, 1)
    else 0
  }

  results <- lapply(all_cl, function(cl) {
    cells <- colnames(obj)[obj$seurat_clusters == cl]
    data.frame(
      cluster      = cl,
      n_cells      = length(cells),
      T_CD4        = (pct_fn("CD3D",cells) + pct_fn("CD4",cells))   / 2,
      T_CD8        = (pct_fn("CD3D",cells) + pct_fn("CD8A",cells))  / 2,
      NK           = (pct_fn("NKG7",cells) + pct_fn("GNLY",cells))  / 2,
      Mono         = (pct_fn("LYZ",cells)  + pct_fn("CD14",cells))  / 2,
      B_cell       = (pct_fn("MS4A1",cells)+ pct_fn("CD79A",cells)) / 2,
      Plasma       = (pct_fn("MZB1",cells) + pct_fn("TNFRSF17",cells)) / 2,
      MZB1_pct     = pct_fn("MZB1",cells),
      TNFRSF17_pct = pct_fn("TNFRSF17",cells),
      Ery          = pct_fn("ALAS2",cells),
      Cycling      = pct_fn("MKI67",cells),
      Platelet     = pct_fn("PPBP",cells),
      pDC          = pct_fn("LILRA4",cells),
      HSC          = pct_fn("CD34",cells)
    )
  })

  df <- bind_rows(results) %>%
    mutate(cell_type = case_when(
      cluster %in% pc_clusters_final[[pid]] ~ "Plasma_cell",
      Platelet > 30                          ~ "Platelet",
      pDC > 30                               ~ "pDC",
      HSC > 25                               ~ "HSC_progenitor",
      Ery > 70                               ~ "Erythroid",
      Cycling > 25 & MZB1_pct < 20         ~ "Cycling",
      MZB1_pct > 50                          ~ "Plasma_cell",
      B_cell > 25 & MZB1_pct < 20          ~ "B_cell",
      NK > 25 & T_CD4 < 15 & T_CD8 < 15   ~ "NK",
      T_CD8 > 15                             ~ "CD8_T",
      T_CD4 > 15                             ~ "CD4_T",
      Mono > 25                              ~ "Monocyte",
      TRUE                                   ~ "Unknown"
    ))

  obj$cell_type <- df$cell_type[
    match(as.character(obj$seurat_clusters),
          as.character(df$cluster))
  ]

  cat(sprintf("  Plasma cells: %d | Ref cells: %d\n",
              sum(obj$cell_type == "Plasma_cell", na.rm = TRUE),
              sum(obj$cell_type %in% ref_types, na.rm = TRUE)))

  saveRDS(obj, ann_path)
  rm(obj, counts); gc()
}

# Annotate HD reference using the five healthy donor samples
hd_ann_path <- "../objects/HD_reference_annotated.rds"

if (!file.exists(hd_ann_path)) {

  cat("Building HD reference\n")

  hd_samples <- c("GSM8369863_HD1","GSM8369864_HD2",
                   "GSM8369865_HD3","GSM8369866_HD4",
                   "GSM8369867_HD5")

  hd_list <- lapply(hd_samples, function(s) {
    obj <- readRDS(paste0("../objects/clean_samples/",
                           s, "_clean.rds"))
    obj$sample_id <- s
    obj
  })

  hd_ref <- merge(hd_list[[1]], hd_list[-1],
                   project = "HD_reference")
  hd_ref <- JoinLayers(hd_ref)
  rm(hd_list); gc()

  hd_ref <- SCTransform(hd_ref, method = "glmGamPoi",
                          verbose = FALSE)
  hd_ref <- RunPCA(hd_ref, npcs = 30, verbose = FALSE)
  hd_ref <- RunHarmony(hd_ref, group.by.vars = "sample_id",
                        reduction = "pca",
                        reduction.save = "harmony",
                        verbose = FALSE)
  hd_ref <- RunUMAP(hd_ref, reduction = "harmony",
                     dims = 1:20, verbose = FALSE)
  hd_ref <- FindNeighbors(hd_ref, reduction = "harmony",
                            dims = 1:20, verbose = FALSE)
  hd_ref <- FindClusters(hd_ref, resolution = 0.8,
                           verbose = FALSE)

  saveRDS(hd_ref, "../objects/HD_reference.rds")

  # Annotate HD clusters using same marker-based approach
  DefaultAssay(hd_ref) <- "RNA"
  hd_ref <- JoinLayers(hd_ref)
  counts_hd <- LayerData(hd_ref, assay = "RNA", layer = "counts")

  pct_hd <- function(g, cells) {
    if (g %in% rownames(counts_hd))
      round(mean(counts_hd[g, cells] > 0) * 100, 1)
    else 0
  }

  results_hd <- lapply(levels(hd_ref$seurat_clusters), function(cl) {
    cells <- colnames(hd_ref)[hd_ref$seurat_clusters == cl]
    data.frame(
      cluster  = cl, n = length(cells),
      CD4_T    = (pct_hd("CD3D",cells)+pct_hd("CD4",cells))/2,
      CD8_T    = (pct_hd("CD3D",cells)+pct_hd("CD8A",cells))/2,
      NK       = (pct_hd("NKG7",cells)+pct_hd("GNLY",cells))/2,
      Mono     = (pct_hd("LYZ",cells) +pct_hd("CD14",cells))/2,
      B_cell   = (pct_hd("MS4A1",cells)+pct_hd("CD79A",cells))/2,
      MZB1     = pct_hd("MZB1",cells),
      Ery      = pct_hd("ALAS2",cells),
      Platelet = pct_hd("PPBP",cells),
      pDC      = pct_hd("LILRA4",cells),
      HSC      = pct_hd("CD34",cells),
      Cycling  = pct_hd("MKI67",cells)
    )
  })

  df_hd <- bind_rows(results_hd) %>%
    mutate(cell_type = case_when(
      Platelet > 30                           ~ "Platelet",
      pDC > 30                                ~ "pDC",
      HSC > 25                                ~ "HSC_progenitor",
      Ery > 70                                ~ "Erythroid",
      Cycling > 25 & MZB1 < 20              ~ "Cycling",
      MZB1 > 50                               ~ "Plasma_cell_HD",
      B_cell > 25 & MZB1 < 20               ~ "B_cell",
      NK > 25 & CD4_T < 15 & CD8_T < 15    ~ "NK",
      CD8_T > 15                              ~ "CD8_T",
      CD4_T > 15                              ~ "CD4_T",
      Mono > 25                               ~ "Monocyte",
      TRUE                                    ~ "Unknown"
    ))

  hd_ref$cell_type <- df_hd$cell_type[
    match(as.character(hd_ref$seurat_clusters),
          as.character(df_hd$cluster))
  ]

  cat("HD cell type distribution:\n")
  print(sort(table(hd_ref$cell_type), decreasing = TRUE))

  saveRDS(hd_ref, hd_ann_path)
  cat("HD reference saved\n")
  rm(hd_ref, counts_hd); gc()
}

# Run inferCNV per patient using patient T/NK/Mono cells plus
# pooled HD T/NK/Mono cells as diploid reference
cat("Running inferCNV\n")

hd_ref    <- readRDS(hd_ann_path)
DefaultAssay(hd_ref) <- "RNA"
hd_ref    <- JoinLayers(hd_ref)

hd_types  <- c("CD4_T","CD8_T","NK","Monocyte")
hd_cells  <- colnames(hd_ref)[
  !is.na(hd_ref$cell_type) & hd_ref$cell_type %in% hd_types
]
set.seed(42)
if (length(hd_cells) > 800) hd_cells <- sample(hd_cells, 800)
cat(sprintf("HD reference cells: %d\n", length(hd_cells)))

counts_hd_full <- GetAssayData(hd_ref[, hd_cells],
                                 assay = "RNA", layer = "counts")
rm(hd_ref); gc()

for (pid in patients) {

  out_dir    <- paste0("../results/inferCNV/", pid)
  final_file <- paste0(out_dir, "/infercnv.observations.txt")

  if (file.exists(final_file)) {
    cat("inferCNV already complete:", pid, "\n")
    next
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cat("Running inferCNV:", pid, "\n")

  obj <- readRDS(paste0("../objects/patient_objects/",
                          pid, "_annotated.rds"))
  DefaultAssay(obj) <- "RNA"
  obj <- JoinLayers(obj)

  pc_cells <- colnames(obj)[
    !is.na(obj$cell_type) & obj$cell_type == "Plasma_cell"
  ]
  pat_ref_cells <- colnames(obj)[
    !is.na(obj$cell_type) & obj$cell_type %in% hd_types
  ]
  set.seed(42)
  if (length(pat_ref_cells) > 400)
    pat_ref_cells <- sample(pat_ref_cells, 400)

  cat(sprintf("  Plasma=%d PatRef=%d HDRef=%d\n",
              length(pc_cells), length(pat_ref_cells),
              length(hd_cells)))

  if (length(pc_cells) < 30) {
    cat("  Too few plasma cells, skipping\n")
    next
  }

  counts_pc  <- GetAssayData(obj[, pc_cells],
                               assay = "RNA", layer = "counts")
  counts_ref <- GetAssayData(obj[, pat_ref_cells],
                               assay = "RNA", layer = "counts")
  counts_hd  <- counts_hd_full

  common_genes <- Reduce(intersect, list(
    rownames(counts_pc),
    rownames(counts_ref),
    rownames(counts_hd)
  ))

  counts_combined <- cbind(
    counts_pc[common_genes, ],
    counts_ref[common_genes, ],
    counts_hd[common_genes, ]
  )

  # Plasma cells labelled by timepoint so T1 and T2 appear
  # as separate observation groups in the inferCNV heatmap
  pc_labels      <- paste0("PC_", obj$timepoint[pc_cells])
  pat_ref_labels <- rep("Patient_ref", length(pat_ref_cells))
  hd_ref_labels  <- rep("HD_ref",      length(hd_cells))

  annotations_df <- data.frame(
    cell_type = c(pc_labels, pat_ref_labels, hd_ref_labels),
    row.names  = c(pc_cells, pat_ref_cells, hd_cells)
  )

  ann_file   <- paste0(out_dir, "/", pid, "_annotations.txt")
  count_file <- paste0(out_dir, "/", pid, "_counts.txt")

  write.table(annotations_df, ann_file,
              sep = "\t", quote = FALSE, col.names = FALSE)
  write.table(as.matrix(counts_combined), count_file,
              sep = "\t", quote = FALSE)

  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = count_file,
    annotations_file  = ann_file,
    delim             = "\t",
    gene_order_file   = gene_order_file,
    ref_group_names   = c("Patient_ref","HD_ref")
  )

  infercnv_obj <- infercnv::run(
    infercnv_obj,
    cutoff                            = 0.1,
    out_dir                           = out_dir,
    cluster_by_groups                 = TRUE,
    denoise                           = TRUE,
    HMM                               = TRUE,
    HMM_type                          = "i6",
    analysis_mode                     = "subclusters",
    tumor_subcluster_partition_method = "leiden",
    num_threads                       = 8,
    plot_steps                        = FALSE,
    no_prelim_plot                    = TRUE
  )

  saveRDS(infercnv_obj,
          paste0(out_dir, "/", pid, "_infercnv_obj.rds"))

  cat("  Done:", pid, "\n")
  rm(obj, infercnv_obj, counts_pc, counts_ref,
     counts_combined); gc()
}

cat("Section C complete\n\n")

# ============================================================
# SECTION D: Clonal evolution and Sankey diagrams
# ============================================================

cat("Section D: Clonal evolution\n")

dirs <- c(
  "../results/clonal_evolution",
  "../results/clonal_evolution/trees",
  "../results/clonal_evolution/CNV_fingerprints",
  "../results/clonal_evolution/DEG_clones",
  "../results/CNV_plots",
  "../results/clone_tracking"
)
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

gene_order <- read.table(
  gene_order_file,
  header           = FALSE,
  col.names        = c("gene","chr","start","end"),
  stringsAsFactors = FALSE
)

evolution_summary <- list()

for (pid in patients) {

  cat("\nClonal evolution:", pid, "\n")

  rds_file <- paste0("../results/inferCNV/", pid,
                      "/", pid, "_infercnv_obj.rds")
  if (!file.exists(rds_file)) {
    cat("  inferCNV object missing, skipping\n")
    next
  }

  infercnv_obj <- readRDS(rds_file)
  cnv_matrix   <- infercnv_obj@expr.data
  pc_cols      <- grep("^PC_", colnames(cnv_matrix), value = TRUE)

  if (length(pc_cols) == 0) {
    cat("  No PC_ columns, skipping\n")
    next
  }

  cnv_pc <- cnv_matrix[, pc_cols, drop = FALSE]

  groupings_file <- paste0(
    "../results/inferCNV/", pid,
    "/infercnv.19_HMM_pred.Bayes_Net.Pnorm_0.5",
    ".observation_groupings.txt"
  )
  if (!file.exists(groupings_file)) {
    cat("  Groupings file missing, skipping\n")
    next
  }

  groupings_raw <- read.table(groupings_file, header = TRUE,
                                sep = "", quote = "\"",
                                stringsAsFactors = FALSE)
  groupings <- data.frame(
    cell     = rownames(groupings_raw),
    subclone = groupings_raw$Dendrogram.Group,
    stringsAsFactors = FALSE
  )

  # Attach timepoint from annotated patient object
  obj     <- readRDS(paste0("../objects/patient_objects/",
                              pid, "_annotated.rds"))
  pc_meta <- obj@meta.data[
    !is.na(obj$cell_type) & obj$cell_type == "Plasma_cell", ]
  tp_map  <- setNames(pc_meta$timepoint,
                       paste0("PC_", rownames(pc_meta)))

  groupings$timepoint <- tp_map[groupings$cell]
  groupings           <- groupings[!is.na(groupings$timepoint), ]

  if (nrow(groupings) == 0) {
    cat("  No timepoint-matched cells, skipping\n")
    next
  }

  # Order genes by chromosome and position
  gene_meta <- gene_order[
    gene_order$gene %in% rownames(cnv_pc) &
      gene_order$chr %in% chr_order, , drop = FALSE
  ]
  gene_meta <- gene_meta[order(
    match(gene_meta$chr, chr_order), gene_meta$start
  ), ]
  cnv_pc    <- cnv_pc[gene_meta$gene, , drop = FALSE]

  # Compute mean CNV profile per subclone
  subclones    <- unique(groupings$subclone)
  profile_list <- list()

  for (sc in subclones) {
    cells_sc <- intersect(
      groupings$cell[groupings$subclone == sc],
      colnames(cnv_pc)
    )
    if (length(cells_sc) == 0) next
    profile_list[[sc]] <- if (length(cells_sc) == 1)
      cnv_pc[, cells_sc]
    else
      rowMeans(cnv_pc[, cells_sc, drop = FALSE])
  }

  if (length(profile_list) == 0) {
    cat("  No subclone profiles, skipping\n")
    next
  }

  profile_mat <- do.call(cbind, profile_list)

  sc_meta <- groupings %>%
    group_by(subclone, timepoint) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    group_by(timepoint) %>%
    mutate(pct = round(n_cells / sum(n_cells) * 100, 2)) %>%
    ungroup()

  cat("  Subclone composition:\n")
  print(as.data.frame(sc_meta))

  # CNV burden per subclone: mean squared deviation from diploid
  sc_burden <- sapply(colnames(profile_mat), function(sc) {
    mean((profile_mat[, sc] - 1)^2)
  })

  sc_tp_lookup <- sc_meta %>%
    group_by(subclone) %>%
    summarise(
      timepoints  = paste(sort(unique(timepoint)), collapse = "+"),
      total_cells = sum(n_cells),
      .groups     = "drop"
    )

  sc_burden_df <- data.frame(
    subclone   = names(sc_burden),
    cnv_burden = sc_burden,
    stringsAsFactors = FALSE
  ) %>% left_join(sc_tp_lookup, by = "subclone")

  # Pairwise CNV correlation between all subclones
  cor_all <- cor(profile_mat, method = "pearson")

  write.csv(as.data.frame(cor_all),
            paste0("../results/clonal_evolution/",
                    pid, "_all_subclone_correlation.csv"))

  # Also save as clone_tracking file for Sankey step
  write.csv(as.data.frame(cor_all),
            paste0("../results/clone_tracking/",
                    pid, "_clone_correlation.csv"))

  # Neighbour-joining phylogenetic tree rooted at lowest burden clone
  dist_mat       <- as.dist(1 - cor_all)
  nj_tree        <- nj(dist_mat)
  root_clone     <- names(which.min(sc_burden))

  if (root_clone %in% nj_tree$tip.label) {
    nj_tree_rooted <- root(nj_tree, outgroup = root_clone,
                            resolve.root = TRUE)
  } else {
    nj_tree_rooted <- nj_tree
  }

  tip_tp <- sapply(nj_tree_rooted$tip.label, function(sc) {
    paste(sort(unique(
      sc_meta$timepoint[sc_meta$subclone == sc]
    )), collapse = "+")
  })
  tip_colors <- ifelse(tip_tp == "T1",    "#4393C3",
                ifelse(tip_tp == "T2",    "#D6604D",
                ifelse(tip_tp == "T1+T2", "#7B3294", "grey50")))

  png(paste0("../results/clonal_evolution/trees/",
              pid, "_NJ_tree.png"),
      width = 900, height = 700, res = 150)
  plot(nj_tree_rooted, type = "phylogram",
       tip.color = tip_colors, cex = 0.9, edge.width = 2,
       main = paste(pid, "Subclone phylogeny"))
  legend("topright",
         legend = c("T1 only","T2 only","Shared T1+T2"),
         fill   = c("#4393C3","#D6604D","#7B3294"),
         bty    = "n", cex = 0.8)
  dev.off()

  # Classify T2 clones by evolutionary mode using CNV correlation
  # threshold >= 0.8 with single T1 match = linear descent
  # threshold >= 0.5 with multiple matches = branching
  # threshold < 0.4 = newly arising with no T1 precursor
  t1_clones <- intersect(
    unique(groupings$subclone[groupings$timepoint == "T1"]),
    colnames(profile_mat)
  )
  t2_clones <- intersect(
    unique(groupings$subclone[groupings$timepoint == "T2"]),
    colnames(profile_mat)
  )

  evolution_calls <- list()

  if (length(t1_clones) > 0 && length(t2_clones) > 0) {

    cor_t1t2 <- cor_all[t1_clones, t2_clones, drop = FALSE]

    for (t2c in t2_clones) {

      cor_to_t1  <- cor_t1t2[, t2c]
      max_cor    <- max(cor_to_t1)
      best_t1    <- names(which.max(cor_to_t1))
      n_high_cor <- sum(cor_to_t1 > 0.7)
      n_mod_cor  <- sum(cor_to_t1 > 0.4 & cor_to_t1 <= 0.7)

      evo_type <- if (max_cor >= 0.8 && n_high_cor == 1) {
        "Linear_descent"
      } else if (max_cor >= 0.5 && n_high_cor > 1) {
        "Branching"
      } else if (max_cor >= 0.4 && n_mod_cor >= 1) {
        "Branching"
      } else if (max_cor < 0.4) {
        "Newly_arising"
      } else {
        "Uncertain"
      }

      t1_pct <- sc_meta$pct[sc_meta$subclone == best_t1 &
                               sc_meta$timepoint == "T1"]
      t2_pct <- sc_meta$pct[sc_meta$subclone == t2c &
                               sc_meta$timepoint == "T2"]
      t1_pct <- if (length(t1_pct) == 0) 0 else t1_pct[1]
      t2_pct <- if (length(t2_pct) == 0) 0 else t2_pct[1]

      evolution_calls[[t2c]] <- data.frame(
        patient       = pid,
        t2_clone      = t2c,
        best_t1_match = best_t1,
        max_cor       = round(max_cor, 3),
        n_high_cor    = n_high_cor,
        evo_type      = evo_type,
        t1_pct        = t1_pct,
        t2_pct        = t2_pct,
        pct_change    = round(t2_pct - t1_pct, 2),
        cnv_burden    = round(sc_burden[t2c], 6),
        stringsAsFactors = FALSE
      )
    }
  }

  evo_df <- bind_rows(evolution_calls)

  if (nrow(evo_df) > 0) {
    cat("  Evolutionary mode per T2 clone:\n")
    print(as.data.frame(evo_df))
    write.csv(evo_df,
              paste0("../results/clonal_evolution/",
                      pid, "_evolution_calls.csv"),
              row.names = FALSE)
  }

  # Identify shared vs private chromosomal events
  # Threshold 0.01 accounts for inferCNV Bayesian denoising
  # which pulls values toward diploid baseline
  threshold <- 0.01

  chr_profiles <- lapply(chr_order, function(ch) {
    g <- intersect(gene_meta$gene[gene_meta$chr == ch],
                   rownames(profile_mat))
    if (length(g) < 5) return(NULL)
    colMeans(profile_mat[g, , drop = FALSE]) - 1
  })
  names(chr_profiles) <- chr_order
  chr_profiles        <- chr_profiles[!sapply(chr_profiles, is.null)]
  chr_mat             <- do.call(rbind, chr_profiles)

  shared_gains  <- rownames(chr_mat)[
    apply(chr_mat, 1, function(x) all(x >  threshold))
  ]
  shared_losses <- rownames(chr_mat)[
    apply(chr_mat, 1, function(x) all(x < -threshold))
  ]

  private_df <- bind_rows(lapply(colnames(chr_mat), function(sc) {
    others     <- setdiff(colnames(chr_mat), sc)
    if (length(others) == 0) return(NULL)
    sc_dev     <- chr_mat[, sc]
    other_mean <- if (length(others) == 1) chr_mat[, others] else
      rowMeans(chr_mat[, others, drop = FALSE])
    data.frame(
      subclone       = sc,
      private_gains  = paste(
        names(sc_dev)[sc_dev > threshold &
                        other_mean < threshold/2],
        collapse = ","
      ),
      private_losses = paste(
        names(sc_dev)[sc_dev < -threshold &
                        other_mean > -threshold/2],
        collapse = ","
      ),
      stringsAsFactors = FALSE
    )
  }))

  write.csv(private_df,
            paste0("../results/clonal_evolution/",
                    pid, "_private_CNV_events.csv"),
            row.names = FALSE)

  # Per-chromosome CNV fingerprint heatmap
  col_labels <- sapply(colnames(chr_mat), function(sc) {
    tp <- sc_tp_lookup$timepoints[sc_tp_lookup$subclone == sc]
    if (length(tp) == 0 || is.na(tp[1])) tp <- "?"
    paste0(sc, "\n(", tp[1], ")")
  })
  colnames(chr_mat) <- col_labels

  col_tp <- ifelse(grepl("T1\\+T2", colnames(chr_mat)), "Shared",
             ifelse(grepl("\\(T1\\)",   colnames(chr_mat)), "T1",
             ifelse(grepl("\\(T2\\)",   colnames(chr_mat)), "T2",
                    "Unknown")))

  ha_col <- HeatmapAnnotation(
    Timepoint = col_tp,
    col = list(Timepoint = c(
      "T1" = "#4393C3", "T2" = "#D6604D",
      "Shared" = "#7B3294", "Unknown" = "grey80"
    )),
    annotation_name_gp = gpar(fontsize = 8),
    simple_anno_size   = unit(3, "mm")
  )

  ht <- Heatmap(
    chr_mat,
    name             = "CNV\ndeviation",
    col              = colorRamp2(c(-0.15,0,0.15),
                                   c("#2166AC","white","#D6604D")),
    cluster_rows     = FALSE,
    cluster_columns  = FALSE,
    top_annotation   = ha_col,
    column_title     = paste(pid, "Subclone CNV fingerprints"),
    column_title_gp  = gpar(fontsize = 12, fontface = "bold"),
    row_names_gp     = gpar(fontsize = 9),
    column_names_gp  = gpar(fontsize = 8),
    column_names_rot = 45
  )

  png(paste0("../results/clonal_evolution/CNV_fingerprints/",
              pid, "_CNV_fingerprint_heatmap.png"),
      width = 1800, height = 1400, res = 180)
  draw(ht, padding = unit(c(5,25,5,5), "mm"))
  dev.off()
  cat("  CNV fingerprint heatmap saved\n")

  # DEG: dominant T2 clone vs other T2 plasma cells
  # Uses largest T2 clone by proportion as dominant
  if (nrow(evo_df) > 0) {

    dominant_clone <- evo_df$t2_clone[which.max(evo_df$t2_pct)]
    dom_cells      <- sub("^PC_", "",
                           groupings$cell[
                             groupings$subclone == dominant_clone &
                               groupings$timepoint == "T2"
                           ])
    other_cells    <- sub("^PC_", "",
                           groupings$cell[
                             groupings$subclone != dominant_clone &
                               groupings$timepoint == "T2"
                           ])

    if (length(dom_cells) >= 10 && length(other_cells) >= 10) {

      obj_pc <- subset(obj, cell_type == "Plasma_cell")
      DefaultAssay(obj_pc) <- "RNA"
      obj_pc <- JoinLayers(obj_pc)
      obj_pc <- NormalizeData(obj_pc, verbose = FALSE)
      obj_t2 <- subset(obj_pc, timepoint == "T2")

      obj_t2$clone_group <- ifelse(
        rownames(obj_t2@meta.data) %in% dom_cells,
        "Dominant", "Other_T2"
      )

      deg <- FindMarkers(obj_t2,
                          ident.1  = "Dominant",
                          ident.2  = "Other_T2",
                          group.by = "clone_group",
                          min.pct  = 0.1,
                          logfc.threshold = 0.2,
                          test.use = "wilcox",
                          verbose  = FALSE)

      deg$gene    <- rownames(deg)
      deg$patient <- pid

      write.csv(deg,
                paste0("../results/clonal_evolution/DEG_clones/",
                        pid, "_dominant_vs_other_T2.csv"),
                row.names = FALSE)

      cat(sprintf("  Dominant clone DEGs: %d up %d down\n",
                  sum(deg$p_val_adj < 0.05 & deg$avg_log2FC > 0,
                      na.rm = TRUE),
                  sum(deg$p_val_adj < 0.05 & deg$avg_log2FC < 0,
                      na.rm = TRUE)))

      rm(obj_pc, obj_t2); gc()
    }
  }

  evolution_summary[[pid]] <- list(
    evo_df        = evo_df,
    shared_gains  = shared_gains,
    shared_losses = shared_losses,
    private_df    = private_df,
    sc_burden_df  = sc_burden_df,
    sc_meta       = sc_meta
  )

  rm(infercnv_obj, cnv_matrix, cnv_pc, obj); gc()
}

# Cross-patient evolution summary
all_evo <- bind_rows(lapply(evolution_summary, function(x) x$evo_df))

if (nrow(all_evo) > 0) {

  write.csv(all_evo,
            "../results/clonal_evolution/all_patients_evolution_calls.csv",
            row.names = FALSE)

  evo_table <- all_evo %>%
    mutate(is_progressor = progressor[patient]) %>%
    group_by(is_progressor, evo_type) %>%
    summarise(n = n(), .groups = "drop")

  cat("\nEvolution type by progression status:\n")
  print(as.data.frame(evo_table))

  shared_summary <- bind_rows(lapply(names(evolution_summary), function(pid) {
    s <- evolution_summary[[pid]]
    data.frame(
      patient         = pid,
      shared_gains    = paste(s$shared_gains,  collapse = ","),
      shared_losses   = paste(s$shared_losses, collapse = ","),
      n_shared_gains  = length(s$shared_gains),
      n_shared_losses = length(s$shared_losses),
      stringsAsFactors = FALSE
    )
  }))

  write.csv(shared_summary,
            "../results/clonal_evolution/shared_CNV_events_summary.csv",
            row.names = FALSE)
}

# Sankey diagrams showing clonal flow from T1 to T2
cat("\nGenerating Sankey diagrams\n")

for (pid in patients) {

  groupings_file <- paste0(
    "../results/inferCNV/", pid,
    "/infercnv.19_HMM_pred.Bayes_Net.Pnorm_0.5",
    ".observation_groupings.txt"
  )
  cor_file <- paste0("../results/clone_tracking/",
                      pid, "_clone_correlation.csv")

  if (!file.exists(groupings_file)) {
    cat("Groupings missing for Sankey:", pid, "\n")
    next
  }

  groupings_raw <- read.table(groupings_file, header = TRUE,
                                sep = "", quote = "\"",
                                stringsAsFactors = FALSE)
  groupings <- data.frame(
    cell     = rownames(groupings_raw),
    subclone = groupings_raw$Dendrogram.Group,
    stringsAsFactors = FALSE
  )

  obj     <- readRDS(paste0("../objects/patient_objects/",
                              pid, "_annotated.rds"))
  pc_meta <- obj@meta.data[
    obj$cell_type == "Plasma_cell" & !is.na(obj$cell_type), ]
  tp_map  <- setNames(pc_meta$timepoint,
                       paste0("PC_", rownames(pc_meta)))
  rm(obj); gc()

  groupings$timepoint <- tp_map[groupings$cell]
  groupings           <- groupings[!is.na(groupings$timepoint), ]

  t2_sizes <- groupings %>%
    filter(timepoint == "T2") %>%
    group_by(subclone) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(desc(n))

  clone_label_map <- setNames(
    paste0("Clone_", seq_len(nrow(t2_sizes))),
    t2_sizes$subclone
  )

  if (file.exists(cor_file)) {
    cor_mat        <- as.matrix(read.csv(cor_file, row.names = 1))
    best_match_idx <- apply(cor_mat, 1, which.max)
    best_match_t2  <- colnames(cor_mat)[best_match_idx]
    best_r         <- apply(cor_mat, 1, max)
    t1_clone_names <- rownames(cor_mat)

    t1_label_map <- setNames(
      sapply(seq_along(t1_clone_names), function(i) {
        t2c <- best_match_t2[i]
        r   <- best_r[i]
        if (r >= 0.3 && t2c %in% names(clone_label_map))
          clone_label_map[t2c]
        else
          paste0("T1_only_", t1_clone_names[i])
      }),
      t1_clone_names
    )
  } else {
    t1_sizes <- groupings %>%
      filter(timepoint == "T1") %>%
      group_by(subclone) %>%
      summarise(n = n(), .groups = "drop") %>%
      arrange(desc(n))
    t1_label_map <- setNames(
      paste0("T1_Clone_", seq_len(nrow(t1_sizes))),
      t1_sizes$subclone
    )
  }

  groupings$clone_label <- ifelse(
    groupings$timepoint == "T2",
    clone_label_map[groupings$subclone],
    t1_label_map[groupings$subclone]
  )

  alluvial_df <- groupings %>%
    filter(!is.na(clone_label)) %>%
    group_by(timepoint, clone_label) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(timepoint) %>%
    mutate(pct = round(n / sum(n) * 100, 2)) %>%
    ungroup() %>%
    mutate(timepoint = factor(timepoint, levels = c("T1","T2")))

  write.csv(alluvial_df,
            paste0("../results/clone_tracking/",
                    pid, "_alluvial_final.csv"),
            row.names = FALSE)

  all_clones <- sort(unique(alluvial_df$clone_label))
  base_cols  <- c(
    "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
    "#A65628","#F781BF","#66C2A5","#FC8D62","#8DA0CB",
    "#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3",
    "#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E"
  )
  clone_cols <- setNames(
    colorRampPalette(base_cols)(length(all_clones)),
    all_clones
  )
  t1_only <- grep("^T1_only", all_clones, value = TRUE)
  clone_cols[t1_only] <- "#CCCCCC"

  p <- ggplot(alluvial_df,
              aes(x = timepoint, y = pct,
                  alluvium = clone_label,
                  stratum  = clone_label,
                  fill     = clone_label,
                  label    = clone_label)) +
    geom_flow(stat = "alluvium", lode.guidance = "frontback",
               alpha = 0.65, color = "white", linewidth = 0.2) +
    geom_stratum(alpha = 0.9, color = "white", linewidth = 0.3) +
    geom_text(stat = "stratum", size = 2.8,
               fontface = "bold", min.y = 4) +
    scale_fill_manual(values = clone_cols) +
    scale_y_continuous(labels = function(x) paste0(x,"%"),
                        expand = c(0,0)) +
    labs(
      title    = paste(pid, "Clonal evolution T1 to T2"),
      subtitle = stage_map[pid],
      x = "Timepoint", y = "Percent of plasma cells",
      fill = "Clone"
    ) +
    theme_classic(base_size = 13) +
    theme(legend.position = "right",
          legend.key.size = unit(0.45, "cm"),
          legend.text     = element_text(size = 7))

  ggsave(paste0("../results/CNV_plots/", pid,
                 "_alluvial_final.png"),
         p, width = 9, height = 6, dpi = 300)
  cat("  Sankey saved:", pid, "\n")
}

cat("\nSection D complete\n")
cat("All outputs in CNV_analysis/results/\n")