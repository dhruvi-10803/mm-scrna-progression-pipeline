# 01_QC_normalization_annotation.R
# Pipeline step 1: Load raw H5 files, perform QC, normalize with SCTransform,
# cluster, and annotate using Azimuth bone marrow reference.
# Input:  Raw H5 files per patient per timepoint in data/
# Output: Per-patient annotated Seurat objects in objects/Px/

library(Seurat)
library(ggplot2)
library(dplyr)
library(scDblFinder)
library(SingleCellExperiment)
library(Azimuth)
library(glmGamPoi)

setwd("~/Dhruvi/TME/scripts")

# Patient manifest mapping each patient to their raw H5 files and metadata
patient_manifest <- list(

  P1 = list(
    T1       = "../data/GSM8369868_MGUS_1.h5",
    T2       = "../data/GSM8369869_MGUS_2.h5",
    T1_name  = "MGUS_1",   T2_name  = "MGUS_2",
    patient  = "P1",       gender   = "Male",
    race     = "White/Caucasian",
    ig_type  = "IgG",      lc       = "Kappa",
    cyto     = "t(14q32)", mm_sub   = "MF",
    T1_stage = "MGUS",     T2_stage = "MGUS",
    group    = "Non_progressor"
  ),

  P2 = list(
    T1       = "../data/GSM8369870_MGUS_3.h5",
    T2       = "../data/GSM8369871_MGUS_4.h5",
    T1_name  = "MGUS_3",   T2_name  = "MGUS_4",
    patient  = "P2",       gender   = "Female",
    race     = "White/Caucasian",
    ig_type  = "IgG",      lc       = "Kappa",
    cyto     = "del(13q)", mm_sub   = "LB",
    T1_stage = "MGUS",     T2_stage = "MGUS",
    group    = "Non_progressor"
  ),

  P3 = list(
    T1       = "../data/GSM8369872_MGUS_5.h5",
    T2       = "../data/GSM8369874_SMM_1.h5",
    T1_name  = "MGUS_5",   T2_name  = "SMM_1",
    patient  = "P3",       gender   = "Male",
    race     = "White/Caucasian",
    ig_type  = "IgG",      lc       = "Lambda",
    cyto     = "NA",       mm_sub   = "CD2",
    T1_stage = "MGUS",     T2_stage = "SMM",
    group    = "Progressor"
  ),

  P4 = list(
    T1       = "../data/GSM8369873_MGUS_6.h5",
    T2       = "../data/GSM8369878_MM_1.h5",
    T1_name  = "MGUS_6",   T2_name  = "MM_1",
    patient  = "P4",       gender   = "Female",
    race     = "White/Caucasian",
    ig_type  = "IgA",      lc       = "Kappa",
    cyto     = "NA",       mm_sub   = "HY",
    T1_stage = "MGUS",     T2_stage = "MM",
    group    = "Progressor"
  ),

  P5 = list(
    T1       = "../data/GSM8369875_SMM_2.h5",
    T2       = "../data/GSM8369879_MM_2.h5",
    T1_name  = "SMM_2",        T2_name  = "MM_2",
    patient  = "P5",           gender   = "Female",
    race     = "White/Caucasian",
    ig_type  = "IgD",          lc       = "Kappa",
    cyto     = "del(17p13.1)", mm_sub   = "HY",
    T1_stage = "SMM",          T2_stage = "MM",
    group    = "Progressor"
  ),

  P6 = list(
    T1       = "../data/GSM8369876_SMM_3.h5",
    T2       = "../data/GSM8369880_MM_3.h5",
    T1_name  = "SMM_3",        T2_name  = "MM_3",
    patient  = "P6",           gender   = "Male",
    race     = "White/Caucasian",
    ig_type  = "IgG",          lc       = "Kappa",
    cyto     = "del(17p13.1)", mm_sub   = "HY",
    T1_stage = "SMM",          T2_stage = "MM",
    group    = "Progressor"
  ),

  P7 = list(
    T1       = "../data/GSM8369877_SMM_4.h5",
    T2       = "../data/GSM8369881_MM_4.h5",
    T1_name  = "SMM_4",   T2_name  = "MM_4",
    patient  = "P7",      gender   = "Female",
    race     = "Black",
    ig_type  = "IgG",     lc       = "Kappa",
    cyto     = "t(11;14)",mm_sub   = "CD2",
    T1_stage = "SMM",     T2_stage = "MM",
    group    = "Progressor"
  )
)

# Marker panel used for cluster inspection dotplots
marker_list <- list(
  T_cell         = c("CD3D","CD3E","CD3G","TRAC","CD2","CD7",
                     "IL7R","CCR7","TCF7","SELL"),
  CD8_T_cell     = c("CD8A","CD8B","GZMB","PRF1","GNLY",
                     "NKG7","CTSW"),
  Treg           = c("FOXP3","IL2RA","CTLA4","TIGIT","IKZF2"),
  NK_cell        = c("NKG7","GNLY","KLRD1","KLRF1","KLRB1",
                     "FCGR3A","GZMB","PRF1"),
  B_cell         = c("MS4A1","CD79A","CD79B","CD19","CD74",
                     "CD37","BANK1","BLK"),
  Plasma_cell    = c("MZB1","SDC1","XBP1","JCHAIN","DERL3",
                     "IGKC","IGLC1","IGHA1"),
  Classical_Mono = c("CD14","LYZ","S100A8","S100A9",
                     "VCAN","LGALS3"),
  NonClass_Mono  = c("FCGR3A","CX3CR1","MS4A7","LILRB1"),
  cDC            = c("CLEC10A","CD1C","FCER1A",
                     "HLA-DRA","HLA-DPB1"),
  pDC            = c("LILRA4","GZMB","IRF7","TCF4"),
  Neutrophil     = c("MPO","ELANE","AZU1","FCGR3B"),
  HSPC           = c("CD34","KIT","PROM1","GATA2","MEIS1"),
  Erythroid      = c("HBB","HBA1","HBA2","ALAS2",
                     "AHSP","GYPB"),
  Megakaryocyte  = c("PF4","PPBP","GP9","ITGA2B","TUBB1"),
  IFN_activated  = c("ISG15","IFI6","IFIT1","IFIT3",
                     "MX1","IRF7")
)

# Log to track runtime per patient
pipeline_log <- data.frame()

for (pid in names(patient_manifest)) {

  meta      <- patient_manifest[[pid]]
  start_time <- Sys.time()

  cat("Processing patient:", pid, "\n")

  obj_dir <- paste0("../objects/", pid)
  res_dir <- paste0("../results/", pid)

  dir.create(obj_dir,                     recursive = TRUE,
             showWarnings = FALSE)
  dir.create(paste0(res_dir, "/QC"),      recursive = TRUE,
             showWarnings = FALSE)
  dir.create(paste0(res_dir, "/figures"), recursive = TRUE,
             showWarnings = FALSE)
  dir.create(paste0(res_dir, "/tables"),  recursive = TRUE,
             showWarnings = FALSE)

  samples <- list(
    T1 = list(path      = meta$T1,
               name      = meta$T1_name,
               timepoint = "TimePoint_1",
               stage     = meta$T1_stage),
    T2 = list(path      = meta$T2,
               name      = meta$T2_name,
               timepoint = "TimePoint_2",
               stage     = meta$T2_stage)
  )

  seurat_list <- list()
  qc_df       <- data.frame()

  for (tp in names(samples)) {

    s <- samples[[tp]]
    cat("Loading:", s$name, "\n")

    mat <- Read10X_h5(s$path)

    obj <- CreateSeuratObject(
      counts       = mat,
      project      = s$name,
      min.cells    = 3,
      min.features = 200
    )

    # Add patient and sample metadata to every cell
    obj$sample    <- s$name
    obj$timepoint <- s$timepoint
    obj$stage     <- s$stage
    obj$patient   <- meta$patient
    obj$group     <- meta$group
    obj$gender    <- meta$gender
    obj$ig_type   <- meta$ig_type
    obj$lc        <- meta$lc
    obj$cyto      <- meta$cyto
    obj$mm_sub    <- meta$mm_sub

    obj[["percent.mt"]] <- PercentageFeatureSet(obj,
                                                 pattern = "^MT-")

    qc_df <- rbind(qc_df, data.frame(
      patient    = pid,
      sample     = s$name,
      timepoint  = s$timepoint,
      cells_raw  = ncol(obj),
      genes_med  = median(obj$nFeature_RNA),
      counts_med = median(obj$nCount_RNA),
      mito_med   = median(obj$percent.mt)
    ))

    # QC filtering: minimum 200 genes, maximum 30% mitochondrial reads
    # No upper gene/UMI threshold to preserve active immune populations
    cells_before <- ncol(obj)
    obj <- subset(obj,
                  subset = nFeature_RNA > 200 & percent.mt < 30)
    cat("After QC filter:", ncol(obj),
        "(removed", cells_before - ncol(obj), "cells)\n")

    # Doublet detection using scDblFinder, only singlets retained
    sce <- as.SingleCellExperiment(obj)
    sce <- scDblFinder(sce)
    obj$scDblFinder.class <- sce$scDblFinder.class
    before_dbl <- ncol(obj)
    obj <- subset(obj,
                  subset = scDblFinder.class == "singlet")
    cat("After doublet removal:", ncol(obj),
        "(removed", before_dbl - ncol(obj), "doublets)\n")

    seurat_list[[tp]] <- obj
    rm(mat, sce)
    gc()
  }

  write.csv(qc_df,
            paste0(res_dir, "/tables/", pid, "_QC_metrics.csv"),
            row.names = FALSE)

  # Merge T1 and T2 without batch correction to preserve
  # longitudinal biological variation between timepoints
  px <- merge(seurat_list[["T1"]], seurat_list[["T2"]],
              project = paste0(pid, "_merged"))

  cat("Merged object:", ncol(px), "total cells\n")
  print(table(px$timepoint))
  rm(seurat_list)
  gc()

  # Normalize using SCTransform with mitochondrial percent regressed out
  # glmGamPoi used for efficiency with single-cell count distributions
  cat("Running SCTransform\n")
  px <- SCTransform(px,
                    vars.to.regress = "percent.mt",
                    method          = "glmGamPoi",
                    verbose         = FALSE)

  saveRDS(px, paste0(obj_dir, "/", pid, "_after_sct.rds"))

  # PCA using top 30 PCs, UMAP and clustering on first 15 PCs
  # Dims 1:15 validated on elbow plot across patients
  cat("Running PCA, UMAP, clustering\n")
  px <- RunPCA(px, npcs = 30, verbose = FALSE)

  ggsave(paste0(res_dir, "/figures/ElbowPlot.png"),
         ElbowPlot(px, ndims = 30),
         width = 6, height = 4, dpi = 300)

  DIMS <- 1:15
  px   <- RunUMAP(px,       dims = DIMS, verbose = FALSE)
  px   <- FindNeighbors(px, dims = DIMS, verbose = FALSE)
  px   <- FindClusters(px,  resolution = 0.5, verbose = FALSE)

  n_clusters <- nlevels(px$seurat_clusters)
  cat("Clusters found:", n_clusters, "\n")

  saveRDS(px, paste0(obj_dir, "/", pid, "_after_clustering.rds"))

  # Annotate using Azimuth human bone marrow reference
  # az_l1 = broad cell type, az_l2 = fine-grained subtype
  # Cells with mapping score below 0.5 excluded downstream
  cat("Running Azimuth annotation\n")
  px <- RunAzimuth(px, reference = "bonemarrowref")

  if ("predicted.celltype.l1" %in% colnames(px@meta.data)) {
    px$az_l1       <- px$predicted.celltype.l1
    px$az_l2       <- px$predicted.celltype.l2
    px$az_l1_score <- px$predicted.celltype.l1.score
    px$az_l2_score <- px$predicted.celltype.l2.score
  } else {
    px$az_l1       <- px$celltype.l1
    px$az_l2       <- px$celltype.l2
    px$az_l1_score <- px$celltype.l1.score
    px$az_l2_score <- px$celltype.l2.score
  }

  if (!"mapping.score" %in% colnames(px@meta.data))
    px$mapping.score <- rep(1, ncol(px))

  cat("Azimuth L1 distribution:\n")
  print(sort(table(px$az_l1), decreasing = TRUE))

  saveRDS(px, paste0(obj_dir, "/", pid, "_after_azimuth.rds"))

  # Save UMAP panels for visual inspection of clusters and annotation
  DefaultAssay(px) <- "SCT"

  pdf(paste0(res_dir, "/figures/", pid, "_UMAP_panels.pdf"),
      width = 14, height = 10)

  print(DimPlot(px, label = TRUE, raster = FALSE) +
          ggtitle(paste0(pid, " clusters")))

  print(DimPlot(px, group.by = "timepoint", raster = FALSE) +
          ggtitle(paste0(pid, " T1 vs T2")))

  print(DimPlot(px, group.by = "az_l1",
                label = TRUE, repel = TRUE, raster = FALSE) +
          ggtitle(paste0(pid, " Azimuth L1")))

  print(DimPlot(px, group.by = "az_l2", raster = FALSE) +
          ggtitle(paste0(pid, " Azimuth L2")) +
          theme(legend.text = element_text(size = 6)))

  print(DimPlot(px, split.by = "timepoint",
                label = TRUE, raster = FALSE) +
          ggtitle(paste0(pid, " clusters by timepoint")))

  print(DimPlot(px, group.by = "az_l1",
                split.by = "timepoint", raster = FALSE) +
          ggtitle(paste0(pid, " Azimuth L1 by timepoint")) +
          theme(legend.text = element_text(size = 7)))

  dev.off()

  # DotPlot of full marker panel per cluster for manual curation
  full_markers <- unique(unlist(
    lapply(marker_list, function(x) x[x %in% rownames(px)])
  ))

  pdf(paste0(res_dir, "/", pid,
             "_DotPlot_clusters_full_markers.pdf"),
      width = 38,
      height = max(8, n_clusters * 0.5 + 4))

  print(DotPlot(px,
                features = full_markers,
                group.by = "seurat_clusters") +
          RotatedAxis() +
          ggtitle(paste0(pid, " clusters vs marker panel")) +
          theme(axis.text.x = element_text(size = 6),
                axis.text.y = element_text(size = 9)))

  dev.off()

  # Save average expression per cluster for marker inspection
  avg <- AverageExpression(
    px,
    features = full_markers,
    assays   = "SCT",
    group.by = "seurat_clusters"
  )$SCT

  write.csv(round(avg, 3),
            paste0(res_dir, "/tables/", pid,
                   "_avg_expression_clusters.csv"))

  # Save Azimuth L1 and L2 per cluster summaries
  az_summary <- px@meta.data %>%
    group_by(seurat_clusters, az_l1) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(seurat_clusters) %>%
    mutate(total = sum(n),
           pct   = round(n / total * 100, 1)) %>%
    ungroup() %>%
    arrange(seurat_clusters, desc(n))

  write.csv(az_summary,
            paste0(res_dir, "/tables/", pid,
                   "_azimuth_l1_cluster_summary.csv"),
            row.names = FALSE)

  az_l2_summary <- px@meta.data %>%
    group_by(seurat_clusters, az_l2) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(seurat_clusters, desc(n))

  write.csv(az_l2_summary,
            paste0(res_dir, "/tables/", pid,
                   "_azimuth_l2_cluster_summary.csv"),
            row.names = FALSE)

  # Final annotated object saved for downstream analysis
  saveRDS(px, paste0(obj_dir, "/", pid, "_annotated.rds"))

  end_time <- Sys.time()
  elapsed  <- round(difftime(end_time, start_time,
                              units = "mins"), 1)

  pipeline_log <- rbind(pipeline_log, data.frame(
    patient     = pid,
    group       = meta$group,
    T1_stage    = meta$T1_stage,
    T2_stage    = meta$T2_stage,
    n_clusters  = n_clusters,
    total_cells = ncol(px),
    time_mins   = as.numeric(elapsed)
  ))

  cat(pid, "complete in", elapsed, "mins\n\n")
  rm(px)
  gc()
}

write.csv(pipeline_log,
          "../results/pipeline_run_log.csv",
          row.names = FALSE)

cat("All patients complete\n")
print(pipeline_log)