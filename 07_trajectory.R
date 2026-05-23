# 07_trajectory.R (optional)
# Pipeline step 7: Pseudotime trajectory analysis using Monocle3.
# Runs per cell type across all patients using az_l2 annotations.
# Cells are ordered from T1 as root along a learned principal graph.
# Pseudotime-correlated genes are identified per patient and
# summarised across patients to find consistent progression signals.
# Input:  Per-patient scored objects from consolidation/module_scoring/Px_scored.rds
# Output: Trajectory figures in consolidation/figures/trajectory/
#         Pseudotime gene tables in consolidation/trajectory_output/

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

obj_dir <- "~/Dhruvi/consolidation/module_scoring/"
fig_dir <- "~/Dhruvi/consolidation/figures/trajectory/"
out_dir <- "~/Dhruvi/consolidation/trajectory_output/"

dir.create(path.expand(fig_dir), recursive = TRUE,
           showWarnings = FALSE)
dir.create(path.expand(out_dir), recursive = TRUE,
           showWarnings = FALSE)

patients <- paste0("P", 1:7)

traj_group <- c(
  P1 = "Stable_MGUS",  P2 = "Stable_MGUS",
  P3 = "MGUS_to_SMM",  P4 = "MGUS_to_MM",
  P5 = "SMM_to_MM",    P6 = "SMM_to_MM",
  P7 = "SMM_to_MM"
)

stable_pts    <- c("P1","P2")
progressor_pts <- c("P3","P4","P5","P6","P7")

MIN_CELLS <- 50

# Cell type groupings based on az_l2 labels
celltype_groups <- list(
  Monocyte  = c("CD14 Mono","CD16 Mono","Macrophage"),
  CD8_T     = c("CD8 Memory","CD8 Effector_1",
                 "CD8 Effector_2","CD8 Effector_3","CD8 Naive"),
  CD4_T     = c("CD4 Memory","CD4 Naive",
                 "CD4 Effector","MAIT"),
  NK        = c("NK","NK CD56+","NK Proliferating"),
  B_cell    = c("Memory B","Naive B","transitional B",
                 "pre B","pro B"),
  Plasma    = c("Plasma"),
  Erythroid = c("Late Eryth","Early Eryth","EMP"),
  HSPC      = c("GMP","HSC","LMPP","CLP",
                 "Prog Mk","BaEoMa")
)

# Module score columns to overlay on trajectory plots
module_overlays <- list(
  Monocyte  = c("Mono_Homeostatic","Mono_Classical_Inflammatory",
                 "Mono_Immunosuppressive","Mono_IFN_Activated",
                 "Mono_AntigenPresenting"),
  CD8_T     = c("CD8_NaiveStemlike","CD8_CytotoxicMemory",
                 "CD8_EffectorCytotoxic","CD8_Exhausted",
                 "CD8_IFN_Activated"),
  CD4_T     = c("CD4_Helper_Memory","CD4_Regulatory",
                 "CD4_Exhausted","CD4_IFN_Activated"),
  NK        = c("NK_Cytotoxic_Mature","NK_Immature_Tissue",
                 "NK_Dysfunctional","NK_IFN_Activated"),
  B_cell    = c("CD4_Helper_Memory"),
  Plasma    = c("Plasma_Normal","Plasma_Myeloma_Proliferating",
                 "Plasma_Myeloma_Survival","Plasma_IFN_Stress"),
  Erythroid = c("Erythroid_Progenitor","Erythroid_Mature",
                 "Erythroid_Stress_IFN"),
  HSPC      = c("HSPC_StemQuiescent","HSPC_GMP_Myeloid",
                 "HSPC_Lymphoid_CLP","HSPC_Megakaryocyte")
)

# Biologically relevant genes to test for pseudotime correlation
genes_of_interest <- list(
  Monocyte  = c("SOCS3","PER1","LUCAT1","BCL6","BTG2","NR4A1",
                 "KLF10","ISG15","CXCL8","JUN","BST2","MX1",
                 "HLA-DRA","HLA-DRB1","CX3CR1","IL1B"),
  CD8_T     = c("GZMK","GZMB","GZMH","PRF1","NKG7",
                 "PDCD1","LAG3","TIGIT","HAVCR2","TOX",
                 "TCF7","IL7R","CCL3","CCL4","CCL5",
                 "HIST1H1D","ISG15","MX1","IFITM1"),
  CD4_T     = c("FOXP3","IL2RA","CTLA4","PDCD1","LAG3",
                 "CD40LG","ICOS","IL7R","TCF7",
                 "HIST1H1D","ISG15","GZMK"),
  NK        = c("NKG7","GZMB","PRF1","GZMK",
                 "CCL3","CCL4","CCL5","CXCR4","XCL1",
                 "HIST1H1D","LAG3","TIGIT","ISG15",
                 "KLRD1","FCGR3A"),
  B_cell    = c("MS4A1","CD19","CD38","SDC1",
                 "MZB1","IGHG1","IGHM","AICDA",
                 "BCL6","CD27","CD24"),
  Plasma    = c("MZB1","XBP1","PRDM1","IRF4",
                 "MKI67","TOP2A","MCL1","BCL2",
                 "DKK1","CXCR4","IFI27","ISG15",
                 "BTG1","BTG2","TXNIP"),
  Erythroid = c("GATA1","KLF1","HBB","HBA1","GYPA",
                 "CD34","TFRC","IFI27","PARP9","ISG15",
                 "SAMD9L","BST2","ENG","ORC1"),
  HSPC      = c("CD34","MPO","ELANE","IL7R","DNTT",
                 "FLT3","CEBPA","SPI1","GATA1","GATA2",
                 "GP1BB","PF4","HLF","MECOM","AVP")
)

# Convert Seurat object to Monocle3 CellDataSet
seurat_to_cds <- function(seu) {
  DefaultAssay(seu) <- "RNA"
  seu      <- JoinLayers(seu)
  expr_mat <- GetAssayData(seu, assay = "RNA", layer = "counts")
  cell_meta <- seu@meta.data
  gene_meta <- data.frame(
    gene_short_name = rownames(expr_mat),
    row.names       = rownames(expr_mat)
  )
  new_cell_data_set(
    expression_data = expr_mat,
    cell_metadata   = cell_meta,
    gene_metadata   = gene_meta
  )
}

# Run full Monocle3 preprocessing and graph learning pipeline
run_monocle <- function(cds, reduction_dims = 30) {
  cds <- preprocess_cds(cds, num_dim = reduction_dims,
                         verbose = FALSE)
  cds <- reduce_dimension(cds, reduction_method = "UMAP",
                           preprocess_method = "PCA",
                           verbose = FALSE)
  cds <- cluster_cells(cds, verbose = FALSE)
  cds <- learn_graph(cds, verbose = FALSE,
                      use_partition = FALSE)
  cds
}

# Order cells in pseudotime using T1 cells as root nodes
order_from_t1 <- function(cds,
                            timepoint_col = "timepoint",
                            t1_label = "TimePoint_1") {
  t1_cells <- rownames(colData(cds))[
    colData(cds)[[timepoint_col]] == t1_label
  ]
  if (length(t1_cells) == 0) {
    cat("  No T1 cells found for rooting, using default\n")
    return(order_cells(cds, verbose = FALSE))
  }
  order_cells(cds, root_cells = t1_cells, verbose = FALSE)
}

# Compute Spearman correlation between pseudotime and gene expression
get_pseudotime_genes <- function(cds, genes_test) {
  pt   <- pseudotime(cds)
  expr <- as.matrix(
    GetAssayData(
      as.Seurat(cds, assay = "RNA"),
      assay = "RNA", layer = "data"
    )
  )
  genes_present <- intersect(genes_test, rownames(expr))
  if (length(genes_present) == 0) return(NULL)

  results <- lapply(genes_present, function(g) {
    gene_expr <- expr[g, names(pt)]
    gene_expr <- gene_expr[!is.na(pt)]
    pt_use    <- pt[!is.na(pt)]
    if (length(unique(gene_expr)) < 3) return(NULL)
    cor_res <- tryCatch(
      cor.test(pt_use, gene_expr, method = "spearman"),
      error = function(e) NULL
    )
    if (is.null(cor_res)) return(NULL)
    data.frame(
      gene      = g,
      rho       = round(cor_res$estimate, 3),
      p_val     = round(cor_res$p.value, 5),
      direction = ifelse(cor_res$estimate > 0,
                          "Increases along pseudotime",
                          "Decreases along pseudotime"),
      stringsAsFactors = FALSE
    )
  })

  bind_rows(Filter(Negate(is.null), results)) %>%
    arrange(p_val)
}

all_pt_genes <- list()
all_cds_meta <- list()

cat("Trajectory analysis started\n")
cat(format(Sys.time()), "\n\n")

for (ct_name in names(celltype_groups)) {

  cat("Cell type:", ct_name, "\n")

  az_labels  <- celltype_groups[[ct_name]]
  overlays   <- module_overlays[[ct_name]]
  genes_test <- genes_of_interest[[ct_name]]

  pt_plots_by_patient <- list()
  ct_pt_genes         <- list()

  for (pt in patients) {

    cat(sprintf("  Patient: %s\n", pt))

    obj_path <- file.path(path.expand(obj_dir),
                           paste0(pt, "_scored.rds"))

    if (!file.exists(obj_path)) {
      cat("  Object not found, skipping\n")
      next
    }

    obj <- tryCatch(readRDS(obj_path),
                    error = function(e) {
                      cat("  Load error:", e$message, "\n")
                      NULL
                    })
    if (is.null(obj)) next

    cells_keep <- WhichCells(obj,
                              expression = az_l2 %in% az_labels)

    if (length(cells_keep) < MIN_CELLS) {
      cat(sprintf("  Only %d cells, skipping\n",
                  length(cells_keep)))
      rm(obj); gc()
      next
    }

    sub_obj <- subset(obj, cells = cells_keep)
    rm(obj); gc()

    cds <- tryCatch(seurat_to_cds(sub_obj),
                    error = function(e) {
                      cat("  CDS error:", e$message, "\n")
                      NULL
                    })
    if (is.null(cds)) { rm(sub_obj); gc(); next }

    cds <- tryCatch(run_monocle(cds),
                    error = function(e) {
                      cat("  Monocle error:", e$message, "\n")
                      NULL
                    })
    if (is.null(cds)) { rm(sub_obj, cds); gc(); next }

    cds <- tryCatch(order_from_t1(cds),
                    error = function(e) {
                      cat("  Ordering error:", e$message, "\n")
                      NULL
                    })
    if (is.null(cds)) { rm(sub_obj, cds); gc(); next }

    # Compute pseudotime gene correlations
    pt_genes <- tryCatch(
      get_pseudotime_genes(cds, genes_test),
      error = function(e) NULL
    )

    if (!is.null(pt_genes) && nrow(pt_genes) > 0) {
      pt_genes$patient    <- pt
      pt_genes$traj_group <- traj_group[pt]
      pt_genes$cell_type  <- ct_name
      ct_pt_genes[[pt]]   <- pt_genes

      cat("  Top pseudotime genes:\n")
      top5 <- head(pt_genes[order(pt_genes$p_val), ], 5)
      for (i in seq_len(nrow(top5))) {
        cat(sprintf("    %-15s rho=%+.2f  %s\n",
                    top5$gene[i], top5$rho[i],
                    top5$direction[i]))
      }
    }

    # Save pseudotime metadata for downstream use
    cds_meta             <- as.data.frame(colData(cds))
    cds_meta$pseudotime  <- pseudotime(cds)
    cds_meta$patient     <- pt
    cds_meta$cell_type   <- ct_name
    all_cds_meta[[paste(ct_name, pt, sep = "_")]] <- cds_meta

    # Timepoint coloured trajectory plot
    p_tp <- plot_cells(
      cds,
      color_cells_by           = "timepoint",
      show_trajectory_graph    = TRUE,
      label_cell_groups        = FALSE,
      label_leaves             = FALSE,
      label_branch_points      = FALSE,
      cell_size                = 0.6,
      trajectory_graph_color   = "grey30",
      trajectory_graph_segment_size = 0.6
    ) +
      scale_colour_manual(
        values = c(TimePoint_1 = "#2196F3",
                   TimePoint_2 = "#F44336"),
        labels = c(TimePoint_1 = "T1 baseline",
                   TimePoint_2 = "T2 follow-up")
      ) +
      ggtitle(sprintf("%s %s", pt, traj_group[pt])) +
      theme(plot.title       = element_text(size = 9,
                                             face = "bold"),
            legend.position  = "bottom",
            legend.text      = element_text(size = 7))

    # Pseudotime gradient plot
    p_pt <- plot_cells(
      cds,
      color_cells_by           = "pseudotime",
      show_trajectory_graph    = TRUE,
      label_cell_groups        = FALSE,
      label_leaves             = FALSE,
      label_branch_points      = FALSE,
      cell_size                = 0.6,
      trajectory_graph_color   = "grey30",
      trajectory_graph_segment_size = 0.6
    ) +
      ggtitle("Pseudotime") +
      theme(plot.title      = element_text(size = 9),
            legend.position = "bottom",
            legend.text     = element_text(size = 7))

    # Module score overlay plots for the top four relevant modules
    overlays_present <- intersect(overlays, colnames(colData(cds)))
    mod_plots <- lapply(head(overlays_present, 4), function(ov) {
      plot_cells(
        cds,
        color_cells_by        = ov,
        show_trajectory_graph = FALSE,
        label_cell_groups     = FALSE,
        cell_size             = 0.5
      ) +
        scale_colour_gradient2(low = "blue", mid = "white",
                                high = "red", midpoint = 0) +
        ggtitle(gsub("_", " ", ov)) +
        theme(plot.title      = element_text(size = 8),
              legend.position = "bottom",
              legend.text     = element_text(size = 6))
    })

    top_row <- p_tp + p_pt
    patient_panel <- if (length(mod_plots) >= 2) {
      top_row / wrap_plots(mod_plots, nrow = 1) +
        plot_layout(heights = c(1.2, 1))
    } else {
      top_row
    }

    pt_plots_by_patient[[pt]] <- patient_panel

    rm(sub_obj, cds); gc()
  }

  # Save individual patient trajectory figures
  for (pt in names(pt_plots_by_patient)) {
    cairo_pdf(
      file.path(path.expand(fig_dir),
                sprintf("%s_%s_trajectory.pdf", ct_name, pt)),
      width = 10, height = 8
    )
    print(pt_plots_by_patient[[pt]])
    dev.off()
  }

  # Save consolidated figure across all patients for this cell type
  if (length(pt_plots_by_patient) > 0) {
    n_pts  <- length(pt_plots_by_patient)
    n_cols <- min(4, n_pts)

    consolidated <- wrap_plots(pt_plots_by_patient,
                                ncol = n_cols) +
      plot_annotation(
        title    = sprintf("%s Trajectory Across All Patients",
                           gsub("_", " ", ct_name)),
        subtitle = "Ordered by pseudotime from T1 root",
        theme    = theme(
          plot.title    = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 10, colour = "grey50")
        )
      )

    cairo_pdf(
      file.path(path.expand(fig_dir),
                sprintf("%s_ALL_PATIENTS_trajectory.pdf",
                        ct_name)),
      width  = n_cols * 10,
      height = ceiling(n_pts / n_cols) * 8
    )
    print(consolidated)
    dev.off()
  }

  all_pt_genes[[ct_name]] <- ct_pt_genes

  # Within cell type consensus: genes with consistent pseudotime
  # direction across patients, with emphasis on stable vs progressor
  # directionality differences
  if (length(ct_pt_genes) >= 2) {

    all_genes_ct <- bind_rows(ct_pt_genes)

    consensus <- all_genes_ct %>%
      group_by(gene) %>%
      summarise(
        n_patients   = n(),
        n_positive   = sum(rho > 0),
        n_negative   = sum(rho < 0),
        mean_rho     = round(mean(rho), 3),
        n_prog       = sum(patient %in% progressor_pts),
        n_prog_pos   = sum(patient %in% progressor_pts & rho > 0),
        n_prog_neg   = sum(patient %in% progressor_pts & rho < 0),
        n_stable     = sum(patient %in% stable_pts),
        n_stable_pos = sum(patient %in% stable_pts & rho > 0),
        n_stable_neg = sum(patient %in% stable_pts & rho < 0),
        patients_sig = paste(sort(patient), collapse = ","),
        .groups      = "drop"
      ) %>%
      mutate(
        direction = ifelse(n_positive >= n_negative,
                            "UP along pseudotime",
                            "DOWN along pseudotime"),
        opposite  = (n_stable > 0) & (
          (n_prog_pos >= n_prog_neg &
             n_stable_neg > n_stable_pos) |
            (n_prog_neg >= n_prog_pos &
               n_stable_pos > n_stable_neg)
        )
      ) %>%
      filter(n_patients >= 2) %>%
      arrange(desc(pmax(n_positive, n_negative)),
              desc(abs(mean_rho)))

    fwrite(consensus,
           file.path(path.expand(out_dir),
                     sprintf("%s_pseudotime_consensus.csv",
                             ct_name)))

    strong <- consensus %>%
      filter(pmax(n_positive, n_negative) >= 3)

    if (nrow(strong) > 0) {
      cat(sprintf("  Strong consensus genes (3+ patients):\n"))
      print(head(strong[, c("gene","n_patients",
                             "mean_rho","direction")], 10))
    }
  }

  cat(ct_name, "complete\n\n")
}

# Cross-cell-type consensus: genes showing consistent pseudotime
# correlation direction across multiple cell types in progressors
cross_ct <- bind_rows(lapply(names(all_pt_genes), function(ct) {
  genes_ct <- all_pt_genes[[ct]]
  if (length(genes_ct) == 0) return(NULL)
  df           <- bind_rows(genes_ct)
  df$cell_type <- ct
  df
}))

if (!is.null(cross_ct) && nrow(cross_ct) > 0) {

  cross_summary <- cross_ct %>%
    filter(patient %in% progressor_pts) %>%
    group_by(gene, direction) %>%
    summarise(
      n_cell_types = n_distinct(cell_type),
      n_patients   = n_distinct(patient),
      cell_types   = paste(sort(unique(cell_type)),
                            collapse = ","),
      mean_rho     = round(mean(rho), 3),
      .groups      = "drop"
    ) %>%
    filter(n_cell_types >= 2) %>%
    arrange(desc(n_cell_types), desc(abs(mean_rho)))

  fwrite(cross_summary,
         file.path(path.expand(out_dir),
                   "cross_celltype_pseudotime_consensus.csv"))

  cat("Cross-cell-type consensus genes (2+ cell types):\n")
  print(head(cross_summary[, c("gene","n_cell_types",
                                "n_patients","direction",
                                "cell_types")], 20))
}

# Save all pseudotime metadata combined across patients and cell types
all_meta_combined <- rbindlist(
  lapply(all_cds_meta, as.data.table),
  fill = TRUE
)
fwrite(all_meta_combined,
       file.path(path.expand(out_dir),
                 "all_pseudotime_metadata.csv"))

cat("\nTrajectory analysis complete\n")
cat(format(Sys.time()), "\n")
cat("Figures saved to:", path.expand(fig_dir), "\n")
cat("Data saved to:",    path.expand(out_dir), "\n")