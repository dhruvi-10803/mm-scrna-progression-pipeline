# 02_module_scoring.R
# Pipeline step 2: Functional state module scoring per cell type.
# Module scores are computed within each broad cell type compartment
# rather than across the full object, ensuring scores reflect
# meaningful variation within the relevant cell population.
# Input:  Per-patient annotated objects from objects/Px/Px_after_azimuth.rds
# Output: Per-patient scored objects in consolidation/module_scoring/Px_scored.rds

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

obj_dir <- "~/Dhruvi/TME/objects/"
out_dir <- "~/Dhruvi/consolidation/module_scoring/"
dir.create(path.expand(out_dir), recursive = TRUE,
           showWarnings = FALSE)

patients <- paste0("P", 1:7)

# Module gene sets curated from published literature.
# Each module is scored within its relevant broad cell type compartment
# to ensure the score reflects variation within that population.
# Sources: Azizi et al. 2018, Zheng et al. 2021, Gayoso et al. 2022,
# van Galen et al. 2019, Tirosh et al. 2016

module_definitions <- list(

  # Monocyte states — scored within Mono compartment
  Mono = list(

    Mono_Homeostatic = list(
      genes = c("CX3CR1","FCGR3A","LST1","MS4A7",
                "LILRB1","LILRB2","LRRC25","TCF7L2",
                "HES4","IFITM3","CDKN1C","CSF1R"),
      compartment = "Mono"
    ),

    Mono_Classical_Inflammatory = list(
      genes = c("CD14","LYZ","S100A8","S100A9","VCAN",
                "LGALS3","IL1B","CXCL8","CCL2","CCL7",
                "THBS1","CD68","PLAUR","G0S2","CLEC4E"),
      compartment = "Mono"
    ),

    Mono_Immunosuppressive = list(
      genes = c("CD274","IDO1","HMOX1","IL10","TGFB1",
                "LILRB2","LILRB4","VSIG4","MRC1","CD163",
                "MARCO","MSR1","CCL13","CCL18","FOLR2"),
      compartment = "Mono"
    ),

    Mono_IFN_Activated = list(
      genes = c("ISG15","ISG20","IFI6","IFI27","IFI44L",
                "IFIT1","IFIT2","IFIT3","IFITM1","MX1",
                "MX2","OAS1","OAS2","OASL","IRF7",
                "RSAD2","HERC5","CXCL10","CXCL9","STAT1"),
      compartment = "Mono"
    ),

    Mono_AntigenPresenting = list(
      genes = c("HLA-DRA","HLA-DRB1","HLA-DPB1","HLA-DPA1",
                "HLA-DQA1","HLA-DQB1","CD74","CIITA",
                "CD86","CD80","CD83","CCR7","LAMP3",
                "FSCN1","MARCKSL1","IL12B"),
      compartment = "Mono"
    )
  ),

  # CD8 T cell states — scored within CD8 T compartment
  CD8 = list(

    CD8_NaiveStemlike = list(
      genes = c("TCF7","CCR7","SELL","LEF1","KLF2",
                "IL7R","S1PR1","FOXO1","MYC","ID3",
                "BACH2","BCL6","SATB1","CD27","CD28"),
      compartment = "CD8 T"
    ),

    CD8_CytotoxicMemory = list(
      genes = c("GZMK","GZMA","CCL5","CCL4","CXCR3",
                "CXCR4","CD44","EOMES","TBX21","IFNG",
                "TNF","IL2","KLRG1","CX3CR1","S1PR5",
                "NKG7","CST7","LYAR","KLRD1","KLRB1"),
      compartment = "CD8 T"
    ),

    CD8_EffectorCytotoxic = list(
      genes = c("GZMB","PRF1","GNLY","FGFBP2","FCGR3A",
                "CX3CR1","KLRG1","S1PR5","TBX21","EOMES",
                "NKG7","CTSW","CD57","KLRD1","ADGRG1"),
      compartment = "CD8 T"
    ),

    CD8_Exhausted = list(
      genes = c("PDCD1","LAG3","TIGIT","HAVCR2","CTLA4",
                "TOX","TOX2","NR4A1","NR4A2","NR4A3",
                "ENTPD1","CXCL13","LAYN","PHLDA1",
                "BATF","IRF4","RBPJ","IKZF2","MYO7A"),
      compartment = "CD8 T"
    ),

    CD8_IFN_Activated = list(
      genes = c("ISG15","MX1","IFIT1","IFIT2","IFIT3",
                "IFI6","IFI27","IRF7","STAT1","OAS1",
                "RSAD2","HERC5","IFI44L","CXCL10","BST2"),
      compartment = "CD8 T"
    )
  ),

  # CD4 T cell states — scored within CD4 T compartment
  CD4 = list(

    CD4_Helper_Memory = list(
      genes = c("IL7R","CCR7","S1PR1","KLF2","SELL",
                "TCF7","LEF1","CD40LG","ICOS","IL2",
                "CXCR5","BCL6","ASCL2","MAF","BATF",
                "PDCD1","TIGIT","CTLA4","CD200","FCRL2"),
      compartment = "CD4 T"
    ),

    CD4_Regulatory = list(
      genes = c("FOXP3","IL2RA","CTLA4","TIGIT","IKZF2",
                "ENTPD1","LAYN","TNFRSF9","CCR8","BATF",
                "IRF4","PRDM1","IL10","TGFB1","TNFRSF18",
                "DUSP4","MAGEH1","RTKN2","FCRL3","LRRC32"),
      compartment = "CD4 T"
    ),

    CD4_Exhausted = list(
      genes = c("PDCD1","LAG3","TIGIT","HAVCR2","TOX",
                "NR4A1","CXCL13","ENTPD1","CTLA4",
                "BATF","IRF4","RBPJ","PHLDA1","LAYN"),
      compartment = "CD4 T"
    ),

    CD4_IFN_Activated = list(
      genes = c("ISG15","MX1","IFIT1","IFIT2","IFIT3",
                "IFI6","IFI27","IRF7","STAT1","OAS1",
                "RSAD2","HERC5","IFI44L","BST2","CXCL10"),
      compartment = "CD4 T"
    )
  ),

  # NK cell states — scored within NK compartment
  NK = list(

    NK_Cytotoxic_Mature = list(
      genes = c("GNLY","GZMB","PRF1","NKG7","FGFBP2",
                "FCGR3A","CX3CR1","KLRD1","KLRF1","KLRB1",
                "S1PR5","TBX21","EOMES","TYROBP","NCAM1",
                "NCR1","NCR3","KIR2DL1","KIR2DL3","IFNG"),
      compartment = "NK"
    ),

    NK_Immature_Tissue = list(
      genes = c("CD56","NCAM1","XCL1","XCL2","SELL",
                "CCR7","IL7R","KIT","GATA3","TCF7",
                "LEF1","CXCR4","CD117","SPINK2","ALDH1A1"),
      compartment = "NK"
    ),

    NK_Dysfunctional = list(
      genes = c("LAG3","TIGIT","PDCD1","HAVCR2","CTLA4",
                "TOX","NR4A1","ENTPD1","CXCL13","BATF",
                "CD96","KIR2DL1","KIR3DL1","LILRB1",
                "SIGLEC7","SIGLEC9"),
      compartment = "NK"
    ),

    NK_IFN_Activated = list(
      genes = c("ISG15","MX1","IFIT1","IFIT2","IFIT3",
                "IFI6","IFI27","IRF7","STAT1","OAS1",
                "RSAD2","CXCL10","BST2","HERC5","IFI44L"),
      compartment = "NK"
    )
  ),

  # Erythroid states — scored within Erythroid compartment
  Erythroid = list(

    Erythroid_Progenitor = list(
      genes = c("CD34","KIT","GATA2","TAL1","FLI1",
                "RUNX1","LMO2","MEIS1","HOXA9","HOXA10",
                "CD117","CD38","CD133","PROM1","MYCN"),
      compartment = "Erythroid"
    ),

    Erythroid_Mature = list(
      genes = c("HBB","HBA1","HBA2","ALAS2","AHSP",
                "GYPB","GYPA","SLC4A1","EKLF","KLF1",
                "EPOR","TFRC","SPTA1","ANK1","EPB42"),
      compartment = "Erythroid"
    ),

    Erythroid_Stress_IFN = list(
      genes = c("IFI27","ISG15","IFI6","MX1","IFIT1",
                "BST2","RSAD2","OAS1","IRF7","STAT1",
                "IFITM1","IFITM2","IFITM3","CXCL10","ISG20"),
      compartment = "Erythroid"
    )
  ),

  # HSPC states — scored within HSPC compartment
  HSPC = list(

    HSPC_StemQuiescent = list(
      genes = c("CD34","MLLT3","MECOM","SPINK2","HOPX",
                "HLF","MEIS1","GATA2","AVP","CRHBP",
                "PROCR","PROM1","KIT","FLT3","ALDH1A1"),
      compartment = "HSPC"
    ),

    HSPC_GMP_Myeloid = list(
      genes = c("MPO","ELANE","AZU1","PRTN3","CTSG",
                "CSF3R","CEBPA","CEBPE","GFI1","SPI1",
                "FCGR3B","CD66B","ITGAM","CXCR2","G0S2"),
      compartment = "HSPC"
    ),

    HSPC_Lymphoid_CLP = list(
      genes = c("IL7R","RAG1","RAG2","DNTT","CD79A",
                "VPREB1","IGLL1","EBF1","PAX5","IKZF1",
                "TCF3","ID2","FLT3","CD19","CD10"),
      compartment = "HSPC"
    ),

    HSPC_Megakaryocyte = list(
      genes = c("PF4","PPBP","GP9","ITGA2B","TUBB1",
                "VWF","SELP","MYH9","GATA1","FOG1",
                "TMEM40","CLU","GNG11","TREML1","SPARC"),
      compartment = "HSPC"
    )
  ),

  # Plasma cell states — scored within Plasma compartment
  Plasma = list(

    Plasma_Normal = list(
      genes = c("MZB1","JCHAIN","DERL3","XBP1","SEC11C",
                "PDIA4","PDIA6","HSPA5","HSP90B1","PPIB",
                "SSR3","FKBP11","IGKC","IGLC2","IGHG1"),
      compartment = "Plasma"
    ),

    Plasma_Myeloma_Proliferating = list(
      genes = c("MKI67","TOP2A","PCNA","MCM2","MCM6",
                "CCNB1","CCNB2","CDK1","BUB1","BUB1B",
                "AURKB","PLK1","TYMS","RRM2","CENPF",
                "UBE2C","HJURP","NUSAP1","TPX2","KIF20A"),
      compartment = "Plasma"
    ),

    Plasma_Myeloma_Survival = list(
      genes = c("MCL1","BCL2","CCND1","CCND2","IRF4",
                "MYC","SDC1","TNFRSF17","TNFRSF13B",
                "TNFRSF13C","IL6R","CXCR4","VLA4",
                "CD44","ITGA4","FUT4","CD200","HM13",
                "PRDM1","ELK4"),
      compartment = "Plasma"
    ),

    Plasma_IFN_Stress = list(
      genes = c("ISG15","IFI6","IFI27","MX1","IFIT1",
                "IFIT3","RSAD2","OAS1","IRF7","BST2",
                "HERC5","IFI44L","CXCL10","STAT1","ISG20"),
      compartment = "Plasma"
    )
  )
)

# Map broad compartment names to Azimuth L1 labels
compartment_to_az_l1 <- list(
  "Mono"      = c("CD14+ Monocytes","CD16+ Monocytes",
                   "Monocyte","Mono"),
  "CD8 T"     = c("CD8 T","CD8T","CD8+ T"),
  "CD4 T"     = c("CD4 T","CD4T","CD4+ T"),
  "NK"        = c("NK","NK cell","Natural Killer"),
  "Erythroid" = c("Erythroid","Early Eryth","Late Eryth",
                   "Erythrocyte"),
  "HSPC"      = c("HSPC","HSC","Progenitor","GMP",
                   "CLP","LMPP","EMP","Prog Mk"),
  "Plasma"    = c("Plasma","Plasmablast","Plasma cell",
                   "PC","Plasma_cell")
)

for (pt in patients) {

  cat("Patient:", pt, "\n")

  obj_path <- file.path(
    path.expand(obj_dir), pt,
    paste0(pt, "_after_azimuth.rds")
  )

  if (!file.exists(obj_path)) {
    cat("  Object not found, skipping\n")
    next
  }

  obj <- readRDS(obj_path)
  cat(sprintf("  Loaded: %d cells\n", ncol(obj)))

  # Initialise score columns as NA in full object
  all_module_names <- unlist(
    lapply(module_definitions, names)
  )
  for (mod in all_module_names) {
    obj@meta.data[[mod]] <- NA_real_
  }

  # Score each module within its relevant compartment
  for (compartment_name in names(module_definitions)) {

    compartment_modules <- module_definitions[[compartment_name]]
    az_labels <- compartment_to_az_l1[[compartment_name]]

    # Subset to the relevant broad cell type
    cells_keep <- rownames(obj@meta.data)[
      obj@meta.data$az_l1 %in% az_labels
    ]

    if (length(cells_keep) < 20) {
      cat(sprintf("  Skipping %s: only %d cells\n",
                  compartment_name, length(cells_keep)))
      next
    }

    cat(sprintf("  Scoring %s compartment: %d cells\n",
                compartment_name, length(cells_keep)))

    sub_obj <- subset(obj, cells = cells_keep)

    for (mod_name in names(compartment_modules)) {

      gene_set <- compartment_modules[[mod_name]]$genes

      # Keep only genes present in the object
      genes_use <- gene_set[gene_set %in% rownames(sub_obj)]

      if (length(genes_use) < 5) {
        cat(sprintf("    Skipping %s: only %d genes found\n",
                    mod_name, length(genes_use)))
        next
      }

      cat(sprintf("    %s: %d/%d genes\n",
                  mod_name, length(genes_use),
                  length(gene_set)))

      # Compute module score within the compartment subset
      sub_obj <- AddModuleScore(
        sub_obj,
        features  = list(genes_use),
        name      = mod_name,
        assay     = "RNA",
        seed      = 42
      )

      # Score column gets a trailing 1 appended by AddModuleScore
      score_col <- paste0(mod_name, "1")

      # Transfer scores back to full object metadata
      obj@meta.data[cells_keep, mod_name] <-
        sub_obj@meta.data[cells_keep, score_col]

      rm(sub_obj)
    }
  }

  # Save scored object
  out_path <- file.path(
    path.expand(out_dir),
    paste0(pt, "_scored.rds")
  )
  saveRDS(obj, out_path)

  # Save metadata with scores as CSV for easy access downstream
  meta_out <- obj@meta.data[, c(
    "patient","timepoint","stage","group",
    "az_l1","az_l2","mapping.score",
    all_module_names
  ), drop = FALSE]

  write.csv(
    meta_out,
    file.path(path.expand(out_dir),
              paste0(pt, "_module_scores.csv")),
    row.names = TRUE
  )

  cat(sprintf("  Saved: %s\n", out_path))

  rm(obj)
  gc()
}

# Combine metadata across all patients into one file
cat("\nCombining metadata across all patients\n")

all_meta <- lapply(patients, function(pt) {
  f <- file.path(path.expand(out_dir),
                  paste0(pt, "_module_scores.csv"))
  if (file.exists(f)) read.csv(f, row.names = 1)
  else NULL
})

all_meta <- do.call(rbind, Filter(Negate(is.null), all_meta))

write.csv(
  all_meta,
  file.path(path.expand(out_dir),
            "all_patients_metadata.csv"),
  row.names = TRUE
)

cat(sprintf("Combined metadata: %d cells across all patients\n",
            nrow(all_meta)))
cat("Module scoring complete\n")