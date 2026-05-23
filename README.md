# mm-scrna-progression-pipeline
Longitudinal single-cell RNA-seq analysis pipeline for multiple myeloma progression
# Longitudinal Single-Cell RNA-seq Pipeline for Multiple Myeloma Progression

## Overview

This repository contains a complete analytical pipeline for longitudinal 
single-cell RNA sequencing of bone marrow samples across the multiple myeloma 
disease continuum. The pipeline was developed for a study profiling seven 
patients with paired samples at two timepoints spanning the MGUS to MM 
progression spectrum.

The pipeline covers all major analytical steps from raw data processing 
through to cell-cell communication inference and clonal evolution analysis, 
and is designed to be adaptable to other longitudinal single-cell studies.

## Dataset

- **GEO accession:** GSE271107
- **Patients:** 7 (2 stable MGUS, 2 MGUS progressors, 3 SMM-MM progressors)
- **Design:** Paired bone marrow aspirates at T1 (baseline) and T2 (follow-up)
- **Sequencing:** 10x Genomics Chromium single-cell RNA-seq

## Key Analytical Steps

**QC and annotation** — Per-patient quality control using adaptive thresholds, 
doublet removal with scDblFinder, normalization with SCTransform, and cell type 
annotation using the Azimuth human bone marrow reference atlas.

**Module scoring** — Functional immune state scores computed within each broad 
cell type compartment using curated gene sets from published literature covering 
monocyte, T cell, NK cell, B cell, plasma cell, erythroid, and HSPC states.

**Longitudinal DEG** — Within-patient differential expression between T1 and T2 
using the MAST hurdle model with cellular detection rate as covariate, run 
independently per patient per cell type.

**Pathway enrichment** — fgsea on MAST results using Hallmark, KEGG, and 
Reactome gene sets, with cross-patient directionality analysis to identify 
pathways consistently changed in progressors versus stable patients.

**Cell-cell communication** — CellChat inference of ligand-receptor 
communication probabilities at T1 and T2 separately, with pathway-level 
log2FC comparison between timepoints.

**Clonal evolution** — inferCNV-based copy number inference in plasma cells 
using patient-matched immune cells and healthy donor cells as diploid reference, 
followed by subclone tracking and evolutionary mode classification.

**Trajectory analysis** — Monocle3 pseudotime analysis per cell type, ordering 
cells from T1 as root to identify genes that change along the progression 
trajectory within individual patients.

## Requirements

**R version:** 4.3.3

**Key packages:**
- Seurat 5.0
- CellChat 1.6.1
- inferCNV
- MAST
- fgsea
- monocle3
- scDblFinder
- Azimuth
- ComplexHeatmap
- harmony

**Environment setup:**
```bash
conda env create -f envs/mm_scrna.yml
conda activate mm_scRNA
```

## Usage

Scripts are numbered and should be run in order. Each script reads from 
outputs of the previous step. Key input and output paths are defined at 
the top of each script and can be modified for a new dataset.

```bash
Rscript pipeline/01_QC_normalization_annotation.R
Rscript pipeline/02_module_scoring.R
Rscript pipeline/03_longitudinal_DEG.R
Rscript pipeline/04_fgsea_longitudinal.R
Rscript pipeline/05_cellchat.R
Rscript pipeline/06_infercnv_clonal_evolution.R
Rscript pipeline/07_trajectory.R
```

Scripts 03 to 07 are computationally intensive and are designed to be run 
with nohup or submitted to a cluster. Approximate runtimes on 8 cores 
with 64GB RAM per patient: DEG 2-4 hours, CellChat 1-2 hours, 
inferCNV 6-12 hours.
