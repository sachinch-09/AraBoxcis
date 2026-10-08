# AraBoxcis: Transcriptional Architecture of Photomorphogenesis in Seedlings

Research project, University of York, 2025 to 2026.

Single-cell RNA-seq analysis of *Arabidopsis* seedlings spanning pre-illumination developmental stages, aimed at inferring the gene regulatory network (GRN) controlling the transition into photomorphogenesis (light-driven development).

## Key result

Identified **PIF5 (PIL6)** as the dominant regulatory hub in the inferred network, alongside interconnected repressors **PIF4, HY5, and ABF3**, via cross-pathway heatmap analysis. This is consistent with PIF5's known role as a central integrator of light and developmental signalling in *Arabidopsis*.

## Pipeline

1. **Data:** single-cell RNA-seq transcripts from the Seedling_3D dataset, spanning pre-illumination developmental stages.
2. **QC and filtering (R):** a 1% minimum cell-expression threshold applied to remove transcript dropout artifacts.
3. **Dimensionality reduction:** PCA-initialized UMAP embeddings generated to visualise cell state structure.
4. **GRN inference:** directed gene regulatory network inferred via tree-based regression (GENIE3).
5. **Hub identification:** regulatory hub genes identified via cross-pathway heatmap analysis of the inferred network.

## Repository structure

```
AraBoxcis/
├── Figures/    UMAP embeddings, GRN heatmaps, regulatory network plots
├── dev/        analysis scripts (QC, UMAP, GENIE3 GRN inference)
└── README.md
```

*(Update this section to match your actual file names inside `dev/`, listing each script and which stage of the pipeline it covers.)*

## Requirements

- R 
- Key packages: Seurat (or equivalent single-cell toolkit), GENIE3, umap, igraph, ggraph

*(Add a package list or `renv.lock` if you want this reproducible by someone else.)*

## Attribution

Research project conducted as part of the MSc Bioinformatics programme, University of York. This repository contains my individual analysis and code.
