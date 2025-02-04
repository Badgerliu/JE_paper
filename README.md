# p63-RUNX1-ODAM Axis in Junctional Epithelium Repair

This repository contains code and analysis scripts for the paper "p63-RUNX1-ODAM Axis is Essential for Junctional Epithelium Repair". The code is organized by figure and analysis type.

## Repository Structure

### ATACSeq_DiiffBind/
- `Gingiva_palate_ATAC_240506.Rmd`: R Markdown for ATAC-seq differential binding analysis
- `Gingiva_palate_ATAC_samplelist.csv`: Sample metadata for ATAC-seq analysis

### ATACSeq_mapping_code/
- `ATAC_oneClick_mm10_xiyou_v1.1.sh`: Shell script for ATAC-seq data processing
- `gingiva_palate_LY_ATAC_0503_output.txt`: Example output from ATAC-seq mapping

### Figure_1_related/
- `Figure_1_plots.Rmd`: R Markdown for generating Figure 1 plots

### Figure_2_related/
- `Figure_2_plots.Rmd`: R Markdown for Figure 2 plots
- `Integrate_JE_epi_harmony_imputation_3D.Rmd`: Integration and imputation of JE epithelial cells
- `Integrate_JE_epi_harmony_imputation_pseudotime_slingshot_2.Rmd`: Pseudotime analysis using Slingshot
- `scVelo_integrated_JE_epi_dynamic_model.ipynb`: Jupyter notebook for scVelo dynamic modeling

### Figure_4_related/
Contains scripts for spatial transcriptomics and trajectory analyses:
- Spatial analysis scripts: `3.9_spatial.R`, `3.14_spatial_Bayes.R`, `4.5_space_bayes_prepare.R`, `4.6_spatialDE.ipynb`
- Trajectory analysis: `4.15_TrajAtlas_GEP.ipynb`
- Gene expression programs: `6.14_aucell.ipynb`, `6.18_wound_go.R`
- Figure generation: `9.23_figure4_revision.R`, `Figure4B.R`, `Figure4C.R`
- Shiny app preparation: `shinyPrepare.R`

### Figure_5_related/
- `Figure_5_scRNA_plots.Rmd`: R Markdown for generating Figure 5 single-cell RNA-seq plots

## Requirements

The code in this repository requires:
- R (≥ 4.0)
- Python (≥ 3.8)
- Jupyter Notebook
- Common bioinformatics packages (Seurat, Scanpy, etc.)

## Usage

1. Clone this repository
2. Install required dependencies
3. Run individual scripts for specific analyses

## Citation

If you use this code in your research, please cite our paper:

[Citation information will be added here once published]

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
